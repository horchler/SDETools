function [Y W] = sde_euler(f,g,tspan,y0,options,varargin)
%SDE_EULER  Solve stochastic differential equations, Euler methods.
%   YOUT = SDE_EULER(FFUN,GFUN,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL]
%   integrates the system of stochastic differential equations dy = f(t,y)*dt +
%   g(t,y)*dW with diagonal noise by default from time T0 to TFINAL (all
%   increasing or all decreasing with arbitrary step size) with initial conditions
%   Y0. FFUN and GFUN are function handles. For a scalar T and a vector Y,
%   FFUN(T,Y) and GFUN(T,Y) return column vectors corresponding to f(t,y) and
%   g(t,y), the deterministic and stochastic parts of the SDE, respectively.
%   Each row in the solution array YOUT corresponds to a time in the input
%   vector TSPAN.
%
%   YOUT = SDE_EULER(FFUN,GFUN,TSPAN,Y0,OPTIONS) solves the above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the SDESET function. See SDESET for details. The type of SDE to be
%   integrated, 'Ito' or the default 'Stratonovich', can be specified via the
%   SDEType property. GFUN can output more general D-dimensional correlated
%   noise by setting the DiagonalNoise property to 'no'. The Euler-Maruyama
%   method is used for Ito SDEs. The Euler-Heun method is used for Stratonovitch
%   SDEs unless either the AdditiveNoise or the ConstGFUN property is set to
%   'yes', in which case the Ito and Stratonovitch SDEs are equivalent. Another 
%   commonly used option is manually specifying the random seed with the
%   RandSeed property and creates a separate random number stream rather than
%   using the default one.
%
%   [YOUT, W] = SDE_EULER(FFUN,GFUN,TSPAN,Y0,...) outputs the matrix W of
%   integrated Weiner increments that were used by the solver. W is LENGTH(Y0)
%   rows by LENGTH(TSPAN) columns, corresponding to [T0 T1 T2 ... TFINAL].
%
%   YOUT = SDE_EULER([],GFUN,TSPAN,Y0,...)
%
%   Example:
%       y = sde_euler(@f1,@g1,0:0.01:1,[1 0]);
%       solves the 2-D Stratonovich SDE system dy = f1(t,y)*dt + g1(t,y)*dW,
%       using the Euler-Heun method and the default random number stream.
%
%   NOTE:
%       SDEs are assumed to be in Stratonovich form by default. SDEs in Ito form
%       can be handled by setting the 'SDEType' OPTIONS property to 'Ito' via
%       the SDESET function. These forms are generally not equivalent and will
%       converge to different solutions, so care should be taken to ensure that
%       the form of SDEFUN matches 'SDEType'.
%
%   See also:
%       Explicit SDE solvers:	SDE_MILSTEIN
%       Implicit SDE solvers:   
%       Option handling:        SDESET, SDEGET
%       SDE demos:           
%   	Other:                  FUNCTION_HANDLE, RANDSTREAM

%   SDE_EULER is an implementation of the order 0.5 strong (order 1.0 weak)
%   explicit Euler-Maruyama and Euler-Heun schemes, which are also order 0.5
%   strong Taylor approximation schemes. In the additive noise case, Euler
%   schemes improve to order 1.0 strong convergence.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, Created 10-28-10
%   Revision: 1.0, 1-2-12


solver = 'SDE_EULER';

% Check inputs and outputs
if nargin < 5
    if nargin < 4
        error(  'SDELab:sde_euler:NotEnoughInputs',...
                'Not enough input arguments.  See %s.',solver);
    end
    if isa(y0,'struct')
        error(  'SDELab:sde_euler:NotEnoughInputsOptions',...
               ['An SDE options structure was provided as the last '...
                'argument,  but one of the first four input arguments is '...
                'missing.  See %s.'],solver);
    end
    options = [];
elseif isempty(options) && (ndims(options) ~= 2 || any(size(options) ~= 0) ...
        || ~isnumeric(options) || ~iscell(options) || ~isstruct(options)) ...
        || (~isempty(options) && ~isa(options,'struct'))
	error(  'SDELab:sde_euler:InvalidSDESETStruct',...
            'Invalid SDE options structure.  See SDESET.');
end

% Handle solver arguments
[N D tspan tdir lt y0 fout gout h ConstStep dataType DiagonalNoise ...
	ScalarNoise ConstFFUN ConstGFUN Stratonovich RandFUN CustomRandFUN]...
	= sdearguments(solver,f,g,tspan,y0,options,varargin);

Y = zeros(lt,N,dataType);   % State array

% Calculate Wiener increments from normal variates, store in state if possible
sh = tdir*sqrt(h);
h = tdir*h;
if CustomRandFUN    % check output of alternative RandFUN
    try
        if D > N
            if nargout == 2 % Store Wiener increments in W
                W = feval(RandFUN,D,lt);
                if ndims(W) ~= 2 || isempty(W) || ~isfloat(W)
                    error(  'SDELab:sde_euler:RandFUNNot2DArray1',...
                           ['RandFUN must return a non-empty matrix of '...
                            'floating point values.  See %s.'],solver);
                end
                [m n] = size(W);
                if m ~= lt || n ~= D
                    error(  'SDELab:sde_euler:RandFUNDimensionMismatch1',...
                           ['The specified alternative RandFUN did not '...
                            'output a %d by %d matrix as requested.'...
                            '  See %s.',D,lt,solver]);
                end
                if ConstStep
                    W = sh*W;
                    W(1,1:D) = zeros(1,D,dataType);
                else
                    W = bsxfun(@times,[0;sh],W);
                end
            else            % Unable to store Wiener increments
                dW = feval(RandFUN,1,D);
                if ndims(dW) ~= 2 || isempty(dW) || ~isfloat(dW)
                    error(  'SDELab:sde_euler:RandFUNNot2DArray2',...
                           ['RandFUN must return a non-empty matrix of '...
                            'floating point values.  See %s.'],solver);
                end
                [m n] = size(dW);
                if m ~= 1 || n ~= D
                    error(  'SDELab:sde_euler:RandFUNDimensionMismatch2',...
                           ['The specified alternative RandFUN did not '...
                            'output a %d by 1 column vector as requested.'...
                            '  See %s.',D,solver]);
                end
                dW = sh(1)*dW;
            end
        elseif D == N       % Store Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,D);
            if ndims(r) ~= 2 || isempty(r) || ~isfloat(r)
                error(  'SDELab:sde_euler:RandFUNNot2DArray3',...
                       ['RandFUN must return a non-empty matrix of floating '...
                        'point values.  See %s.'],solver);
            end
            [m n] = size(r);
            if m ~= lt-1 || n ~= D
                error(  'SDELab:sde_euler:RandFUNDimensionMismatch3',...
                       ['The specified alternative RandFUN did not output a '...
                        '%d by %d matrix as requested.'...
                        '   See %s.',D,lt-1,solver]);
            end
            if ConstStep
                Y(2:end,1:D) = sh*r;
            else
                Y(2:end,1:D) = bsxfun(@times,sh,r);
            end
            clear r;    % remove large temporary variable to save memory
        end
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error(  'SDELab:sde_euler:RandFUNTooFewInputs',...
                       ['RandFUN must have at least two inputs.'...
                        '  See %s.'],solver);
            case 'MATLAB:TooManyOutputs'
                error(  'SDELab:sde_euler:RandFUNNoOutput',...
                       ['The output of RandFUN was not specified. RandFUN '...
                        'must return a non-empty matrix.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error(  'SDELab:sde_euler:RandFUNUnassignedOutput',...
                       ['The first output of RandFUN was not assigned.'...
                        '  See %s.'],solver);
            case 'MATLAB:minrhs'
                error(  'SDELab:sde_euler:RandFUNTooManyInputs',...
                       ['RandFUN must not require more than two inputs.'...
                        '  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
else    % No error checking needed if default RANDN used
    if D > N
        if nargout == 2 % Store Wiener increments in W
            if ConstStep
                W = sh*feval(RandFUN,lt,D);
                W(1,1:D) = zeros(1,D,dataType);
            else
                W = bsxfun(@times,[0;sh],feval(RandFUN,lt,D));
            end
        else            % Unable to store Wiener increments
            dW = sh(1)*feval(RandFUN,1,D);
        end
    else                % Store Wiener increments in Y
        if ConstStep
            Y(2:end,1:D) = sh*feval(RandFUN,lt-1,D);
        else
            Y(2:end,1:D) = bsxfun(@times,sh,feval(RandFUN,lt-1,D));
        end
    end
end

% Only allocate W matrix if requested as output
if nargout == 2 && D <= N
	W = cumsum(Y(:,1:D),1);
end

% Integrate
if ConstFFUN && ConstGFUN && (D <= N || nargout == 2)   % no FOR loop needed
    if ScalarNoise
        if ConstStep
            Y(2:end,:) = bsxfun(@plus,y0',cumsum(bsxfun(@plus,fout'*h,Y(2:end,ones(1,N))*gout),1));
        else
            Y(2:end,:) = bsxfun(@plus,y0',cumsum(h*fout'+Y(2:end,ones(1,N))*gout,1));
        end
    elseif DiagonalNoise
        if ConstStep
            Y(2:end,:) = bsxfun(@plus,y0',cumsum(bsxfun(@plus,h*fout',bsxfun(@times,gout',Y(2:end,:))),1));
        else
            Y(2:end,:) = bsxfun(@plus,y0',cumsum(h*fout'+bsxfun(@times,gout',Y(2:end,:)),1));
        end
    else
        if ConstStep
            if D > N
                Y(2:end,:) = bsxfun(@plus,y0',cumsum(bsxfun(@plus,fout'*h,W(2:end,:)*gout'),1));
            else
                Y(2:end,:) = bsxfun(@plus,y0',cumsum(bsxfun(@plus,fout'*h,Y(2:end,1:D)*gout'),1));
            end
        else
            if D > N
                Y(2:end,:) = bsxfun(@plus,y0',cumsum(h*fout'+W(2:end,:)*gout',1));
            else
                Y(2:end,:) = bsxfun(@plus,y0',cumsum(bsxfun(@plus,h*fout',Y(2:end,1:D)*gout'),1));
            end
        end
    end
else
    % step size and Wiener increment for first time step
    dt=h(1);
    if D > N
        if nargout == 2      	% dW already exists if D > N && nargin == 1
            dW = W(2,:);
            Wi = W(1,:)+dW;
            W(2,:) = Wi;        % integrate Wiener increments
        end
    else
        dW = Y(2,1:D);
    end
    
    Y(1,:) = y0;	% Set initial conditions
    
    % Use existing fout and gout to calculate the first time step, Y(2,:)
    if ConstGFUN
        if D <= N        	% All values of gout*dW can be stored in Y array
            if ScalarNoise
                Y(2:end,:) = Y(2:end,ones(1,N))*gout;
            elseif DiagonalNoise
                Y(2:end,:) = bsxfun(@times,gout',Y(2:end,:));
            else
                Y(2:end,:) = Y(2:end,1:D)*gout';
            end
            Yi = y0+fout*dt+Y(2,:)';
        elseif nargout == 2	% All values of gout*dW can be stored in Y array
            if ScalarNoise
                Y(2:end,:) = W(2:end,ones(1,N))*gout;
            elseif DiagonalNoise
                Y(2:end,:) = bsxfun(@times,gout',W(2:end,:));
            else
                Y(2:end,:) = W(2:end,:)*gout';
            end
            Yi = y0+fout*dt+Y(2,:)';
        else               	% Wiener increments were not precalculated
            if DiagonalNoise
                Yi = y0+fout*dt+gout.*dW(:);
            else
                Yi = y0+fout*dt+gout*dW(:);
            end
        end
    else
        if Stratonovich	% Use Euler-Heun step
            if DiagonalNoise
                gout = 0.5*(gout+feval(g,tspan(1),y0+gout.*dW(:),varargin{:}));
                Yi = y0+fout*dt+gout.*dW(:);
            else
                gout = 0.5*(gout+feval(g,tspan(1),y0+gout*dW(:),varargin{:}));
                Yi = y0+fout*dt+gout*dW(:);
            end
        else
            if DiagonalNoise
                Yi = y0+fout*dt+gout.*dW(:);
            else
                Yi = y0+fout*dt+gout*dW(:);
            end
        end
    end
    Y(2,:) = Yi;
    
    % Integration loop using cached state, Yi, and increment, Wi
    for i=2:lt-1
        % Step size and Wiener increment
        if ConstStep
            if D > N
                if nargout == 2
                    dW = W(i+1,:);
                    Wi = Wi+dW;
                    W(i+1,:) = Wi;              % Integrate Wiener increments
                else
                    dW = sh*feval(RandFUN,1,D);     % Generate Wiener increments
                end
            else
                dW = Y(i+1,1:D);
            end
        else    % Variable time step case
            dt=h(i);
            if D > N
                if nargout == 2
                    dW = W(i+1,:);
                    Wi = Wi+dW;
                    W(i+1,:) = Wi;              % Integrate Wiener increments
                else
                    dW = sh(i)*feval(RandFUN,1,D);  % Generate Wiener increments
                end
            else
                dW = Y(i+1,1:D);
            end
        end
        
        % Calculate next time step
        if ConstGFUN
            if D > N && nargout == 1
                if DiagonalNoise
                    Yi = Yi+feval(f,tspan(i),Yi,varargin{:})*dt+gout.*dW(:);
                else
                    Yi = Yi+feval(f,tspan(i),Yi,varargin{:})*dt+gout*dW(:);
                end
            else
                Yi = Yi+feval(f,tspan(i),Yi,varargin{:})*dt+Y(i+1,:)';
            end
        else
            if ~ConstFFUN
                fout = feval(f,tspan(i),Yi,varargin{:})*dt;
            end
            
            if Stratonovich	% Use Euler-Heun step
                if DiagonalNoise
                    gout = 0.5*(gout+feval(g,tspan(i),Yi+gout.*dW(:),varargin{:}));
                    Yi = Yi+fout+gout.*dW(:);
                else
                    gout = 0.5*(gout+feval(g,tspan(i),Yi+gout*dW(:),varargin{:}));
                    Yi = Yi+fout+gout*dW(:);
                end
            else
                if DiagonalNoise
                    Yi = Yi+fout+feval(g,tspan(i),Yi,varargin{:}).*dW(:);
                else
                    Yi = Yi+fout+feval(g,tspan(i),Yi,varargin{:})*dW(:);
                end
            end
        end
        Y(i+1,:) = Yi;
    end
end