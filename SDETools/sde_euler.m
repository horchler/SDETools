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
%   [YOUT, W] = SDE_EULER(FFUN,GFUN,TSPAN,Y0,...) outputs the matrix W of
%   integrated Weiner increments that were used by the solver. W is LENGTH(Y0)
%   rows by LENGTH(TSPAN) columns, corresponding to [T0 T1 T2 ... TFINAL].
%
%   [...] = SDE_EULER(FFUN,GFUN,TSPAN,Y0,OPTIONS) solves the above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the SDESET function. See SDESET for details. The type of SDE to be
%   integrated, 'Ito' or the default 'Stratonovich', can be specified via the
%   SDEType property. GFUN can output more general D-dimensional correlated
%   noise by setting the DiagonalNoise property to 'no'. The Euler-Maruyama
%   method is used for Ito SDEs. The Euler-Heun method is used for Stratonovitch
%   SDEs unless either the AdditiveNoise or the ConstGFUN property is set to
%   'yes', in which case the Ito and Stratonovitch SDEs are equivalent. Another
%   commonly used option is to manually specify the random seed via the RandSeed
%   property, which creates a new random number stream, instead of using the
%   default stream, to generate the Wiener increments.
%
%   YOUT = SDE_EULER([],GFUN,TSPAN,Y0,...)
%
%   Example:
%       % Solve 2-D Stratonovich SDE using Euler-Heun method
%       mu = 1; sig = [0.1;0.5]; dt = 1e-2; t = 0:dt:1;
%       f = @(t,y)mu.*y; g = @(t,y)sig.*y;
%       y = sde_euler(f,g,t,[1 1],opts); plot(t,y);
%       title(['Euler-Heun Method, dt = ' num2str(dt) ', \mu = ' num2str(mu)]);
%       txt = {['\sigma = ' num2str(sig(1))],['\sigma = ' num2str(sig(2))]};
%       legend(txt,2); legend boxoff; xlabel('t'); ylabel('y(t)');
%
%   Note:
%       SDEs are assumed to be in Stratonovich form by default. SDEs in Ito form
%       can be handled by setting the SDEType OPTIONS property to 'Ito' via the
%       SDESET function. These forms are generally not equivalent and will
%       converge to different solutions, so care should be taken to ensure that
%       the form of SDEFUN matches SDEType.
%
%       In the case of additve noise, i.e., when diffusion term, g(t,y), is not
%       a function of the state variables, both Ito and Stratonovich
%       interpretations are equivalent. If GFUN(T,Y) is specified as a floating
%       point matrix rather than a function handle then it is automatically
%       assumed to be constant. Alternatively, the ConstGFUN property can be set
%       to 'yes' to indicate that GFUN(T,Y) is constant and only needs to be
%       evaluated once.
%
%   See also:
%       Explicit SDE solvers:	SDE_MILSTEIN
%       Implicit SDE solvers:   
%       Stochastic processes:	SDE_GBM, SDE_OU
%       Option handling:        SDESET, SDEGET
%       SDE demos/validation:   SDE_EULER_VALIDATE, SDE_MILSTEIN_VALIDATE
%   	Other:                  FUNCTION_HANDLE, RANDSTREAM

%   SDE_EULER is an implementation of the order 0.5 strong (order 1.0 weak)
%   explicit Euler-Maruyama and Euler-Heun schemes, which are also order 0.5
%   strong Taylor approximation schemes. In the additive noise case, Euler
%   schemes improve to order 1.0 strong convergence.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, Created 10-28-10
%   Revision: 1.0, 4-21-12


solver = 'SDE_EULER';

% Check inputs and outputs
if nargin < 5
    if nargin < 4
        error(  'SDETools:sde_euler:NotEnoughInputs',...
                'Not enough input arguments.  See %s.',solver);
    end
    if isa(y0,'struct')
        error(  'SDETools:sde_euler:NotEnoughInputsOptions',...
               ['An SDE options structure was provided as the last '...
                'argument, but one of the first four input arguments is '...
                'missing.  See %s.'],solver);
    end
    options = [];
elseif isempty(options) && (ndims(options) ~= 2 || ...
        any(size(options) ~= 0) || ~(isstruct(options) || iscell(options) || ...
        isnumeric(options))) || ~isempty(options) && ~isstruct(options)	%#ok<*ISMAT>
	error(  'SDETools:sde_euler:InvalidSDESETStruct',...
            'Invalid SDE options structure.  See SDESET.');
end

% Handle solver arguments
[N D tspan tdir lt y0 fout gout h ConstStep dataType idxNonNegative ...
    NonNegative DiagonalNoise ScalarNoise ConstFFUN ConstGFUN Stratonovich ...
    RandFUN CustomRandFUN ResetAntithetic] ...
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
                    error(  'SDETools:sde_euler:RandFUNNot2DArray1',...
                           ['RandFUN must return a non-empty matrix of '...
                            'floating point values.  See %s.'],solver);
                end
                [m n] = size(W);
                if m ~= lt || n ~= D
                    error(  'SDETools:sde_euler:RandFUNDimensionMismatch1',...
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
                    error(  'SDETools:sde_euler:RandFUNNot2DArray2',...
                           ['RandFUN must return a non-empty matrix of '...
                            'floating point values.  See %s.'],solver);
                end
                [m n] = size(dW);
                if m ~= 1 || n ~= D
                    error(  'SDETools:sde_euler:RandFUNDimensionMismatch2',...
                           ['The specified alternative RandFUN did not '...
                            'output a %d by 1 column vector as requested.'...
                            '  See %s.',D,solver]);
                end
                dW = sh(1)*dW;
            end
        elseif D == N       % Store Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,D);
            if ndims(r) ~= 2 || isempty(r) || ~isfloat(r)
                error(  'SDETools:sde_euler:RandFUNNot2DArray3',...
                       ['RandFUN must return a non-empty matrix of floating '...
                        'point values.  See %s.'],solver);
            end
            [m n] = size(r);
            if m ~= lt-1 || n ~= D
                error(  'SDETools:sde_euler:RandFUNDimensionMismatch3',...
                       ['The specified alternative RandFUN did not output a '...
                        '%d by %d matrix as requested.'...
                        '   See %s.',D,lt-1,solver]);
            end
            if N == 1 || ConstStep
                Y(2:end,1:D) = sh.*r;
            else
                Y(2:end,1:D) = bsxfun(@times,sh,r);
            end
            clear r;    % remove large temporary variable to save memory
        end
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error(  'SDETools:sde_euler:RandFUNTooFewInputs',...
                       ['RandFUN must have at least two inputs.'...
                        '  See %s.'],solver);
            case 'MATLAB:TooManyOutputs'
                error(  'SDETools:sde_euler:RandFUNNoOutput',...
                       ['The output of RandFUN was not specified. RandFUN '...
                        'must return a non-empty matrix.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error(  'SDETools:sde_euler:RandFUNUnassignedOutput',...
                       ['The first output of RandFUN was not assigned.'...
                        '  See %s.'],solver);
            case 'MATLAB:minrhs'
                error(  'SDETools:sde_euler:RandFUNTooManyInputs',...
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
        if N == 1 && DiagonalNoise || ConstStep
            Y(2:end,1:D) = sh.*feval(RandFUN,lt-1,D);
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
if ConstFFUN && ConstGFUN && (D <= N || nargout == 2) && ~NonNegative	% no FOR loop needed
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
    % Optimized integration loop for 1-D scalar noise case
    if N == 1 && DiagonalNoise
        dt = h(1);  % Step size for first time step and fixed step-size
        Y(1) = y0;	% Set initial conditions
        
        % Use existing fout and gout to calculate the first time step, Y(2)
        if ConstGFUN
            Y(2) = y0+fout*dt+gout*Y(2);
        else
            if Stratonovich     % Use Euler-Heun step
                Y(2) = y0+fout*dt+0.5*(gout+feval(g,tspan(1),y0+gout*Y(2),varargin{:}))*Y(2);
            else
                Y(2) = y0+fout*dt+gout*Y(2);
            end
        end
        if NonNegative
            Y(2) = max(Y(2),0);
        end
        
        % Integration loop using Weiner increments stored in Y(i+1)
        for i = 2:lt-1
            if ~ConstStep
                dt = h(i);  % Step size
            end

            % Calculate next time step
            if ConstGFUN
              	Y(i+1) = Y(i)+feval(f,tspan(i),Y(i),varargin{:})*dt+gout*Y(i+1);
            else
                if ~ConstFFUN
                    fout = feval(f,tspan(i),Y(i),varargin{:})*dt;
                end
                gout = feval(g,tspan(i),Y(i),varargin{:});

                if Stratonovich     % Use Euler-Heun step
                    Y(i+1) = Y(i)+fout+0.5*(gout+feval(g,tspan(i),Y(i)+gout*Y(i+1),varargin{:}))*Y(i+1);
                else
                    Y(i+1) = Y(i)+fout+gout*Y(i+1);
                end
            end
            if NonNegative
                Y(i+1) = max(Y(i+1),0);
            end
        end
    else
        % step size and Wiener increment for first time-step and fixed step-size
        dt = h(1);	% Step size for first time step and fixed step-size
        if D > N
            if nargout == 2         % dW already exists if D > N && nargin == 1
                dW = W(2,:);
                W(2,:) = W(1,:)+dW;	% integrate Wiener increments
            else
                sdt = sh(1);
            end
        else
            dW = Y(2,1:D);
        end
        dW = dW(:);

        Y(1,:) = y0;	% Set initial conditions

        % Use existing fout and gout to calculate the first time step, Y(2,:)
        if ConstGFUN
            if DiagonalNoise
                Yi = y0+fout*dt+gout.*dW;
            else
                Yi = y0+fout*dt+gout*dW;
            end
        else
            if Stratonovich	% Use Euler-Heun step
                if DiagonalNoise
                    Yi = y0+fout*dt+0.5*(gout+feval(g,tspan(1),y0+gout.*dW,varargin{:})).*dW;
                else
                    Yi = y0+fout*dt+0.5*(gout+feval(g,tspan(1),y0+gout*dW,varargin{:}))*dW;
                end
            else
                if DiagonalNoise
                    Yi = y0+fout*dt+gout.*dW;
                else
                    Yi = y0+fout*dt+gout*dW;
                end
            end
        end
        if NonNegative
            Yi(idxNonNegative) = max(Yi(idxNonNegative),0);
        end
        Y(2,:) = Yi;

        % Integration loop using cached state, Yi, and increment, dW
        for i = 2:lt-1
            % Step size and Wiener increment
            if ~ConstStep
                dt = h(i);
                sdt = sh(i);
            end
            if D > N
                if nargout == 2
                    dW = W(i+1,:);
                    W(i+1,:) = W(i,:)+dW;      	% Integrate Wiener increments
                else
                    dW = sdt*feval(RandFUN,1,D);	% Generate Wiener increments
                end
            else
                dW = Y(i+1,1:D);
            end
            dW = dW(:);

            % Calculate next time step
            if ConstGFUN
                if DiagonalNoise
                    Yi = Yi+feval(f,tspan(i),Yi,varargin{:})*dt+gout.*dW;
                else
                    Yi = Yi+feval(f,tspan(i),Yi,varargin{:})*dt+gout*dW;
                end
            else
                if ~ConstFFUN
                    fout = feval(f,tspan(i),Yi,varargin{:})*dt;
                end
                gout = feval(g,tspan(i),Yi,varargin{:});

                if Stratonovich     % Use Euler-Heun step
                    if DiagonalNoise
                        Yi = Yi+fout+0.5*(gout+feval(g,tspan(i),Yi+gout.*dW,varargin{:})).*dW;
                    else
                        Yi = Yi+fout+0.5*(gout+feval(g,tspan(i),Yi+gout*dW,varargin{:}))*dW;
                    end
                else
                    if DiagonalNoise
                        Yi = Yi+fout+gout.*dW;
                    else
                        Yi = Yi+fout+gout*dW;
                    end
                end
            end
            if NonNegative
                Yi(idxNonNegative) = max(Yi(idxNonNegative),0);
            end
            Y(i+1,:) = Yi;
        end
    end
end

% Reset antihetic property if global stream was used
if ResetAntithetic
    try
        Stream = RandStream.getGlobalStream;
    catch                                       %#ok<CTCH>
        Stream = RandStream.getDefaultStream;	%#ok<GETRS>
    end
    set(Stream,'Antithetic',~Stream.Antithetic);
end