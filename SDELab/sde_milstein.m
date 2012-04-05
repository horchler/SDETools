function [Y W] = sde_milstein(f,g,tspan,y0,options,varargin)
%SDE_MILSTEIN  Solve stochastic differential equations, Milstein method.
%   YOUT = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL]
%   integrates the N-dimensional system of stochastic differential equations
%   dy = f(t,y)*dt + g(t,y)*dW with N-dimensional diagonal noise from time T0 to
%   TFINAL (all increasing or all decreasing) with initial conditions Y0. FFUN
%   and GFUN are function handles or floating point matrices. For scalar T and
%   vector Y, FFUN(T,Y) and GFUN(T,Y) return column vectors corresponding to
%   f(t,y) and g(t,y), the deterministic and stochastic parts of the SDE,
%   respectively. GFUN(T,Y) can also return a scalar, in which case this value
%   is used across all N dimensions. TSPAN is a length M vector. Y0 is a length
%   N vector. Each row in the M-by-N solution array YOUT corresponds to a time
%   in TSPAN.
%
%   [YOUT, W] = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0) outputs the M-by-N matrix W of
%   integrated Weiner increments that were used by the solver. Each row of W to
%   a time in TSPAN.
%
%   YOUT = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0,OPTIONS) solves the above with
%   default integration properties replaced by values in OPTIONS, an argument
%   created with the SDESET function. See SDESET for details. The type of SDE to
%   be integrated, 'Ito' or the default 'Stratonovich', can be specified via the
%   SDEType property. By default, the derivative-free Milstein method is
%   applied, but if the if the DGFUN property is set to a function handle or
%   floating point matrix to specify the derivative of the noise function
%   g(t,y), then the general Milstein method is used. Another commonly used
%   option is manually specifying the random seed with the RandSeed property and
%   creates a separate random number stream rather than using the default one.
%
%   Example:
%       y = sde_milstein(@f1,@g1,0:0.01:1,[1 0]);
%       solves the 2-D Stratonovich SDE system dy = f1(t,y)*dt + g1(t,y)*dW,
%       using the derivative-free Milstein method and the default random number
%       stream.
%
%   Notes:
%       SDEs are assumed to be in Stratonovich form by default. SDEs in Ito form
%       can be handled by setting the 'SDEType' OPTIONS property to 'Ito' via
%       the SDESET function. These forms are generally not equivalent and will
%       converge to different solutions, so care should be taken to ensure that
%       the form of SDEFUN matches 'SDEType'.
%
%       In the case of additve noise, i.e., when g(t,y) is constant, both Ito
%       and Stratonovich interpretations are equivalent and the Milstein method
%       reduces to the Euler-Maruyama method implemented in SDE_EULER. If
%       GFUN(T,Y) is specified as a floating point matrix rather than a function
%       handle then it is automatically assumed to be constant. Alternatively,
%       the ConstGFUN property can be set to 'yes' to indicate that GFUN(T,Y) is
%       constant and only needs to be evaluated once.
%
%   See also:
%       other SDE solvers:  SDE_EULER
%       implicit SDEs:      
%       options handling:   SDESET, SDEGET
%       SDE examples:       
%       function handles:   FUNCTION_HANDLE

%   SDE_MILSTEIN is an implementation of the order 1.0 strong (order 1.0 weak)
%   explicit Milstein scheme, which is also an order 1.0 strong Taylor
%   approximation scheme.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, 10-25-10
%   Revision: 1.0, 4-4-12


solver = 'SDE_MILSTEIN';

% Check inputs and outputs
if nargin < 5
    if nargin < 4
        error('SDELab:sde_milstein:NotEnoughInputs',...
              'Not enough input arguments.  See %s.',solver);
    end
    if isa(y0,'struct')
        error('SDELab:sde_milstein:NotEnoughInputsOptions',...
             ['An SDE options structure was provided as the last argument, '...
              'but one of the first four input arguments is missing.'...
              '  See %s.'],solver);
    end
    options = [];
elseif isempty(options) && (ndims(options) ~= 2 || ...
        any(size(options) ~= 0) || ~(isstruct(options) || iscell(options) || ...
        isnumeric(options))) || ~isempty(options) && ~isstruct(options)
	error('SDELab:sde_milstein:InvalidSDESETStruct',...
          'Invalid SDE options structure.  See SDESET.');
end

% Handle solver arguments
[N D tspan tdir lt y0 fout gout h ConstStep dataType idxNonNegative ...
    NonNegative DiagonalNoise ScalarNoise ConstFFUN ConstGFUN Stratonovich ...
    RandFUN CustomRandFUN] = sdearguments(solver,f,g,tspan,y0,options,varargin);

% If optional derivative function property is set
if ~ConstGFUN
    dg = sdeget(options,'DGFUN',[],'fast');

    % Ensure DGFUN is function handle, or matrix for constant function
    if ~isa(dg,'function_handle')
        if isempty(dg) && ndims(dg) == 2 && all(size(dg) == 0) && isnumeric(dg)
            Derivative = false;
        elseif ~isempty(dg) && ndims(dg) == 2 && isfloat(dg)
            dgout = dg;
            ConstDGFUN = true;
            Derivative = true;
            if all(dgout == 0)
                ConstGFUN = true;
            end
        else
            error(  'SDELab:sde_milstein:InvalidDGFUN',...
                   ['The input DGFUN must be a function handle, a matrix of '...
                    'singles or doubles, or an empty numeric matrix, [], '...
                    'denoting no stochastic derivative function.'...
                    '  See %s.'],solver);
        end
    else
        % Check output of DGFUN and save it
        try
            dgout = feval(dg,t0,y0,varargin{:});
        catch err
            switch err.identifier
                case 'MATLAB:TooManyInputs'
                    error(  'SDELab:sde_milstein:DGFUNTooFewInputs',...
                           ['DGFUN must have at least two inputs.'...
                            '  See %s.'],solver);
                case 'MATLAB:TooManyOutputs'
                    error(  'SDELab:sde_milstein:DGFUNNoOutput',...
                           ['The output of DGFUN was not specified. DGFUN '...
                            'must return a non-empty matrix.  See %s.'],solver);
                case 'MATLAB:unassignedOutputs'
                    error(  'SDELab:sde_milstein:DGFUNUnassignedOutput',...
                           ['The first output of DGFUN was not assigned.'...
                            '  See %s.'],solver);
                case 'MATLAB:minrhs'
                    error(  'SDELab:sde_milstein:DGFUNTooManyInputs',...
                           ['DGFUN requires one or more input arguments '...
                            '(parameters) that were not supplied.'...
                            '  See %s.'],solver);
                otherwise
                    rethrow(err);
            end
        end
        % If stochastic derivative function specified as constant
        ConstDGFUN = strcmp(sdeget(options,'ConstDGFUN','no','flag'),'yes');
        Derivative = true;
    end

    % Ensure that size and type of output of DGFUN is consistent
    if Derivative
        if ndims(dgout) ~= 2 || isempty(dgout) || ~isfloat(dgout)
            error(  'SDELab:sde_milstein:DGFUNNot2DArray',...
                   ['DGFUN must return a non-empty matrix of floating point '...
                    'values.  See %s.'],solver);
        end
        [m n] = size(dgout);
        if DiagonalNoise
            if n ~= 1 || ~(m == N || m == 1)
                error(  'SDELab:sde_milstein:DGFUNNotColumnVector',...
                       ['For diagonal noise, DGFUN must return a scalar '...
                        'value or a non-empty column vector the same length '...
                        'as Y0.  See %s.'],solver);
            end
        else
            if m ~= N
                error(  'SDELab:sde_milstein:DGFUNDimensionMismatchNonDiag',...
                       ['For non-diagonal noise, DGFUN must return a scalar '...
                        'value or a non-empty matrix with the same number '...
                        'of rows as the length as Y0.  See %s.'],solver);
            end
        end
    end
end

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
                    error(  'SDELab:sde_milstein:RandFUNNot2DArray1',...
                           ['RandFUN must return a non-empty matrix of '...
                            'floating point values.  See %s.'],solver);
                end
                [m n] = size(W);
                if m ~= lt || n ~= D
                    error(  'SDELab:sde_milstein:RandFUNDimensionMismatch1',...
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
                    error(  'SDELab:sde_milstein:RandFUNNot2DArray2',...
                           ['RandFUN must return a non-empty matrix of '...
                            'floating point values.  See %s.'],solver);
                end
                [m n] = size(dW);
                if m ~= 1 || n ~= D
                    error(  'SDELab:sde_milstein:RandFUNDimensionMismatch2',...
                           ['The specified alternative RandFUN did not '...
                            'output a %d by 1 column vector as requested.'...
                            '  See %s.',D,solver]);
                end
                dW = sh(1)*dW;
            end
        elseif D == N       % Store Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,D);
            if ndims(r) ~= 2 || isempty(r) || ~isfloat(r)
                error(  'SDELab:sde_milstein:RandFUNNot2DArray3',...
                       ['RandFUN must return a non-empty matrix of floating '...
                        'point values.  See %s.'],solver);
            end
            [m n] = size(r);
            if m ~= lt-1 || n ~= D
                error(  'SDELab:sde_milstein:RandFUNDimensionMismatch3',...
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
                error(  'SDELab:sde_milstein:RandFUNTooFewInputs',...
                       ['RandFUN must have at least two inputs.'...
                        '  See %s.'],solver);
            case 'MATLAB:TooManyOutputs'
                error(  'SDELab:sde_milstein:RandFUNNoOutput',...
                       ['The output of RandFUN was not specified. RandFUN '...
                        'must return a non-empty matrix.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error(  'SDELab:sde_milstein:RandFUNUnassignedOutput',...
                       ['The first output of RandFUN was not assigned.'...
                        '  See %s.'],solver);
            case 'MATLAB:minrhs'
                error(  'SDELab:sde_milstein:RandFUNTooManyInputs',...
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

% Integration loop
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

        % Use existing fout and gout to calculate the first time step, Y(2,:)
        if ConstGFUN
            Y(2) = y0+fout*dt+gout*Y(2);
        else
            if Derivative
                if Stratonovich
                    dW2 = 0.5*Y(2)^2;
                else
                    dW2 = 0.5*(Y(2)^2-dt);
                end
                Y(2) = y0+fout*dt+gout*(Y(2)+dgout*dW2);
            else
                sdt = sh(1);
                if Stratonovich
                    dW2 = (0.5/sdt)*Y(2)^2;
                else
                    dW2 = (0.5/sdt)*(Y(2)^2-dt);
                end
              	Ybar = y0+gout*Y(2);
                Y(2) = Ybar+fout*dt+(feval(g,tspan(1),Ybar,varargin{:})-gout)*dW2;
            end
        end
        if NonNegative
            Y(2) = max(Y(2),0);
        end
        
        % Integration loop using Weiner increments stored in Y(i+1)
        for i = 2:lt-1
            % Step size and Wiener increments
            if ~ConstStep
                dt = h(i);
                sdt = sh(i);
            end

            % Calculate next time step
            if ConstGFUN
              	Y(i+1) = Y(i)+feval(f,tspan(i),Y(i),varargin{:})*dt+gout*Y(i+1);
            else
                if ~ConstFFUN
                    fout = feval(f,tspan(i),Y(i),varargin{:})*dt;
                end

                gout = feval(g,tspan(i),Y(i),varargin{:});
                if Derivative   % Use Milstein with derivative if specified
                    if Stratonovich
                        dW2 = 0.5*Y(i+1)^2;
                    else
                        dW2 = 0.5*(Y(i+1)^2-dt);
                    end
                    if ~ConstDGFUN
                        dgout = feval(dg,tspan(i),Y(i),varargin{:});
                    end
                    Y(i+1) = Y(i)+fout+gout*(Y(i+1)+dgout*dW2);
                else            % Otherwise use derivative-free Milstein method
                    if Stratonovich
                        dW2 = (0.5/sdt)*Y(i+1)^2;
                    else
                        dW2 = (0.5/sdt)*(Y(i+1)^2-dt);
                    end
                 	Ybar = Y(i)+gout*Y(i+1);
                	Y(i+1) = Ybar+fout+(feval(g,tspan(i),Ybar,varargin{:})-gout)*dW2;
                end
            end
            if NonNegative
                Y(i+1) = max(Y(i+1),0);
            end
        end
    else
        % step size and Wiener increment for first time step
        dt = h(1);
        sdt = sh(1);
        if D > N 
            if nargout == 2         % dW already exists if D > N && nargin == 1
                dW = W(2,:);
                W(2,:) = W(1,:)+dW;	% integrate Wiener increments
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
            if Derivative
                if Stratonovich
                    dW2 = 0.5*dW.^2;
                else
                    dW2 = 0.5*(dW.^2-dt);
                end
                if DiagonalNoise
                    Yi = y0+fout*dt+gout.*(dW+dgout.*dW2);
                else
                    Yi = y0+fout*dt+gout*(dW+dgout*dW2);
                end
            else
                if Stratonovich
                    dW2 = (0.5/sdt)*dW.^2;
                else
                    dW2 = (0.5/sdt)*(dW.^2-dt);
                end
                if DiagonalNoise
                    Ybar = y0+gout.*dW;
                    Yi = Ybar+fout*dt+(feval(g,tspan(1),Ybar,varargin{:})-gout).*dW2;
                else
                    Ybar = y0+gout*dW;
                    Yi = Ybar+fout*dt+(feval(g,tspan(1),Ybar,varargin{:})-gout)*dW2;
                end
            end
        end
        if NonNegative
            Yi(idxNonNegative) = max(Yi(idxNonNegative),0);
        end
        Y(2,:) = Yi;

        % Integration loop using cached state, Yi, and increment, dW
        for i = 2:lt-1
            % Step size and Wiener increments
            if ~ConstStep
                dt = h(i);
                sdt = sh(i);
            end
            if D > N
                if nargout == 2
                    dW = W(i+1,:);
                    W(i+1,:) = W(i,:)+dW;     	% Integrate Wiener increments
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
                if Derivative   % Use Milstein with derivative function if specified
                    if Stratonovich
                        dW2 = 0.5*dW.^2;
                    else
                        dW2 = 0.5*(dW.^2-dt);
                    end
                    if ~ConstDGFUN
                        dgout = feval(dg,tspan(i),Yi,varargin{:});
                    end
                    if DiagonalNoise
                        Yi = Yi+fout+gout.*(dW+dgout.*dW2);
                    else
                        Yi = Yi+fout+gout*(dW+dgout*dW2);
                    end
                else            % Otherwise use derivative-free Milstein method
                    if Stratonovich
                        dW2 = (0.5/sdt)*dW.^2;
                    else
                        dW2 = (0.5/sdt)*(dW.^2-dt);
                    end
                    if DiagonalNoise
                        Ybar = Yi+gout.*dW;
                        Yi = Ybar+fout+(feval(g,tspan(i),Ybar,varargin{:})-gout).*dW2;
                    else
                        Ybar = Yi+gout*dW;
                        Yi = Ybar+fout+(feval(g,tspan(i),Ybar,varargin{:})-gout)*dW2;
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