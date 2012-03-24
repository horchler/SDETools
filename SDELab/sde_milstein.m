function [Y W] = sde_milstein(f,g,tspan,y0,options,varargin)
%SDE_MILSTEIN  Solve stochastic differential equations, Milstein method.
%   YOUT = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL]
%   integrates the system of stochastic differential equations dy = f(t,y)*dt +
%   g(t,y)*dW with diagonal noise from time T0 to TFINAL (all increasing or all
%   decreasing with fixed step size) with initial conditions Y0. FFUN and GFUN
%   are function handles. For a scalar T and a vector Y, FFUN(T,Y) and GFUN(T,Y)
%   return column vectors corresponding to f(t,y) and g(t,y), the deterministic
%   and stochastic parts of the SDE, respectively. Each row in the solution
%   array YOUT corresponds to a time in the input vector TSPAN.
%
%   YOUT = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0,OPTIONS) solves the above with
%   default integration properties replaced by values in OPTIONS, an argument
%   created with the SDESET function. See SDESET for details. The type of SDE to
%   be integrated, 'Ito' or the default 'Stratonovich', can be specified via the
%   SDEType property. By default, the derivative-free Milstein method is 
%   applied, but if the Derivative property is set to a function handle that
%   refers to the derivative of the noise function g(t,y), the full Milstein
%   method is used. Another commonly used option is manually specifying the
%   random seed with the RandSeed property and creates a separate random number
%   stream rather than using the default one.
%
%   [YOUT W] = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0) outputs the matrix W of
%   integrated Weiner increments that were used by the solver. W is LENGTH(Y0)
%   rows by LENGTH(TSPAN) columns, corresponding to [T0 T1 T2 ... TFINAL].
%
%   Example
%       y = sde_milstein(@f1,@g1,0:0.01:1,[1 0]);
%       solves the 2-D Stratonovich SDE system dy = f1(t,y)*dt + g1(t,y)*dW,
%       using the derivative-free Milstein method and the default random number
%       stream.
%
%   NOTE:
%       SDEs are assumed to be in Stratonovich form by default. SDEs in Ito form
%       can be handled by setting the 'SDEType' OPTIONS property to 'Ito' via
%       the SDESET function. These forms are generally not equivalent and will
%       converge to different solutions, so care should be taken to ensure that
%       the form of SDEFUN matches 'SDEType'.
%
%   See also
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
%   Revision: 1.0


% Check inputs
if nargin < 5
	options = [];
    if nargin < 4
        error(  'MATLAB:sde_milstein:NotEnoughInputs',...
                'Not enough input arguments.  See SDE_MILSTEIN.');
    end
end

% Ensure first two inputs are function handles
if ~isa(f,'function_handle')
    error(  'MATLAB:sde_milstein:NotAFunctionHandle',...
            'Input FFUN must be a function handle.  See SDE_MILSTEIN.');
end
if ~isa(g,'function_handle')
    error(  'MATLAB:sde_milstein:NotAFunctionHandle',...
            'Input GFUN must be a function handle.  See SDE_MILSTEIN.');
end

lt = length(tspan);     % number of time steps
if lt < 2
    error(  'MATLAB:sde_milstein:BadInputSize',...
            'Input vector TSPAN must have length >= 2.  See SDE_MILSTEIN.');
end
N = length(y0);         % number of state variables
if N < 1
    error(  'MATLAB:sde_milstein:BadInputSize',...
            'Input vector Y0 must have length >= 1.  See SDE_MILSTEIN.');
end

h = tspan(2)-tspan(1);  % length of fixed step size
sh = sqrt(h);

% Check output of SDE functions and save it
fout = feval(f,tspan(1),y0,varargin{:});
if ~isvector(fout) || length(fout) ~= N
    error(  'MATLAB:sde_milstein:DimensionMismatch',...
            'FFUN must return vector of LENGTH(Y0).  See SDE_MILSTEIN.');
end
yout = fout*h;

% Dimension of stochastic function
DiagonalNoise = strcmp(sdeget(options,'DiagonalNoise','yes','fast'),'yes');

gout = feval(g,tspan(1),y0,varargin{:});
if DiagonalNoise
    if ~isvector(gout) || length(gout) ~= N
        error(  'MATLAB:sde_milstein:DimensionMismatch',...
                'GFUN must return vector of LENGTH(Y0) if noise is diagonal.  See SDE_MILSTEIN.');
    end
else
    if size(gout,1) ~= N || size(gout,2) ~= N
        error(  'MATLAB:sde_milstein:DimensionMismatch',...
                'GFUN must return square matrix of dimension LENGTH(Y0).  See SDE_MILSTEIN.');
    end
end

% If optional derivative function property is set
dg = sdeget(options,'DGFUN',[],'fast');
Derivative = ~isempty(dg);
if Derivative
    if ~isa(dg,'function_handle')
        error(  'MATLAB:sde_milstein:NotAFunctionHandle',...
                'DGFUN must be a function handle.  See SDE_MILSTEIN.');
    end
    dgout = feval(dg,tspan(1),y0,varargin{:});
    if DiagonalNoise
        if ~isvector(dgout) || length(dgout) ~= N
            error(  'MATLAB:sde_milstein:DimensionMismatch',...
                    'DGFUN must return vector of LENGTH(Y0).  See SDE_MILSTEIN.');
        end
    else
        if size(dgout,1) ~= N || size(dgout,2) ~= N
            error(  'MATLAB:sde_milstein:DimensionMismatch',...
                    'DGFUN must return square matrix of dimension LENGTH(Y0).  See SDE_MILSTEIN.');
        end
    end
end

% If noise is specified as scalar, 1-D
ScalarNoise = strcmp(sdeget(options,'ScalarNoise','no','fast'),'yes');

% Calculate Wiener increments from normal variates and store them in state array
RandFUN = sdeget(options,'RandFUN',[],'fast');
if ~isempty(RandFUN)	% Use alternative random number generator
    if ~isa(RandFUN,'function_handle')
        error(  'MATLAB:sde_milstein:NotAFunctionHandle',...
                'RANDFUN must be a function handle.  See SDE_MILSTEIN.');
    end
    if ScalarNoise
        dW = [0 sh*feval(RandFUN,1,lt-1)];
        if size(dW,1) ~= 1 || size(dW,2) ~= lt
            error(  'MATLAB:sde_milstein:DimensionMismatch',...
                   ['Output RANDFUN of must be 1 row by LENGTH(TSPAN)-1 '...
                    'columns.  See SDE_MILSTEIN.']);
        end
        Y = zeros(N,lt);
        Y(1,:) = dW;
    else
        Y = [zeros(N,1) sh*feval(RandFUN,N,lt-1)];
        if size(Y,1) ~= N || size(Y,2) ~= lt
            error(  'MATLAB:sde_milstein:DimensionMismatch',...
                   ['Output RANDFUN of must be LENGTH(Y0) rows by LENGTH(TSPAN)-1 '...
                    'columns.  See SDE_MILSTEIN.']);
        end
    end
else	% Use Matlab's random number generator for normal variates
    RandSeed = sdeget(options,'RandSeed',[],'fast');
    if ~isempty(RandSeed)
        if RandSeed >= 2^32 || RandSeed < 0
            error(	'MATLAB:sde_milstein:ImproperValue',...
                    'RANDSEED must be positive and less than 2^32.  See SDE_MILSTEIN.');
        end
        Stream = RandStream.create('mt19937ar','seed',RandSeed);    % Create stream
    else
        Stream = RandStream.getDefaultStream;   % Use default stream
    end

    % Pre-calculate Wiener increments, optionally using antithetic variates method
    if N > 1 && strcmp(sdeget(options,'Antithetic','no','fast'),'yes')
        M = round(0.5*N);
        Y = zeros(N,lt);
        Y(1:M,2:end) = sh*randn(Stream,M,lt-1);
        Y(M+1:end,2:end) = -Y(1:M-mod(N,2),2:end);
    else
        if ScalarNoise
            Y = zeros(N,lt);
            Y(1,2:end) = sh*randn(Stream,1,lt-1);
        else
            Y = [zeros(N,1) sh*randn(Stream,N,lt-1)];
        end
    end
end

% Only allocate W matrix if requested as output
if nargout > 1
    if ScalarNoise
        W = cumsum(Y(1,:));
    else
        W = cumsum(Y,2);
    end
end

% If deterministic or stochastic functions functions specified as constant
ConstFFUN = strcmp(sdeget(options,'ConstFFUN','no','fast'),'yes');
ConstGFUN = strcmp(sdeget(options,'ConstGFUN','no','fast'),'yes');

% Square of noise term is dependent on if SDE is Stratonovich or Ito form
isIto = strcmp(sdeget(options,'SDEType','Stratonovich','fast'),'Ito');

% Integration loop
Y(:,1) = y0;
if ConstFFUN && ConstGFUN
    if ScalarNoise
        Y(:,2:end) = bsxfun(@plus,y0,cumsum(bsxfun(@plus,yout,gout*Y(1,2:end)),2));
    elseif DiagonalNoise
        Y(:,2:end) = bsxfun(@plus,y0,cumsum(bsxfun(@plus,yout,bsxfun(@times,gout,Y(:,2:end))),2));
    else
        Y(:,2:end) = bsxfun(@plus,y0,cumsum(bsxfun(@plus,yout,gout*Y(:,2:end)),2));
    end
else
    if ConstGFUN
        if ScalarNoise
            Y(:,2:end) = gout*Y(1,2:end);
        elseif DiagonalNoise
            Y(:,2:end) = bsxfun(@times,gout,Y(:,2:end));
        else
            Y(:,2:end) = gout*Y(:,2:end);
        end
        Y(:,2) = y0 + yout + Y(:,2);
    else
        if Derivative
            if DiagonalNoise
                dgout = gout.*dgout;
            else
                dgout = gout*dgout;
            end
        else
            if DiagonalNoise
                dgout = feval(g,tspan(1),y0+yout+gout*sh,varargin{:})-gout;
            else
                dgout = feval(g,tspan(1),y0+yout+diag(gout,0)*sh,varargin{:})-gout;%wrong
            end
        end
        if ScalarNoise
            dW2 = 0.5*(Y(1,2).^2-h*isIto)*(~Derivative/sh+Derivative);
            Y(:,2) = y0 + yout + gout*Y(1,2) + dgout*dW2;
        else
            dW2 = 0.5*(Y(:,2).^2-h*isIto)*(~Derivative/sh+Derivative);
            if DiagonalNoise
                Y(:,2) = y0 + yout + gout.*Y(:,2) + dgout.*dW2;
            else
                Y(:,2) = y0 + yout + gout*Y(:,2) + dgout*dW2;%wrong
            end
        end
    end
    for i=2:lt-1
        if ~ConstFFUN
            yout = feval(f,tspan(i),Y(:,i),varargin{:})*h;
        end
        if ConstGFUN
            Y(:,i+1) = Y(:,i) + yout + Y(:,i+1);
        else
            gout = feval(g,tspan(i),Y(:,i),varargin{:});
            if Derivative   % Use Milstein with derivative function if specified
                if DiagonalNoise
                    dgout = gout.*feval(dg,tspan(i),Y(:,i),varargin{:});
                else
                    dgout = gout*feval(dg,tspan(i),Y(:,i),varargin{:});
                end
            else            % Otherwise use derivative-free Milstein method
                if DiagonalNoise
                    dgout = feval(g,tspan(i),Y(:,i)+yout+gout*sh,varargin{:})-gout;
                else
                    dgout = feval(g,tspan(1),y0+yout+diag(gout,0)*sh,varargin{:})-gout;%wrong
                end
            end
            if ScalarNoise
                dW2 = 0.5*(Y(1,i+1).^2-h*isIto)*(~Derivative/sh+Derivative);
                Y(:,i+1) = Y(:,i) + yout + gout*Y(1,i+1) + dgout*dW2;
            else
                dW2 = 0.5*(Y(:,i+1).^2-h*isIto)*(~Derivative/sh+Derivative);
                if DiagonalNoise
                    Y(:,i+1) = Y(:,i) + yout + gout.*Y(:,i+1) + dgout.*dW2;
                else
                    Y(:,i+1) = Y(:,i) + yout + gout*Y(:,i+1) + dgout*dW2;%wrong
                end
            end
        end
    end
end