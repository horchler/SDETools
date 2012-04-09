function [Y W] = sde_gbm(mu,sig,tspan,y0,options)
%SDE_GBM  Geometric Brownian motion process, analytic solution.
%   YOUT = SDE_GBM(MU,SIG,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL] returns the
%   analytic solution of the system of stochastic differential equations for
%   geometric Brownian motion, dY = MU*Y*dt + SIG*Y*dW, with diagonal noise from
%   time T0 to TFINAL (all increasing or all decreasing with arbitrary step
%   size) with initial conditions Y0. The drift parameter MU and the diffusion
%   (volatility) parameter SIG are scalars or vectors of LENGTH(Y0). Each row in
%   the solution array YOUT corresponds to a time in the input vector TSPAN.
%
%   [YOUT, W] = SDE_GBM(MU,SIG,TSPAN,Y0,...) outputs the matrix W of integrated
%   Weiner increments that were used. W is LENGTH(Y0) rows by LENGTH(TSPAN)
%   columns, corresponding to [T0 T1 T2 ... TFINAL].
%
%   [...] = SDE_GBM(MU,SIG,TSPAN,Y0,OPTIONS) returns the above with default
%   properties replaced by values in OPTIONS, an argument created with the
%   SDESET function. See SDESET for details. The type of SDE to be solved, 'Ito'
%   or the default 'Stratonovich', can be specified via the SDEType property.
%   Another commonly used option is to manually specify the random seed via the
%   RandSeed property, which creates a new random number stream, instead of
%   using the default stream, to generate the Wiener increments.
%
%   Example:
%       % Compare analytical and simulated Geometric Brownian motion
%       npaths = 10; dt = 1e-2; t = 0:dt:1; y0 = ones(1,npaths);
%       mu = 1.5; sig = 0.1:0.1:0.1*npaths;
%       opts = sdeset('RandSeed',1,'SDEType','Ito');
%       y1 = sde_gbm(mu,sig,t,y0,opts);
%       y2 = sde_euler(@(t,y)(mu+sig'.^2/2).*y,@(t,y)sig'.*y,t,y0,opts);
%       h = plot(t,y1,'b',t,y2,'r'); xlabel('t'); ylabel('y(t)');
%       mustr = num2str(mu); npstr = num2str(npaths); dtstr = num2str(dt);
%       txt = {'Analytical solution',['Numerical solution, dt = ' dtstr]};
%       legend(h([1 end]),txt,2); legend boxoff;
%       title(['Geometric Brownian motion, ' npstr ' paths, \mu = ' mustr]);
%
%   Note:
%       The solution assumes the Stratonovich form by default. The Ito form can
%       be handled by setting the SDEType OPTIONS property to 'Ito' via the
%       SDESET function. These forms are generally not equivalent and will
%       converge to different solutions, so care should be taken to ensure that
%       the form of SDEType matches the that of the desired solution.
%
%   See also:
%       Explicit SDE solvers:	SDE_EULER, SDE_MILSTEIN
%       Implicit SDE solvers:   
%       Stochastic processes:	SDE_OU
%       Option handling:        SDESET, SDEGET
%       SDE demos/validation:   SDE_EULER_VALIDATE, SDE_MILSTEIN_VALIDATE
%   	Other:                  FUNCTION_HANDLE, RANDSTREAM

%   The conditional analytic solutions Y = Y0*exp((MU-SIG^2/2)*t+SIG*W) and
%   Y = Y0*exp(MU*t+SIG*W) are used for Ito and Stratonovich cases,
%   respectively, where W is a standard Wiener process.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9@case.edu, Created 4-4-12
%   Revision: 1.0, 4-9-12


func = 'SDE_GBM';

% Check inputs and outputs
if nargin < 5
    if nargin < 4
        error(  'SDELab:sde_gbm:NotEnoughInputs',...
                'Not enough input arguments.  See %s.',func);
    end
    if isa(y0,'struct')
        error(  'SDELab:sde_gbm:NotEnoughInputsOptions',...
               ['An SDE options structure was provided as the last '...
                'argument, but one of the first four input arguments is '...
                'missing.  See %s.'],func);
    end
    options = [];
elseif isempty(options) && (ndims(options) ~= 2 || ...
        any(size(options) ~= 0) || ~(isstruct(options) || iscell(options) || ...
        isnumeric(options))) || ~isempty(options) && ~isstruct(options)
	error(  'SDELab:sde_gbm:InvalidSDESETStruct',...
            'Invalid SDE options structure.  See SDESET.');
end

% Check mu and sig types
if isempty(mu) || ~isfloat(mu) || ~isvector(mu)
    error(  'SDELab:sde_gbm:MuEmptyOrNotFloatVector',...
           ['The drift parameter, MU, must be non-empty vector of singles '...
            'or doubles.  See %s.'],func);
end
if isempty(sig) || ~isfloat(sig) || ~isvector(sig)
    error(  'SDELab:sde_gbm:SigEmptyOrNotFloatVector',...
           ['The diffusion (volatility) parameter, SIG, must be non-empty '...
            'vector of singles or doubles.  See %s.'],func);
end

% Determine the dominant data type, single or double
dataType = superiorfloat(mu,sig,tspan,y0);
if ~all(strcmp(dataType,{class(mu),class(sig),class(tspan),class(y0)}))
    warning( 'SDELab:sde_gbm:InconsistentDataType',...
            ['Mixture of single and double data for inputs MU, SIG, TSPAN, '...
             'and Y0.']);
end

% Handle function arguments
[N tspan tdir lt y0 h ConstStep Stratonovich RandFUN CustomRandFUN] ...
	= sdearguments_special(func,tspan,y0,options,dataType);

% Check mu and sig
if ~any(length(mu) == [1 N])
    error(  'SDELab:sde_gbm:MuDimensionMismatch',...
           ['The drift parameter, MU, must be a scalar or a vector the same '...
            'length as Y0.  See %s.'],func);
end
if ~any(length(sig) == [1 N])
    error(  'SDELab:sde_gbm:SigDimensionMismatch',...
           ['The diffusion parameter, SIG, must be a scalar or a vector the '...
            'same length as Y0.  See %s.'],func);
end
if any(sig < 0)
    error(  'SDELab:sde_gbm:SigNegative',...
           ['The diffusion (volatility) parameter, SIG, must be greater '...
            'than or equal to zero.  See %s.'],func);
end

Y = zeros(lt,N,dataType);   % State array

% Diffusion parameters aren't all zero
if ~all(sig == 0)
    % Calculate Wiener increments from normal variates, store in state if possible
    if CustomRandFUN    % check output of alternative RandFUN
        try
            % Store Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,N);
            if ndims(r) ~= 2 || isempty(r) || ~isfloat(r)
                error(  'SDELab:sde_gbm:RandFUNNot2DArray3',...
                       ['RandFUN must return a non-empty matrix of floating '...
                        'point values.  See %s.'],solver);
            end
            [m n] = size(r);
            if m ~= lt-1 || n ~= N
                error(  'SDELab:sde_gbm:RandFUNDimensionMismatch3',...
                       ['The specified alternative RandFUN did not output a '...
                        '%d by %d matrix as requested.'...
                        '   See %s.',N,lt-1,solver]);
            end
            if N == 1 || ConstStep
                Y(2:end,:) = tdir*sqrt(h).*r;
            else
                Y(2:end,:) = bsxfun(@times,tdir*sqrt(h),r);
            end
            clear r;    % remove large temporary variable to save memory
        catch err
            switch err.identifier
                case 'MATLAB:TooManyInputs'
                    error(  'SDELab:sde_gbm:RandFUNTooFewInputs',...
                           ['RandFUN must have at least two inputs.'...
                            '  See %s.'],solver);
                case 'MATLAB:TooManyOutputs'
                    error(  'SDELab:sde_gbm:RandFUNNoOutput',...
                           ['The output of RandFUN was not specified. '...
                            'RandFUN must return a non-empty matrix.'...
                            '  See %s.'],solver);
                case 'MATLAB:unassignedOutputs'
                    error(  'SDELab:sde_gbm:RandFUNUnassignedOutput',...
                           ['The first output of RandFUN was not assigned.'...
                            '  See %s.'],solver);
                case 'MATLAB:minrhs'
                    error(  'SDELab:sde_gbm:RandFUNTooManyInputs',...
                           ['RandFUN must not require more than two inputs.'...
                            '  See %s.'],solver);
                otherwise
                    rethrow(err);
            end
        end
    else    % No error checking needed if default RANDN used
        % Store Wiener increments in Y
        if N == 1 || ConstStep
            Y(2:end,:) = tdir*sqrt(h).*feval(RandFUN,lt-1,N);
        else
            Y(2:end,:) = bsxfun(@times,tdir*sqrt(h),feval(RandFUN,lt-1,N));
        end
    end
    
    if Stratonovich
        % Only allocate W matrix if requested as output
        if nargout == 2
            W = cumsum(Y,1);
            if N == 1
                Y = exp(tspan*(mu-0.5*sig^2)+sig*W)*y0;
            else
                if isscalar(mu) && isscalar(sig)
                    Y = bsxfun(@times,y0',exp(bsxfun(@plus,tspan*(mu-0.5*sig^2),sig*W)));
                elseif isscalar(sig)
                    Y = bsxfun(@times,y0',exp(tspan*(mu(:)'-0.5*sig^2)+sig*W));
                else
                    sig = sig(:)';
                    Y = bsxfun(@times,y0',exp(tspan*(mu(:)'-0.5*sig.^2)+bsxfun(@times,sig,W)));
                end
            end
        else
            if N == 1
                Y = exp(tspan*(mu-0.5*sig^2)+sig*cumsum(Y,1))*y0;
            else
                if isscalar(mu) && isscalar(sig)
                    Y = bsxfun(@times,y0',exp(bsxfun(@plus,tspan*(mu-0.5*sig^2),sig*cumsum(Y,1))));
                elseif isscalar(sig)
                    Y = bsxfun(@times,y0',exp(tspan*(mu(:)'-0.5*sig^2)+sig*cumsum(Y,1)));
                else
                    sig = sig(:)';
                    Y = bsxfun(@times,y0',exp(tspan*(mu(:)'-0.5*sig.^2)+bsxfun(@times,sig,cumsum(Y,1))));
                end
            end
        end
    else
        % Only allocate W matrix if requested as output
        if nargout == 2
            W = cumsum(Y,1);
            if N == 1
                Y = exp(tspan*mu+sig*W)*y0;
            else
                if isscalar(mu) && isscalar(sig)
                    Y = bsxfun(@times,y0',exp(bsxfun(@plus,tspan*mu,sig*W)));
                elseif isscalar(mu)
                    Y = bsxfun(@times,y0',exp(bsxfun(@plus,tspan*mu,bsxfun(@times,sig(:)',W))));
                elseif isscalar(sig)
                    Y = bsxfun(@times,y0',exp(tspan*mu(:)'+sig*W));
                else
                    Y = bsxfun(@times,y0',exp(tspan*mu(:)'+bsxfun(@times,sig(:)',W)));
                end
            end
        else
            if N == 1
                Y = exp(tspan*mu+sig*cumsum(Y,1))*y0;
            else
                if isscalar(mu) && isscalar(sig)
                    Y = bsxfun(@times,y0',exp(bsxfun(@plus,tspan*mu,sig*cumsum(Y,1))));
                elseif isscalar(mu)
                    Y = bsxfun(@times,y0',exp(bsxfun(@plus,tspan*mu,bsxfun(@times,sig(:)',cumsum(Y,1)))));
                elseif isscalar(sig)
                    Y = bsxfun(@times,y0',exp(tspan*mu(:)'+sig*cumsum(Y,1)));
                else
                    Y = bsxfun(@times,y0',exp(tspan*mu(:)'+bsxfun(@times,sig(:)',cumsum(Y,1))));
                end
            end
        end
    end
else
    % Only allocate W matrix if requested as output
    if nargout == 2
        W = zeros(lt,N,datatype);
    end
    
    % Solution not a function of sig, Ito and Stratonovich coincide
    if N == 1 || isscalar(mu)
        Y = exp(tspan*mu)*y0';
    else
        Y = bsxfun(@times,y0',exp(tspan*mu(:)'));
    end
end