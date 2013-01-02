function [Y,W,TE,YE,IE] = sde_ou(th,mu,sig,tspan,y0,options,varargin)
%SDE_OU  Ornstein-Uhlenbeck process, analytic solution.
%   YOUT = SDE_OU(THETA,MU,SIG,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL] returns
%   the analytic solution of the N-dimensional system of stochastic differential
%   equations for the Ornstein-Uhlenbeck process, dY = THETA*(MU-Y)*dt + SIG*dW,
%   with N-dimensional diagonal noise from time T0 to TFINAL (all increasing or 
%   all decreasing with arbitrary step size) with initial conditions Y0. TSPAN
%   is a length M vector. Y0 is a length N vector. The drift rate parameter
%   THETA, the drift mean parameter MU, and the diffusion parameter SIG may be
%   scalars or length N vectors. Each row in the M-by-N solution array YOUT
%   corresponds to a time in TSPAN.
%
%   [YOUT, W] = SDE_OU(THETA,MU,SIG,TSPAN,Y0,...) outputs the M-by-N matrix W of
%   integrated Weiner increments that were used by the solver. Each row of W
%   corresponds to a time in TSPAN.
%
%   [...] = SDE_OU(THETA,MU,SIG,TSPAN,Y0,OPTIONS) returns as above with default
%   properties replaced by values in OPTIONS, an argument created with the
%   SDESET function. See SDESET for details. A commonly used option is to
%   manually specify the random seed via the RandSeed property, which creates a
%   new random number stream, instead of using the default stream, to generate
%   the Wiener increments. If the DiagonalNoise property is set to 'no' to
%   specify the more general correlated noise case, SIG may be an N-by-D matrix,
%   and a numerical solution will be returned.
%
%   [YOUT, W, TE, YE, IE] = SDE_OU(THETA,MU,SIG,TSPAN,Y0,OPTIONS) with the
%   EventsFUN property set to a function handle, in order to specify an events
%   function, solves as above while also finding zero-crossings. The
%   corresponding function, must take at least two inputs and output three 
%   vectors: [Value, IsTerminal, Direction] = EventsFUN(T,Y). The scalar input T
%   is the current integration time and the vector Y is the current state. For
%   the i-th event, Value(i) is the value of the zero-crossing function and
%   IsTerminal(i) = 1 specifies that integration is to terminate at a zero or to
%   continue if IsTerminal(i) = 0. If Direction(i) = 1, only zeros where
%   Value(i) is increasing are found, if Direction(i) = -1, only zeros where
%   Value(i) is decreasing are found, otherwise if Direction(i) = 0, all zeros
%   are found. If Direction is set to the empty matrix, [], all zeros are found
%   for all events. Direction and IsTerminal may also be scalars.
%
%   Example:
%       % Compare analytical and simulated Ornstein-Uhlenbeck processes
%       npaths = 10; dt = 1e-2; t = 0:dt:1; y0 = -1:2/(npaths-1):1;
%       th = 4; mu = 0; sig = 0.25; opts = sdeset('RandSeed',1);
%       y1 = sde_ou(th,mu,sig,t,y0,opts);
%       y2 = sde_euler(@(t,y)th.*(mu-y),sig,t,y0,opts);
%       h = plot(t([1 end]),[0 0],'k-.',t,y1,'b',t,y2,'r');
%       mustr = num2str(mu); npstr = num2str(npaths); dtstr = num2str(dt);
%       txt = {'Analytical solution',['Numerical solution, dt = ' dtstr]};
%       legend(h([2 end]),txt,1); legend boxoff; xlabel('t'); ylabel('y(t)');
%       title(['Ornstein-Uhlenbeck processes, ' npstr ' paths, \mu = ' mustr]);
%
%   Note:
%       The Ornstein-Uhlenbeck process is based on additive noise, i.e., the
%       diffusion term, g(t,y) = SIG, is not a function of the state variables.
%       In this case the Ito and Stratonovich interpretations are equivalent.
%
%   See also:
%       Explicit SDE solvers:	SDE_EULER, SDE_MILSTEIN
%       Implicit SDE solvers:   
%       Stochastic processes:	SDE_GBM
%       Option handling:        SDESET, SDEGET
%       SDE demos/validation:   SDE_EULER_VALIDATE, SDE_MILSTEIN_VALIDATE
%   	Other:                  FUNCTION_HANDLE, RANDSTREAM

%   The conditional analytic solution used is
%       Y = Y0*exp(-THETA*t)+MU*(1-exp(-THETA*t))
%           +(SIG/sqrt(2*THETA))*exp(-THETA*t)*W(exp(2*THETA*t)-1),
%   where W() is a scaled time-transformed Wiener process.
%
%   From: J. L. Doob, "The Brownian Movement and Stochastic Equations," Annals
%   of Mathematics, Vol. 43, No. 2, pp. 351-369, April 1942.

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-8-12
%   Revision: 1.0, 1-1-13


func = 'SDE_OU';

% Check inputs and outputs
if nargin < 6
    if nargin < 5
        error('SDETools:sde_ou:NotEnoughInputs',...
              'Not enough input arguments.  See %s.',func);
    end
    if isa(y0,'struct')
        error('SDETools:sde_ou:NotEnoughInputsOptions',...
             ['An SDE options structure was provided as the last argument, '...
              'but one of the first four input arguments is missing.'...
              '  See %s.'],func);
    end
    options = [];
elseif isempty(options) && (~sde_ismatrix(options) ...
        || any(size(options) ~= 0) || ~(isstruct(options) || iscell(options) ...
        || isnumeric(options))) || ~isempty(options) && ~isstruct(options)
	error('SDETools:sde_ou:InvalidSDESETStruct',...
          'Invalid SDE options structure.  See SDESET.');
end

% Check th, mu, and sig types
if isempty(th) || ~isfloat(th) || ~isvector(th)
    error('SDETools:sde_ou:ThetaEmptyOrNotFloatVector',...
         ['The drift rate parameter, THETA, must be non-empty vector of '...
          'singles or doubles.  See %s.'],func);
end
if isempty(mu) || ~isfloat(mu) || ~isvector(mu)
    error('SDETools:sde_ou:MuEmptyOrNotFloatVector',...
         ['The drift mean parameter, MU, must be non-empty vector of '...
          'singles or doubles.  See %s.'],func);
end
if isempty(sig) || ~isfloat(sig) || ~isvector(sig)
    error('SDETools:sde_ou:SigEmptyOrNotFloatVector',...
         ['The diffusion parameter, SIG, must be non-empty vector of '...
          'singles or doubles.  See %s.'],func);
end

% Determine the dominant data type, single or double
dataType = superiorfloat(th,mu,sig,tspan,y0);
if ~all(strcmp(dataType,{class(th),class(mu),class(sig),class(tspan),...
        class(y0)}))
    warning('SDETools:sde_ou:InconsistentDataType',...
           ['Mixture of single and double data for inputs THETA, MU, SIG, '...
            'TSPAN, and Y0.']);
end
isDouble = strcmp(dataType,'double');

% Handle function arguments (NOTE: ResetStream is called by onCleanup())
[N,tspan,tdir,lt,y0,h,ConstStep,Stratonovich,RandFUN,CustomRandFUN,...
    ResetStream,EventsFUN,EventsValue]...
	= sdearguments_special(func,tspan,y0,dataType,options,varargin);	%#ok<ASGLU>

% Check th, mu, and sig sizes
if ~any(length(th) == [1 N])
    error('SDETools:sde_ou:ThetaDimensionMismatch',...
         ['The drift rate parameter, THETA, must be a scalar or a vector '...
          'the same length as Y0.  See %s.'],func);
end
if ~any(length(mu) == [1 N])
    error('SDETools:sde_ou:MuDimensionMismatch',...
         ['The drift mean parameter, MU, must be a scalar or a vector the '...
          'same length as Y0.  See %s.'],func);
end
if ~any(length(sig) == [1 N])
    error('SDETools:sde_ou:SigDimensionMismatch',...
         ['The diffusion parameter, SIG, must be a scalar or a vector the '...
          'same length as Y0.  See %s.'],func);
end

% Check sign of sig
if any(sig < 0)
    error('SDETools:sde_ou:SigNegative',...
         ['The diffusion parameter, SIG, must be greater than or equal to '...
          'zero.  See %s.'],func);
end

% Initialize outputs for zero-crossing events
isEvents = ~isempty(EventsFUN);
if isEvents
    if nargout > 5
        error('SDETools:sde_ou:EventsTooManyOutputs',...
              'Too many output arguments.  See %s.',func);
    else
        if nargout >= 3
            TE = [];
            if nargout >= 4
                YE = [];
                if nargout >= 5
                    IE = [];
                end
            end
        end
    end
else
    if nargout > 2
        if nargout <= 5
            error('SDETools:sde_ou:NoEventsTooManyOutputs',...
                 ['Too many output arguments. An events function has not '...
                  'been specified.  See %s.'],func);
        else
            error('SDETools:sde_ou:TooManyOutputs',...
                  'Too many output arguments.  See %s.',func);
        end
    end
end

% State array
if isDouble
    Y(lt,N) = 0;
else
    Y(lt,N) = single(0);
end

% Diffusion parameters aren't all zero
if ~all(sig == 0)
    % Wiener increments from normal variates, store in state if possible
    if CustomRandFUN    % check output of alternative RandFUN
        try
            % Store scaled time-transformed Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,N);
            if ~sde_ismatrix(r) || isempty(r) || ~isfloat(r)
                error('SDETools:sde_ou:RandFUNNot2DArray3',...
                     ['RandFUN must return a non-empty matrix of floating '...
                      'point values.  See %s.'],solver);
            end
            [m,n] = size(r);
            if m ~= lt-1 || n ~= N
                error('SDETools:sde_ou:RandFUNDimensionMismatch3',...
                     ['The specified alternative RandFUN did not output a '...
                      '%d by %d matrix as requested.   See %s.',N,lt-1,solver]);
            end
            if N == 1 || isscalar(th)
                tt = -tspan*th;
                Y(2:end,:) = tdir*sqrt(diff(expm1(-2*tt),1,1)).*r;
            else
                th = th(:)';
                tt = -tspan*th;
                Y(2:end,:) = bsxfun(@times,tdir*sqrt(diff(expm1(-2*tt),1,1)),r);
            end
            clear r;    % remove large temporary variable to save memory
        catch err
            switch err.identifier
                case 'MATLAB:TooManyInputs'
                    error('SDETools:sde_ou:RandFUNTooFewInputs',...
                          'RandFUN must have at least two inputs.  See %s.',...
                          solver);
                case 'MATLAB:TooManyOutputs'
                    error('SDETools:sde_ou:RandFUNNoOutput',...
                         ['The output of RandFUN was not specified. RandFUN '...
                          'must return a non-empty matrix.  See %s.'],solver);
                case 'MATLAB:unassignedOutputs'
                    error('SDETools:sde_ou:RandFUNUnassignedOutput',...
                         ['The first output of RandFUN was not assigned.'...
                          '  See %s.'],solver);
                case 'MATLAB:minrhs'
                    error('SDETools:sde_ou:RandFUNTooManyInputs',...
                         ['RandFUN must not require more than two inputs.'...
                          '  See %s.'],solver);
                otherwise
                    rethrow(err);
            end
        end
    else    % No error checking needed if default RANDN used
        % Store scaled time-transformed Wiener increments in Y
        if N == 1 || isscalar(th)
            tt = -tspan*th;
            Y(2:end,:) = tdir*sqrt(diff(expm1(-2*tt),1,1)).*feval(RandFUN,lt-1,N);
        else
            th = th(:)';
            tt = -tspan*th;
            Y(2:end,:) = bsxfun(@times,tdir*sqrt(diff(expm1(-2*tt),1,1)),feval(RandFUN,lt-1,N));
        end
    end
    
    % Only allocate W matrix if requested as output
    ett = exp(tt);
    if nargout >= 2
        W = cumsum(Y,1);
        if N == 1 || isscalar(th) && ~isscalar(mu) && ~isscalar(sig)
            Y = ett*y0'-expm1(tt)*mu(:)'+ett*(sig(:)'/sqrt(2*th)).*W;
        else
            if isscalar(th) && isscalar(mu) && isscalar(sig)
                Y = bsxfun(@minus,ett*y0',expm1(tt)*mu)+bsxfun(@times,ett*(sig/sqrt(2*th)),W);
            elseif isscalar(th) && isscalar(mu)
                Y = bsxfun(@minus,ett*y0',expm1(tt)*mu)+ett*(sig(:)'/sqrt(2*th)).*W;
            elseif isscalar(th) && isscalar(sig)
                Y = ett*y0'-expm1(tt)*mu(:)'+bsxfun(@times,ett*(sig/sqrt(2*th)),W);
            elseif isscalar(mu)
                Y = bsxfun(@times,ett,y0')-expm1(tt)*mu+bsxfun(@times,ett,(sig(:)'./sqrt(2*th))).*W;
            else
                Y = bsxfun(@times,ett,y0')-bsxfun(@times,expm1(tt),mu(:)')+bsxfun(@times,ett,(sig(:)'./sqrt(2*th))).*W;
            end
        end
    else
        if N == 1 || isscalar(th) && ~isscalar(mu) && ~isscalar(sig)
            Y = ett*y0'-expm1(tt)*mu(:)'+ett*(sig(:)'/sqrt(2*th)).*cumsum(Y,1);
        else
            if isscalar(th) && isscalar(mu) && isscalar(sig)
                Y = bsxfun(@minus,ett*y0',expm1(tt)*mu)+bsxfun(@times,ett*(sig/sqrt(2*th)),cumsum(Y,1));
            elseif isscalar(th) && isscalar(mu)
                Y = bsxfun(@minus,ett*y0',expm1(tt)*mu)+ett*(sig(:)'/sqrt(2*th)).*cumsum(Y,1);
            elseif isscalar(th) && isscalar(sig)
                Y = ett*y0'-expm1(tt)*mu(:)'+bsxfun(@times,ett*(sig/sqrt(2*th)),cumsum(Y,1));
            elseif isscalar(mu)
                Y = bsxfun(@times,ett,y0')-expm1(tt)*mu+bsxfun(@times,ett,(sig(:)'./sqrt(2*th))).*cumsum(Y,1);
            else
                Y = bsxfun(@times,ett,y0')-bsxfun(@times,expm1(tt),mu(:)')+bsxfun(@times,ett,(sig(:)'./sqrt(2*th))).*cumsum(Y,1);
            end
        end
    end
else
    % Only allocate W matrix if requested as output
    if nargout >= 2
        if isDouble
            W(lt,N) = 0;
        else
            W(lt,N) = single(0);
        end
    end
    
    % Solution not a function of sig
    if N == 1 || isscalar(th) && ~isscalar(mu)
        tt = -tspan*th;
        Y = exp(tt)*y0'-expm1(tt)*mu(:)';
    else
        if isscalar(th) && isscalar(mu)
            tt = -tspan*th;
            Y = bsxfun(@minus,exp(tt)*y0',expm1(tt)*mu);
        elseif isscalar(mu)
            tt = -tspan*th(:)';
            Y = bsxfun(@minus,bsxfun(@times,y0',exp(tt)),expm1(tt)*mu);
        else
            tt = -tspan*th(:)';
            Y = bsxfun(@times,y0',exp(tt))-bsxfun(@times,expm1(tt),mu(:)');
        end
    end
end

% Check for and handle zero-crossing events
if isEvents
    for i = 2:lt
        [te,ye,ie,EventsValue,IsTerminal] = sdezero(EventsFUN,tspan(i),Y(i,:),EventsValue,varargin);
        if ~isempty(te)
            if nargout >= 3
                TE = [TE;te];           %#ok<AGROW>
                if nargout >= 4
                    YE = [YE;ye];       %#ok<AGROW>
                    if nargout >= 5
                        IE = [IE;ie];	%#ok<AGROW>
                    end
                end
            end
            if IsTerminal
                Y = Y(1:i,:);
                if nargout >= 2
                    W = W(1:i,:);
                end
                return;
            end
        end
    end
end