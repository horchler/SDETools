function [Y,W,TE,YE,WE,IE] = sde_ou(th,mu,sig,tspan,y0,options)
%SDE_OU  Ornstein-Uhlenbeck mean-reverting process, analytic solution.
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
%   integrated scaled time-transformed Wiener increments that were used. Each
%   row of W corresponds to a time in TSPAN. W = W(exp(2*THETA*t)-1).
%
%   [...] = SDE_OU(THETA,MU,SIG,TSPAN,Y0,OPTIONS) returns as above with default
%   properties replaced by values in OPTIONS, an argument created with the
%   SDESET function. See SDESET for details. A commonly used option is to
%   manually specify the random seed via the RandSeed property, which creates a
%   new random number stream, instead of using the default stream, to generate
%   the Wiener increments.
%
%   [YOUT, W, TE, YE, WE, IE] = SDE_OU(THETA,MU,SIG,TSPAN,Y0,OPTIONS) with the
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
%       figure; h = plot(t([1 end]),[0 0],'k-.',t,y1,'b',t,y2,'r');
%       mustr = num2str(mu); npstr = num2str(npaths); dtstr = num2str(dt);
%       txt = {'Analytical solution',['Numerical solution, dt = ' dtstr]};
%       legend(h([2 end]),txt,'Location','NorthEast'); legend boxoff;
%       xlabel('t'); ylabel('y(t)');
%       title(['Ornstein-Uhlenbeck processes, ' npstr ' paths, \mu = ' mustr]);
%
%   Note:
%       The Ornstein-Uhlenbeck process is based on additive noise, i.e., the
%       diffusion term, g(t,y) = SIG, is not a function of the state variables.
%       In this case the Ito and Stratonovich interpretations are equivalent and
%       the SDEType OPTIONS property will have no effect.
%
%       Only diagonal noise is supported by this function. Setting the
%       DiagonalNoise OPTIONS property to 'no' to specify the more general
%       correlated noise case will result in an error. A numerical SDE solver
%       such as SDE_EULER should be used in this case or for other
%       generalizations, e.g., time-varying parameters.
%
%       In finance, the one-factor Hull-White and Vasicek (HMV) models are
%       examples of Ornstein-Uhlenbeck mean-reverting processes.
%
%   See also:
%       Explicit SDE solvers:	SDE_EULER, SDE_MILSTEIN
%       Implicit SDE solvers:   
%       Stochastic processes:	SDE_BM, SDE_GBM
%       Option handling:        SDESET, SDEGET, SDEPLOT
%       SDE demos/validation:   SDE_EULER_VALIDATE, SDE_MILSTEIN_VALIDATE
%   	Other:                  FUNCTION_HANDLE, RANDSTREAM

%   The conditional analytic solution used is for non-zero THETA
%       Y = Y0*exp(-THETA*t)+MU*(1-exp(-THETA*t))
%           +(SIG/sqrt(2*THETA))*exp(-THETA*t)*W(exp(2*THETA*t)-1),
%   where W() is a scaled time-transformed Wiener process. If THETA is zero, the
%   analytic solution for a driftless Wiener process is used: Y = Y0+SIG*W(t).
%
%   From: J. L. Doob, "The Brownian Movement and Stochastic Equations," Annals
%   of Mathematics, Vol. 43, No. 2, pp. 351-369, April 1942.

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-8-12
%   Revision: 1.2, 4-8-16


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
elseif nargin == 6
    if isempty(options) && (~sde_ismatrix(options) ...
            || any(size(options) ~= 0) || ~(isstruct(options) ...
            || iscell(options) || isnumeric(options))) ...
            || ~isempty(options) && ~isstruct(options)
        error('SDETools:sde_ou:InvalidSDESETStruct',...
              'Invalid SDE options structure.  See SDESET.');
    end
else
    error('SDETools:sde_ou:TooManyInputs',...
          'Too many input arguments.  See %s.',func);
end

% Check th, mu, and sig types
if isempty(th) || ~isfloat(th) || ~isvector(th)
    error('SDETools:sde_ou:ThetaEmptyOrNotFloatVector',...
         ['The drift rate parameter, THETA, must be non-empty '...
          'floating-point vector.  See %s.'],func);
end
if isempty(mu) || ~isfloat(mu) || ~isvector(mu)
    error('SDETools:sde_ou:MuEmptyOrNotFloatVector',...
         ['The drift mean parameter, MU, must be non-empty floating-point '...
          'vector.  See %s.'],func);
end
if isempty(sig) || ~isfloat(sig) || ~isvector(sig)
    error('SDETools:sde_ou:SigEmptyOrNotFloatVector',...
         ['The diffusion parameter, SIG, must be non-empty floating-point '...
          'vector.  See %s.'],func);
end

% Determine the dominant data type, single or double
dataType = superiorfloat(th,mu,sig,tspan,y0);
if ~all(strcmp(dataType,{class(th),class(mu),class(sig),class(tspan),...
        class(y0)}))
    warning('SDETools:sde_ou:InconsistentDataType',...
           ['Mixture of single and double data for inputs THETA, MU, SIG, '...
            'TSPAN, and Y0.']);
end

% Handle function arguments (NOTE: ResetStream is called by onCleanup())
[N,tspan,tdir,lt,y0,h,ConstStep,Stratonovich,RandFUN,ResetStream,EventsFUN,...
    EventsValue,OutputFUN,WSelect] ...
	= sdearguments_process(func,tspan,y0,dataType,options);	%#ok<ASGLU>

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
    if nargout > 6
        error('SDETools:sde_ou:EventsTooManyOutputs',...
              'Too many output arguments.  See %s.',func);
    else
        if nargout >= 3
            TE = [];
            if nargout >= 4
                YE = [];
                if nargout >= 5
                    WE = [];
                    if nargout == 6
                        IE = [];
                    end
                end
            end
        end
    end
else
    if nargout > 2
        if nargout <= 6
            error('SDETools:sde_ou:NoEventsTooManyOutputs',...
                 ['Too many output arguments. An events function has not '...
                  'been specified.  See %s.'],func);
        else
            error('SDETools:sde_ou:TooManyOutputs',...
                  'Too many output arguments.  See %s.',func);
        end
    end
end

% Initialize output function
isOutput = ~isempty(OutputFUN);
isW = (nargout >= 2 || isEvents || WSelect);

% State array
Y(lt,N) = cast(0,dataType);

% Expand and orient parameter and y0 vectors, find non-zero values
if N > 1
    if isscalar(sig)
        sig = zeros(1,N)+sig;
    else
        sig = sig(:).';
    end
    th = th(:).';
    mu = mu(:).';
    y0 = y0.';
end
th0 = (th ~= 0);

% Check if alternative RandFUN function or W matrix is present
CustomRandFUN = (isempty(RandFUN) && isfield(options,'RandFUN'));
CustomWMatrix = (CustomRandFUN && ~isa(options.RandFUN,'function_handle'));

% Reduce effective dimension of noise if integrated Wiener increments not needed
if isW || CustomWMatrix
    sig0 = true(1,N);
    D = N;
else
    sig0 = (sig ~= 0);
    D = nnz(sig0);
end

if CustomRandFUN
    if CustomWMatrix
        % User-specified integrated Wiener increments
        Y = sdeget(options,'RandFUN',[],'flag');
        if ~isfloat(Y) || ~sde_ismatrix(Y) || any(size(Y) ~= [lt D])
            error('SDETools:sde_ou:RandFUNInvalidW',...
                 ['RandFUN must be a function handle or a '...
                  'LENGTH(TSPAN)-by-D (%d by %d) floating-point matrix of '...
                  'integrated Wiener increments.  See %s.'],lt,D,func);
        end
    else
        % User-specified function handle
        RandFUN = sdeget(options,'RandFUN',[],'flag');    

        try
            % Store Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,D);
            if ~sde_ismatrix(r) || isempty(r) || ~isfloat(r)
                error('SDETools:sde_ou:RandFUNNot2DArray3',...
                     ['RandFUN must return a non-empty matrix of '...
                      'floating-point values.  See %s.'],func);
            end
            [m,n] = size(r);
            if m ~= lt-1 || n ~= D
                error('SDETools:sde_ou:RandFUNDimensionMismatch3',...
                     ['The specified alternative RandFUN did not output '...
                      'a %d by %d matrix as requested.   See %s.'],lt-1,D,func);
            end
            
            Y(2:end,sig0) = r;
            clear r;            % Remove large temporary variable to save memory
        catch err
            switch err.identifier
                case 'MATLAB:TooManyInputs'
                    error('SDETools:sde_ou:RandFUNTooFewInputs',...
                          'RandFUN must have at least two inputs.  See %s.',...
                          func);
                case 'MATLAB:TooManyOutputs'
                    error('SDETools:sde_ou:RandFUNNoOutput',...
                         ['The output of RandFUN was not specified. RandFUN '...
                          'must return a non-empty matrix.  See %s.'],func);
                case 'MATLAB:unassignedOutputs'
                    error('SDETools:sde_ou:RandFUNUnassignedOutput',...
                         ['The first output of RandFUN was not assigned.  '...
                          'See %s.'],func);
                case 'MATLAB:minrhs'
                    error('SDETools:sde_ou:RandFUNTooManyInputs',...
                         ['RandFUN must not require more than two inputs.  '...
                          'See %s.'],func);
                otherwise
                    rethrow(err);
            end
        end
    end
else
    % Store random variates in Y, no error checking needed for default RANDN
    Y(2:end,sig0) = feval(RandFUN,lt-1,D);
end

% Diffusion parameters are not all zero
if D > 0
    % Store scaled time-transformed Wiener increments in Y
    if all(th0)
     	tt = -tspan*th;
        if N == 1 || ~isscalar(th)
            Y(2:end,sig0) = tdir*sqrt(diff(expm1(-2*tt),1,1)).*Y(2:end,sig0);
        else
            Y(2:end,sig0) = bsxfun(@times,tdir*sqrt(diff(expm1(-2*tt),1,1)),Y(2:end,sig0));
        end
	elseif all(~th0)
        if N == 1 || ConstStep
            Y(2:end,sig0) = tdir*sqrt(h).*Y(2:end,sig0);
        else
            Y(2:end,sig0) = bsxfun(@times,tdir*sqrt(h),Y(2:end,sig0));
        end
    else
        i = th0 & sig0;
        D = nnz(i);
        
        tt = -tspan*th(i);
        if D == 1
            Y(2:end,i) = tdir*sqrt(diff(expm1(-2*tt),1,1)).*Y(2:end,i);
        else
            Y(2:end,i) = bsxfun(@times,tdir*sqrt(diff(expm1(-2*tt),1,1)),Y(2:end,i));
        end
        
        i = ~th0 & sig0;
        if N-D == 1 || ConstStep
            Y(2:end,i) = tdir*sqrt(h).*Y(2:end,i);
        else
            Y(2:end,i) = bsxfun(@times,tdir*sqrt(h),Y(2:end,i));
        end
    end
    
    % Integrate Wiener increments
    if ~CustomWMatrix
        Y(2:end,sig0) = cumsum(Y(2:end,sig0),1);
    end
    
    % Only allocate W matrix if requested as output or needed
    if isW
        W = Y;
    end
    
    % Evaluate analytic solution
    if all(th0)
        % All th ~= 0
        ett = exp(tt);
        if N == 1
            Y = mu+ett.*(y0-mu+(sig/sqrt(2*th))*Y);
        elseif isscalar(th)
            if isscalar(mu)
                Y = ett*(y0-mu)+mu+ett*(sig/sqrt(2*th)).*Y;
            else
                Y = ett*y0-expm1(tt)*mu+ett*(sig/sqrt(2*th)).*Y;
            end
        else
            if isscalar(mu)
                Y = mu+ett.*(bsxfun(@plus,y0-mu,bsxfun(@times,sig./sqrt(2*th),Y)));
            else
                Y = ett.*(bsxfun(@plus,y0,bsxfun(@times,sig./sqrt(2*th),Y)))-bsxfun(@times,expm1(tt),mu);
            end
        end
	elseif all(~th0)
        % All th = 0, driftless noise
        if N == 1
            Y = y0+sig*Y;
        else
            if isscalar(sig)
                Y = bsxfun(@plus,y0,sig*Y);
            else
                Y = bsxfun(@plus,y0,bsxfun(@times,sig,Y));
            end
        end
    else
        % Some th ~= 0
        i = th0;
        D = nnz(i);
        th = th(i);
        if isscalar(mu)
            mu = mu(ones(1,N));
        else
            mu = mu(:).';
        end
        if D == 1
            Y(:,i) = exp(-tspan*th).*(y0(i)-mu(i)+(sig(i)/sqrt(2*th))*Y(:,i));
        else
            Y(:,i) = exp(-tspan*th).*(bsxfun(@plus,y0(i)-mu(i),bsxfun(@times,sig(i)./sqrt(2*th),Y(:,i))));
        end
        
        % Some th = 0, driftless noise
        i = ~i;
        if N-D == 1 && isscalar(sig)
            Y(:,i) = y0(i)+sig*Y(:,i);
        else
            if isscalar(sig)
                Y(:,i) = bsxfun(@plus,y0(i),sig*Y(:,i));
            else
                Y(:,i) = bsxfun(@plus,y0(i),bsxfun(@times,sig(i),Y(:,i)));
            end
        end
    end
else
    % Solution not a function of sig
    if all(th0)
        % All th ~= 0, pure drift, noise magnitude, sig, is zero
        tt = -tspan*th;
        ett = exp(tt);
        if N == 1
            Y = ett.*(y0-mu);
        elseif isscalar(th)
            if isscalar(mu)
                Y = ett*(y0-mu)+mu;
            else
                Y = ett*y0-expm1(tt)*mu;
            end
        else
            Y = bsxfun(@times,ett,y0-mu);
        end
    elseif all(~th0)
        % All th = 0, driftless noise, but noise magnitude, sig, is zero
        if N == 1
            Y = Y+y0;
        else
            Y = bsxfun(@plus,y0,Y);
        end
    else
        % Some th ~= 0, pure drift, noise magnitude, sig, is zero
        i = th0;
        if isscalar(mu)
            mu = mu(ones(1,N));
        else
            mu = mu(:).';
        end
        if nnz(i) == 1
            Y(:,i) = exp(-tspan*th(i)).*(y0(i)-mu(i));
        else
            Y(:,i) = bsxfun(@times,exp(-tspan*th(i)),y0(i)-mu(i));
        end
        
        % Some th = 0, driftless noise, but noise magnitude, sig, is zero
        i = ~i;
        Y(:,i) = bsxfun(@plus,y0(i),Y(:,i));
    end
end

% Check for and handle zero-crossing events, and output function
if isEvents
    for i = 2:lt
        [te,ye,we,ie,EventsValue,IsTerminal] ...
            = sdezero(EventsFUN,tspan(i),Y(i,:).',W(i,:).',EventsValue);
        if ~isempty(te)
            if nargout >= 3
                TE = [TE;te];               %#ok<AGROW>
                if nargout >= 4
                    YE = [YE;ye];           %#ok<AGROW>
                    if nargout >= 5
                        WE = [WE;we];       %#ok<AGROW>
                        if nargout == 6
                            IE = [IE;ie];	%#ok<AGROW>
                        end
                    end
                end
            end
            if IsTerminal
                Y = Y(1:i,:);
                if nargout >= 2
                    W = W(1:i,:);
                end
                break;
            end
        end
        
        if isOutput
            OutputFUN(tspan(i),Y(i,:).','',W(i,:).');
        end
    end
elseif isOutput
    if isW
        for i = 2:lt
            OutputFUN(tspan(i),Y(i,:).','',W(i,:).');
        end
    else
        for i = 2:lt
            OutputFUN(tspan(i),Y(i,:).','',[]);
        end
    end
end

% Finalize output
if isOutput
    OutputFUN([],[],'done',[]);
end