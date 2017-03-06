function [Y,W,TE,YE,WE,IE] = sde_euler(f,g,tspan,y0,options)
%SDE_EULER  Solve stochastic differential equations, Euler methods.
%   YOUT = SDE_EULER(FFUN,GFUN,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL]
%   integrates the N-dimensional system of stochastic differential equations
%   dy = f(t,y)*dt + g(t,y)*dW with N-dimensional diagonal noise from time T0 to
%   TFINAL (all increasing or all decreasing with arbitrary step size) with
%   initial conditions Y0. TSPAN is a length M vector. Y0 is a length N vector.
%   Each row in the M-by-N solution array YOUT corresponds to a time in TSPAN.
%   FFUN and GFUN are function handles or floating-point matrices. For scalar T
%   and vector Y, FFUN(T,Y) and GFUN(T,Y) return length N column vectors
%   corresponding to f(t,y) and g(t,y), the deterministic and stochastic parts
%   of the SDE, respectively. GFUN may also return a scalar to be used across
%   all N dimensions. FFUN may also be a vector the same length as Y0, a scalar,
%   or the empty matrix, []. GFUN may also be a matrix, a scalar, or the empty
%   matrix, []. If FFUN or GFUN is not a function handle, the corresponding
%   function (and its derivative) are assumed constant. If FFUN or GFUN is a
%   vector or matrix of all zeros or the empty matrix, [], the corresponding
%   function is neglected entirely. Each row in the M-by-N solution array YOUT
%   corresponds to a time in TSPAN.
%
%   [YOUT, W] = SDE_EULER(FFUN,GFUN,TSPAN,Y0,...) outputs the M-by-N matrix W of
%   integrated Wiener increments that were used by the solver. Each row of W
%   corresponds to a time in TSPAN.
%
%   [...] = SDE_EULER(FFUN,GFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the SDESET function. See SDESET for details. The type of SDE to be
%   integrated, 'Ito' or the default 'Stratonovich', can be specified via the
%   SDEType property. GFUN can output more general D-dimensional correlated
%   noise by setting the DiagonalNoise property to 'no'. The Euler-Maruyama
%   method is used for Ito SDEs. The Euler-Heun method is used for Stratonovitch
%   SDEs unless the ConstGFUN property is set to 'yes', in which case the Ito
%   and Stratonovitch SDEs are equivalent. The DiagonalNoise property must be
%   set to 'no' in order to apply the more general correlated noise case where
%   GFUN may return (or is) an N-by-D matrix, where D is the dimension of the
%   noise. Another commonly used option is to manually specify the random seed
%   via the RandSeed property, which creates a new random number stream, instead
%   of using the default stream, to generate the Wiener increments.
%
%   [YOUT, W, TE, YE, WE, IE] = SDE_EULER(FFUN,GFUN,TSPAN,Y0,OPTIONS) with the
%   EventsFUN property set to a function handle, in order to specify an events
%   function, solves as above while also finding zero-crossings. The
%   corresponding function must take at least two inputs and output three
%   vectors: [Value, IsTerminal, Direction] = Events(T,Y). The scalar input T is
%   the current integration time and the vector Y is the current state. For the
%   i-th event, Value(i) is the value of the zero-crossing function and
%   IsTerminal(i) = 1 specifies that integration is to terminate at a zero or to
%   continue if IsTerminal(i) = 0. If Direction(i) = 1, only zeros where
%   Value(i) is increasing are found, if Direction(i) = -1, only zeros where
%   Value(i) is decreasing are found, otherwise if Direction(i) = 0, all zeros
%   are found. If Direction is set to the empty matrix, [], all zeros are found
%   for all events. Direction and IsTerminal may also be scalars.
%
%   Example:
%       % Solve 2-D Stratonovich SDE using Euler-Heun method
%       mu = 1; sig = [0.1;0.5]; dt = 1e-2; t = 0:dt:1;
%       f = @(t,y)mu.*y; g = @(t,y)sig.*y; opts = sdeset('RandSeed',1);
%       y = sde_euler(f,g,t,[1 1],opts);
%       figure; plot(t,y);
%       title(['Euler-Heun Method, dt = ' num2str(dt) ', \mu = ' num2str(mu)]);
%       txt = {['\sigma = ' num2str(sig(1))],['\sigma = ' num2str(sig(2))]};
%       legend(txt,'Location','NorthEast'); legend boxoff;
%       xlabel('t'); ylabel('y(t)');
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
%       Stochastic processes:	SDE_BM, SDE_GBM, SDE_OU
%       Option handling:        SDESET, SDEGET, SDEPLOT
%       SDE demos/validation:   SDE_EULER_VALIDATE, SDE_MILSTEIN_VALIDATE
%   	Other:                  FUNCTION_HANDLE, RANDSTREAM

%   SDE_EULER is an implementation of the order 0.5 strong (order 1.0 weak)
%   explicit Euler-Maruyama and Euler-Heun schemes, which are also order 0.5
%   strong Taylor approximation schemes. In the additive noise case (where the
%   diffusion is not a function of the state, but may be a function of time),
%   Euler schemes improve to order 1.0 strong convergence.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, horchler @ gmail . com, Created 10-28-10
%   Revision: 1.3, 6-8-16


solver = 'SDE_EULER';

% Check inputs and outputs
if nargin < 5
    if nargin < 4
        error('SDETools:sde_euler:NotEnoughInputs',...
              'Not enough input arguments.  See %s.',solver);
    end
    if isa(y0,'struct')
        error('SDETools:sde_euler:NotEnoughInputsOptions',...
             ['An SDE options structure was provided as the last argument, '...
              'but one of the first four input arguments is missing.'...
              '  See %s.'],solver);
    end
    options = [];
elseif nargin == 5
    if isempty(options) && (~sde_ismatrix(options) ...
            || any(size(options) ~= 0) || ~(isstruct(options) ...
            || iscell(options) || isnumeric(options))) ...
            || ~isempty(options) && ~isstruct(options)
        error('SDETools:sde_euler:InvalidSDESETStruct',...
              'Invalid SDE options structure.  See SDESET.');
    end
else
    error('SDETools:sde_euler:TooManyInputs',...
          'Too many input arguments.  See %s.',solver);
end

% Handle solver arguments (NOTE: ResetStream is called by onCleanup())
[N,D,tspan,tdir,lt,y0,fout,gout,dgout,dg,h,ConstStep,dataType,NonNegative,...
    idxNonNegative,DiagonalNoise,ScalarNoise,OneDNoise,ConstFFUN,ConstGFUN,...
    ConstDGFUN,Stratonovich,RandFUN,ResetStream,EventsFUN,EventsValue,...
    OutputFUN,WSelect] = sdearguments(solver,f,g,tspan,y0,options);	%#ok<ASGLU>

% Initialize outputs for zero-crossing events
isEvents = ~isempty(EventsFUN);
if isEvents
    if nargout > 6
        error('SDETools:sde_euler:EventsTooManyOutputs',...
              'Too many output arguments.  See %s.',solver);
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
            error('SDETools:sde_euler:NoEventsTooManyOutputs',...
                 ['Too many output arguments. An events function has not '...
                  'been specified.  See %s.'],solver);
        else
            error('SDETools:sde_euler:TooManyOutputs',...
                  'Too many output arguments.  See %s.',solver);
        end
    end
end

% Initialize output function
isOutput = ~isempty(OutputFUN);

% If drift or diffusion functions, FFUN and GFUN, exist
isDrift = ~(ConstFFUN && isscalar(fout) && fout == 0);
isDiffusion = ~(ConstGFUN && isscalar(gout) && gout == 0);

isW = isDiffusion && (isEvents || WSelect);
if isW
    Wi = 0;
else
    Wi = [];
end

% Is Y allocated and output
isYOutput = (nargout > 0);

% Check if alternative RandFUN function or W matrix is present
if isempty(RandFUN) && isfield(options,'RandFUN')
    CustomRandFUN = isa(options.RandFUN,'function_handle');
    CustomWMatrix = ~CustomRandFUN;
else
    CustomRandFUN = false;
    CustomWMatrix = false;
end

% Location of stored Wiener increments, if they are needed and pre-calculated
dWinY = (D <= N && isYOutput && ~CustomWMatrix || ~isDiffusion);  	% Store in Y
dWinW = (isDiffusion && D > N && nargout >= 2 || CustomWMatrix);	% Store in W

% Allocate state array, Y, if needed (may be allocated in place below)
if isYOutput && (~(CustomRandFUN && dWinY) ...
        || ~(ConstFFUN && ConstGFUN && dWinW && ~NonNegative))
    Y(lt,N) = cast(0,dataType);
end
    
% Calculate Wiener increments from normal variates, store in Y if possible, or W
sh = tdir*sqrt(h);
h = tdir*h;
if isDiffusion || isW
    if CustomRandFUN                            % Check alternative RandFUN
        if CustomWMatrix
            W = sdeget(options,'RandFUN',[],'flag');
            if ~isfloat(W) || ~sde_ismatrix(W) || any(size(W) ~= [lt D])
                error('SDETools:sde_euler:RandFUNInvalidW',...
                     ['RandFUN must be a function handle or a '...
                      'LENGTH(TSPAN)-by-D (%d by %d) floating-point matrix '...
                      'of integrated Wiener increments.  See %s.'],lt,D,solver);
            end
            error('SDETools:sde_euler:RandFUNWMatrixNotSupportedYet',...
                 ['Custom W matrices specified via RandFUN not supported '...
                  'yet in SDE_EULER and SDE_MILSTEIN.']);
        else
            % User-specified function handle
            RandFUN = sdeget(options,'RandFUN',[],'flag');
            
            try
                if dWinY                        % Store Wiener increments in Y
                    Y = feval(RandFUN,lt-1,D);
                    if ~sde_ismatrix(Y) || isempty(Y) || ~isfloat(Y)
                        error('SDETools:sde_euler:RandFUNNot2DArray3',...
                             ['RandFUN must return a non-empty matrix of '...
                              'floating-point values.  See %s.'],solver);
                    end
                    [m,n] = size(Y);
                    if m ~= lt-1 || n ~= D
                        error('SDETools:sde_euler:RandFUNDimensionMismatch3',...
                             ['The specified alternative RandFUN did not '...
                              'output a %d by %d matrix as requested.  '...
                              'See %s.'],lt-1,D,solver);
                    end

                    if ScalarNoise || ConstStep
                        Y = [zeros(1,D,dataType);
                             sh.*Y zeros(lt-1,N-D,dataType)];
                    else
                        Y = [zeros(1,D,dataType);
                             bsxfun(@times,sh,Y) zeros(lt-1,N-D,dataType)];
                    end
                    if nargout >= 2
                        W = cumsum(Y(:,1:D),1); % Integrated Wiener increments
                    end
                elseif dWinW                    % Store Wiener increments in W
                    W = feval(RandFUN,lt-1,D);
                    if ~sde_ismatrix(W) || isempty(W) || ~isfloat(W)
                        error('SDETools:sde_euler:RandFUNNot2DArray1',...
                             ['RandFUN must return a non-empty matrix of '...
                              'floating-point values.  See %s.'],solver);
                    end
                    [m,n] = size(W);
                    if m ~= lt-1 || n ~= D
                        error('SDETools:sde_euler:RandFUNDimensionMismatch1',...
                             ['The specified alternative RandFUN did not '...
                              'output a %d by %d matrix as requested.  '...
                              'See %s.'],lt-1,D,solver);
                    end

                    if ConstStep
                        W = [zeros(1,D,dataType);sh.*W];
                    else
                        W = [zeros(1,D,dataType);bsxfun(@times,sh,W)];
                    end
                else                            % Cannot store Wiener increments
                    dW = feval(RandFUN,1,D);
                    if ~isvector(dW) || isempty(dW) || ~isfloat(dW)
                        error('SDETools:sde_euler:RandFUNNot2DArray2',...
                             ['RandFUN must return a non-empty matrix of '...
                              'floating-point values.  See %s.'],solver);
                    end
                    [m,n] = size(dW);
                    if m ~= 1 || n ~= D
                        error('SDETools:sde_euler:RandFUNDimensionMismatch2',...
                             ['The specified alternative RandFUN did not '...
                              'output a 1 by %d column vector as '...
                              'requested.  See %s.'],D,solver);
                    end

                    dW = sh(1)*dW;
                    if ConstStep
                        RandFUN = @(i)sh*feval(RandFUN,1,D);
                    else
                        RandFUN = @(i)sh(i)*feval(RandFUN,1,D);
                    end
                end
            catch err
                switch err.identifier
                    case 'MATLAB:TooManyInputs'
                        error('SDETools:sde_euler:RandFUNTooFewInputs',...
                             ['RandFUN must have at least two inputs.  '...
                              'See %s.'],solver);
                    case 'MATLAB:TooManyOutputs'
                        error('SDETools:sde_euler:RandFUNNoOutput',...
                             ['The output of RandFUN was not specified. '...
                              'RandFUN must return a non-empty matrix.  '...
                              'See %s.'],solver);
                    case 'MATLAB:unassignedOutputs'
                        error('SDETools:sde_euler:RandFUNUnassignedOutput',...
                             ['The first output of RandFUN was not '...
                              'assigned.  See %s.'],solver);
                    case 'MATLAB:minrhs'
                        error('SDETools:sde_euler:RandFUNTooManyInputs',...
                             ['RandFUN must not require more than two '...
                              'inputs.  See %s.'],solver);
                    otherwise
                        rethrow(err);
                end
            end
        end
    else
        % No error checking needed if default RANDN used
        if dWinY                            % Store Wiener increments in Y
            if ScalarNoise || ConstStep
                Y(2:end,1:D) = sh.*feval(RandFUN,lt-1,D);
            else
                Y(2:end,1:D) = bsxfun(@times,sh,feval(RandFUN,lt-1,D));
            end
            if nargout >= 2
                W = cumsum(Y(:,1:D),1);     % Integrated Wiener increments
            end
        elseif dWinW                        % Store Wiener increments in W
            if ConstStep
                W = [zeros(1,D,dataType);sh.*feval(RandFUN,lt-1,D)];
            else
                W = [zeros(1,D,dataType);...
                     bsxfun(@times,sh,feval(RandFUN,lt-1,D))];
            end
        else                                % Unable to store Wiener increments
            if ConstStep
                RandFUN = @(i)sh*feval(RandFUN,1,D);
            else
                RandFUN = @(i)sh(i)*feval(RandFUN,1,D);
            end
            dW = RandFUN(1);
        end
    end
elseif ~isDiffusion && nargout >= 2
    W = zeros(lt,0,dataType);
end

% Integrate
if ConstFFUN && ConstGFUN && ((~isDiffusion && isYOutput) || dWinY || dWinW) ...
        && ~NonNegative
    % No FOR loop needed
    if dWinY
     	Y = cumsum(Y,1);                    % Integrate Wiener increments in Y
        if isW
            W = Y(:,1:D);                   % If needed for events or output
        end
        if isDrift
            if OneDNoise                    % 1-D scalar (and diagonal) noise
                Y = y0+tspan*fout+Y*gout;
            elseif ScalarNoise
                Y = bsxfun(@plus,y0.',bsxfun(@plus,tspan*fout.',Y(:,1)*gout.'));
            elseif DiagonalNoise
                Y = bsxfun(@plus,y0.',tspan*fout.'+bsxfun(@times,Y,gout.'));
            else
                Y = bsxfun(@plus,y0.',tspan*fout.'+Y(:,1:D)*gout.');
            end
        else
            if OneDNoise                    % 1-D scalar (and diagonal) noise
                Y = y0+Y*gout;
            elseif DiagonalNoise
                Y = bsxfun(@plus,y0.',bsxfun(@times,Y,gout.'));
            else
                Y = bsxfun(@plus,y0.',Y(:,1:D)*gout.');
            end
        end
    elseif dWinW
        W = cumsum(W,1);                    % Integrate Wiener increments in W
        if isDrift
            Y = bsxfun(@plus,y0.',tspan*fout.'+W*gout.');
        else
            Y = bsxfun(@plus,y0.',W*gout.');
        end
    else
        if isW                              % If needed for events or output
            if strcmp(dataType,'double')
                W(lt,D) = 0;
            else
                W(lt,D) = single(0);
            end
        end
        if isDrift
            if N == 1
                Y = y0+tspan*fout;
            else
                Y = bsxfun(@plus,y0.',tspan*fout.');
            end
        else
            Y = ones(lt,1,dataType)*y0.';
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
else
    dt = h(1);                              % Fixed step size
    Ti = tspan(1);                          % Current time
  	Yi = y0;                                % Set initial conditions
    
    if OneDNoise
        % Optimized 1-D scalar (and diagonal) noise case
        if isYOutput
            Y(1) = Yi;                      % Store initial conditions
        end
        
        % Integration loop using Wiener increments stored in Y(i+1)
        for i = 1:lt-1
            if ~ConstStep
                dt = h(i);                  % Step size
            end
            if dWinY
                dW = Y(i+1);                % Wiener increment
            elseif i > 1
                dW = RandFUN(i);            % Generate Wiener increment
            end
            
            % Calculate next time step
            if ConstGFUN
              	Yi = Yi+f(Ti,Yi)*dt+gout*dW;
            else
                if ~ConstFFUN
                    fout = f(Ti,Yi)*dt;
                end
                gout = g(Ti,Yi);
                
                if Stratonovich             % Use Euler-Heun step
                    Yi = Yi+fout+0.5*(gout+g(Ti,Yi+gout*dW))*dW;
                else
                    Yi = Yi+fout+gout*dW;
                end
            end
            
            % Force solution to be >= 0
            if NonNegative
                Yi = abs(Yi);
            end
            
            Ti = tspan(i+1);                % Increment current time
            if isYOutput
                Y(i+1) = Yi;              	% Store solution
            end
            
            % Integrated Wiener increments for events and output functions
            if isW
                if nargout >= 2
                    Wi = W(i+1);            % Use stored W
                else
                    Wi = Wi+dW;             % Integrate Wiener increment
                end
            end
            
            % Check for and handle zero-crossing events
            if isEvents
                [te,ye,we,ie,EventsValue,IsTerminal] ...
                    = sdezero(EventsFUN,Ti,Yi,Wi,EventsValue);
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
                        Y = Y(1:i+1);
                        if nargout >= 2
                            W = W(1:i+1);
                        end
                        break;
                    end
                end
            end
            
            % Check for and handle output function
            if isOutput
                OutputFUN(Ti,Yi,'',Wi);
            end
        end
    else
        % General case
        if isYOutput
            Y(1,:) = Yi;                	% Store initial conditions
        end
        
        % Integration loop using cached state, Yi, and Wiener increments, dW
        for i = 1:lt-1
            if ~ConstStep
                dt = h(i);                  % Variable step size
            end
            if dWinY
                dW = Y(i+1,1:D);            % Wiener increments stored in Y
            elseif dWinW
                dW = W(i+1,:);              % Wiener increments stored in W
                W(i+1,:) = W(i,:)+dW;       % Integrate Wiener increments
            elseif i > 1
                dW = RandFUN(i);            % Generate Wiener increments
            end
            dW = dW(:);                     % Wiener increments
            
            % Calculate next time step
            if ConstGFUN
                if DiagonalNoise
                    Yi = Yi+f(Ti,Yi)*dt+gout.*dW;
                else
                    Yi = Yi+f(Ti,Yi)*dt+gout*dW;
                end
            else
                if ~ConstFFUN
                    fout = f(Ti,Yi)*dt;
                end
                gout = g(Ti,Yi);

                if Stratonovich	% Use Euler-Heun step
                    if DiagonalNoise
                        Yi = Yi+fout+0.5*(gout+g(Ti,Yi+gout.*dW)).*dW;
                    else
                        Yi = Yi+fout+0.5*(gout+g(Ti,Yi+gout*dW))*dW;
                    end
                else
                    if DiagonalNoise
                        Yi = Yi+fout+gout.*dW;
                    else
                        Yi = Yi+fout+gout*dW;
                    end
                end
            end
            
            % Force specified solution indices to be >= 0
            if NonNegative
                Yi(idxNonNegative) = abs(Yi(idxNonNegative));
            end
            
            Ti = tspan(i+1);                % Increment current time
            if isYOutput
                Y(i+1,:) = Yi;           	% Store solution
            end
            
            % Integrated Wiener increments for events and output functions
            if isW
                if nargout >= 2
                    Wi = W(i+1,:);          % Use stored W
                    Wi = Wi(:);
                else
                    Wi = Wi+dW;             % Integrate Wiener increments
                end
            end
            
            % Check for and handle zero-crossing events
            if isEvents
                [te,ye,we,ie,EventsValue,IsTerminal] ...
                    = sdezero(EventsFUN,Ti,Yi(:),Wi,EventsValue);
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
                        Y = Y(1:i+1,:);
                        if nargout >= 2
                            W = W(1:i+1,:);
                        end
                        break;
                    end
                end
            end
            
            % Check for and handle output function
            if isOutput
                OutputFUN(Ti,Yi(:),'',Wi);
            end
        end
    end
end

% Finalize output
if isOutput
    OutputFUN([],[],'done',[]);
end