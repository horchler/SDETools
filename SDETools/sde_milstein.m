function [Y,W,TE,YE,WE,IE] = sde_milstein(f,g,tspan,y0,options)
%SDE_MILSTEIN  Solve stochastic differential equations, Milstein method.
%   YOUT = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0) with TSPAN = [T0 T1 ... TFINAL]
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
%   [YOUT, W] = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0) outputs the M-by-N matrix W of
%   integrated Weiner increments that were used by the solver. Each row of W
%   corresponds to a time in TSPAN.
%
%   [...] = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0,OPTIONS) solves as above with
%   default integration properties replaced by values in OPTIONS, an argument
%   created with the SDESET function. See SDESET for details. The type of SDE to
%   be integrated, 'Ito' or the default 'Stratonovich', can be specified via the
%   SDEType property. By default, the derivative-free Milstein method is
%   applied, but if the if the DGFUN property is set to a function handle or
%   floating point matrix to specify the derivative of the noise function
%   g(t,y), then the general Milstein method is used. The DiagonalNoise property
%   must be set to 'no' in order to apply the more general correlated noise case
%   where GFUN may return (or is) an N-by-D matrix, where D is the dimension of
%   the noise. Another commonly used option is to manually specify the random
%   seed via the RandSeed property, which creates a new random number stream,
%   instead of using the default stream, to generate the Wiener increments.
%
%   [YOUT, W, TE, YE, WE, IE] = SDE_MILSTEIN(FFUN,GFUN,TSPAN,Y0,OPTIONS) with
%   the EventsFUN property set to a function handle, in order to specify an
%   events function, solves as above while also finding zero-crossings. The
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
%       % Solve 2-D Stratonovich SDE using Milstein method with derivative
%       mu = 1; sig = [0.1;0.5]; dt = 1e-2; t = 0:dt:1;
%       f = @(t,y)mu.*y; g = @(t,y)sig.*y; opts = sdeset('DGFUN',sig);
%       y = sde_milstein(f,g,t,[1 1]); plot(t,y);
%       title(['Milstein Method, dt = ' num2str(dt) ', \mu = ' num2str(mu)]);
%       txt = {['\sigma = ' num2str(sig(1))],['\sigma = ' num2str(sig(2))]};
%       legend(txt,2); legend boxoff; xlabel('t'); ylabel('y(t)');
%
%   Notes:
%       SDEs are assumed to be in Stratonovich form by default. SDEs in Ito form
%       can be handled by setting the SDEType OPTIONS property to 'Ito' via
%       the SDESET function. These forms are generally not equivalent and will
%       converge to different solutions, so care should be taken to ensure that
%       the form of SDEFUN matches SDEType.
%
%       In the case of additve noise, i.e., when diffusion term, g(t,y), is not
%       a function of the state variables, both Ito and Stratonovich
%       interpretations are equivalent and the Milstein method reduces to the
%       Euler-Maruyama method implemented in SDE_EULER. If GFUN(T,Y) is
%       specified as a floating point matrix rather than a function handle then
%       it is automatically assumed to be constant. Alternatively, the ConstGFUN
%       property can be set to 'yes' to indicate that GFUN(T,Y) is constant and
%       only needs to be evaluated once.
%
%   See also:
%       Other SDE solvers:      SDE_EULER
%       Implicit SDE solvers:      
%       Stochastic processes:	SDE_BM, SDE_GBM, SDE_OU
%       Option handling:        SDESET, SDEGET
%       SDE demos/validation:	SDE_EULER_VALIDATE, SDE_MILSTEIN_VALIDATE
%       Other:                  FUNCTION_HANDLE, RANDSTREAM

%   SDE_MILSTEIN is an implementation of the order 1.0 strong (order 1.0 weak)
%   explicit Milstein scheme, which is also an order 1.0 strong Taylor
%   approximation scheme.

%   For details of this integration method, see: Peter E. Kloeden and Eckhard
%   Platen, "Numerical solution of Stochastic Differential Equations,"
%   Springer-Verlag, 1992.

%   Andrew D. Horchler, adh9 @ case . edu, 10-25-10
%   Revision: 1.0, 1-12-13


solver = 'SDE_MILSTEIN';

% Check inputs and outputs
if nargin < 5
    if nargin < 4
        error('SDETools:sde_milstein:NotEnoughInputs',...
              'Not enough input arguments.  See %s.',solver);
    end
    if isa(y0,'struct')
        error('SDETools:sde_milstein:NotEnoughInputsOptions',...
             ['An SDE options structure was provided as the last argument, '...
              'but one of the first four input arguments is missing.'...
              '  See %s.'],solver);
    end
    options = [];
elseif isempty(options) && (~sde_ismatrix(options) ...
        || any(size(options) ~= 0) || ~(isstruct(options) || iscell(options) ...
        || isnumeric(options))) || ~isempty(options) && ~isstruct(options)
	error('SDETools:sde_milstein:InvalidSDESETStruct',...
          'Invalid SDE options structure.  See SDESET.');
end

% Handle solver arguments (NOTE: ResetStream is called by onCleanup())
[N,D,D0,tspan,tdir,lt,y0,fout,gout,h,ConstStep,dataType,idxNonNegative,...
    NonNegative,DiagonalNoise,ScalarNoise,idxConstFFUN,ConstFFUN,...
    idxConstGFUN,ConstGFUN,Stratonovich,RandFUN,CustomRandFUN,ResetStream,...
    EventsFUN,EventsValue] = sdearguments(solver,f,g,tspan,y0,options);

% If optional derivative function property is set
if ~ConstGFUN
    dg = sdeget(options,'DGFUN',[],'fast');

    % Ensure DGFUN is function handle, or matrix for constant function
    if ~isa(dg,'function_handle')
        if isempty(dg) && sde_ismatrix(dg) && all(size(dg) == 0) ...
                && isnumeric(dg)
            Derivative = false;
        elseif ~isempty(dg) && sde_ismatrix(dg) && isfloat(dg)
            dgout = dg;
            ConstDGFUN = true;
            Derivative = true;
            if all(dgout == 0)
                ConstGFUN = true;
            end
        else
            error('SDETools:sde_milstein:InvalidDGFUN',...
                 ['The input DGFUN must be a function handle, a matrix of '...
                  'floating point values, or an empty numeric matrix, [], '...
                  'denoting no stochastic derivative function.  See %s.'],...
                  solver);
        end
    else
        % Check output of DGFUN and save it
        try
            dgout = dg(t0,y0);
        catch err
            switch err.identifier
                case 'MATLAB:TooManyInputs'
                    error('SDETools:sde_milstein:DGFUNTooFewInputs',...
                          'DGFUN must have at least two inputs.  See %s.',...
                          solver);
                case 'MATLAB:TooManyOutputs'
                    error('SDETools:sde_milstein:DGFUNNoOutput',...
                         ['The output of DGFUN was not specified. DGFUN '...
                          'must return a non-empty matrix.  See %s.'],solver);
                case 'MATLAB:unassignedOutputs'
                    error('SDETools:sde_milstein:DGFUNUnassignedOutput',...
                         ['The first output of DGFUN was not assigned.'...
                          '  See %s.'],solver);
                case 'MATLAB:minrhs'
                    error('SDETools:sde_milstein:DGFUNTooManyInputs',...
                         ['DGFUN requires one or more input arguments '...
                          '(parameters) that were not supplied.  See %s.'],...
                          solver);
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
        if ~sde_ismatrix(dgout) || isempty(dgout) || ~isfloat(dgout)
            error('SDETools:sde_milstein:DGFUNNot2DArray',...
                 ['DGFUN must return a non-empty matrix of floating point '...
                  'values.  See %s.'],solver);
        end
        [m,n] = size(dgout);
        if DiagonalNoise
            if n ~= 1 || ~(m == N || m == 1)
                error('SDETools:sde_milstein:DGFUNNotColumnVector',...
                     ['For diagonal noise, DGFUN must return a scalar value '...
                      'or a non-empty column vector the same length as Y0.'...
                      '  See %s.'],solver);
            end
        else
            if m ~= N
                error('SDETools:sde_milstein:DGFUNDimensionMismatchNonDiag',...
                     ['For non-diagonal noise, DGFUN must return a scalar '...
                      'value or a non-empty matrix with the same number of '...
                      'rows as the length as Y0.  See %s.'],solver);
            end
        end
    end
end

% Initialize outputs for zero-crossing events
isEvents = ~isempty(EventsFUN);
if isEvents
    if nargout > 6
        error('SDETools:sde_milstein:EventsTooManyOutputs',...
              'Too many output arguments.  See %s.',solver);
    else
        if nargout >= 3
            TE = [];
            if nargout >= 4
                YE = [];
                if nargout >= 5
                    WE = [];
                    if nargout >= 6
                        IE = [];
                    end
                end
            end
        end
    end
else
    if nargout > 2
        if nargout <= 6
            error('SDETools:sde_milstein:NoEventsTooManyOutputs',...
                 ['Too many output arguments. An events function has not '...
                  'been specified.  See %s.'],solver);
        else
            error('SDETools:sde_milstein:TooManyOutputs',...
                  'Too many output arguments.  See %s.',solver);
        end
    end
end

% State array
if strcmp(dataType,'double')
    Y(lt,N) = 0;
else
    Y(lt,N) = single(0);
end

% Calculate Wiener increments from normal variates, store in state if possible
sh = tdir*sqrt(h);
h = tdir*h;
if CustomRandFUN    % check output of alternative RandFUN
    try
        if D > N
            if nargout >= 2 % Store Wiener increments in W
                W = feval(RandFUN,lt-1,D);
                if ~sde_ismatrix(W) || isempty(W) || ~isfloat(W)
                    error('SDETools:sde_milstein:RandFUNNot2DArray1',...
                         ['RandFUN must return a non-empty matrix of '...
                          'floating point values.  See %s.'],solver);
                end
                [m,n] = size(W);
                if m ~= lt-1 || n ~= D
                    error('SDETools:sde_milstein:RandFUNDimensionMismatch1',...
                         ['The specified alternative RandFUN did not output '...
                          'a %d by %d matrix as requested.  See %s.'],lt-1,D,...
                          solver);
                end
                if ConstStep
                    W = [zeros(1,D);sh*W];
                else
                    W = [zeros(1,D);bsxfun(@times,sh,W)];
                end
            else            % Unable to store Wiener increments
                dW = feval(RandFUN,1,D);
                if ~isvector(dW) || isempty(dW) || ~isfloat(dW)
                    error('SDETools:sde_milstein:RandFUNNot2DArray2',...
                         ['RandFUN must return a non-empty matrix of '...
                          'floating point values.  See %s.'],solver);
                end
                [m,n] = size(dW);
                if m ~= 1 || n ~= D
                    error('SDETools:sde_milstein:RandFUNDimensionMismatch2',...
                         ['The specified alternative RandFUN did not output '...
                          'a 1 by %d column vector as requested.  See %s.'],...
                          D,solver);
                end
                dW = sh(1)*dW;
            end
        elseif D == N       % Store Wiener increments in Y indirectly
            r = feval(RandFUN,lt-1,D);
            if ~sde_ismatrix(r) || isempty(r) || ~isfloat(r)
                error('SDETools:sde_milstein:RandFUNNot2DArray3',...
                     ['RandFUN must return a non-empty matrix of floating '...
                      'point values.  See %s.'],solver);
            end
            [m,n] = size(r);
            if m ~= lt-1 || n ~= D
                error('SDETools:sde_milstein:RandFUNDimensionMismatch3',...
                     ['The specified alternative RandFUN did not output a '...
                      '%d by %d matrix as requested.   See %s.',lt-1,D,solver]);
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
                error('SDETools:sde_milstein:RandFUNTooFewInputs',...
                      'RandFUN must have at least two inputs.  See %s.',solver);
            case 'MATLAB:TooManyOutputs'
                error('SDETools:sde_milstein:RandFUNNoOutput',...
                     ['The output of RandFUN was not specified. RandFUN '...
                      'must return a non-empty matrix.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error('SDETools:sde_milstein:RandFUNUnassignedOutput',...
                     ['The first output of RandFUN was not assigned.'...
                      '  See %s.'],solver);
            case 'MATLAB:minrhs'
                error('SDETools:sde_milstein:RandFUNTooManyInputs',...
                     ['RandFUN must not require more than two inputs.'...
                      '  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
else    % No error checking needed if default RANDN used
    if D > N
        if nargout >= 2 % Store Wiener increments in W
            if ConstStep
                W(2:end,:) = sh*feval(RandFUN,lt-1,D);
            else
                W(2:end,:) = bsxfun(@times,sh,feval(RandFUN,lt-1,D));
            end
        else            % Unable to store Wiener increments
            dW = sh(1)*feval(RandFUN,1,D);
        end
    else                % Store Wiener increments in Y
        if D > 0 && ~isempty(D0)
            if N == 1 && DiagonalNoise || ConstStep
                Y(2:end,D0) = sh.*feval(RandFUN,lt-1,D);
            else
                Y(2:end,D0) = bsxfun(@times,sh,feval(RandFUN,lt-1,D));
            end
        end
    end
end

% Only allocate W matrix if requested as output
if nargout >= 2 && D <= N
	W = cumsum(Y(:,D0),1);
end

% Integration loop
if ConstFFUN && ConstGFUN && (D <= N || nargout >= 2) && ~NonNegative	% no FOR loop needed
    if D > N
        W = cumsum(W,1);
    elseif nargout < 2
        Y = cumsum(Y(:,D0),1);
        if isEvents
            W = Y;
        end
    end
    
    % Evaluate analytic solution
    if N == 1 && DiagonalNoise
        Y = y0+tspan*fout+Y*gout;
    else
        if ScalarNoise
            Y = bsxfun(@plus,y0.',bsxfun(@plus,tspan*fout.',Y*gout));
        elseif DiagonalNoise
            Y = bsxfun(@plus,y0.',tspan*fout.'+bsxfun(@times,Y(:,D0),gout(D0).'));
        else
            if D > N
                Y = bsxfun(@plus,y0.',tspan*fout.'+W*gout.');
            else
                Y = bsxfun(@plus,y0.',tspan*fout.'+Y*gout.');
            end
        end
    end
    
    % Check for and handle zero-crossing events
    if isEvents
        for i = 2:lt
            [te,ye,we,ie,EventsValue,IsTerminal] = sdezero(EventsFUN,tspan(i),Y(i,:),W(i,:),EventsValue);
            if ~isempty(te)
                if nargout >= 3
                    TE = [TE;te];               %#ok<AGROW>
                    if nargout >= 4
                        YE = [YE;ye];           %#ok<AGROW>
                        if nargout >= 5
                            WE = [WE;we];       %#ok<AGROW>
                            if nargout >= 6
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
                    return;
                end
            end
        end
    end
else
    % Optimized integration loop for 1-D scalar noise case
    if N == 1 && DiagonalNoise
        dt = h(1);  % Step size for first time step and fixed step-size
        Y(1) = y0;	% Set initial conditions
        dW = Y(2);  % Wiener increment
        
        % Use existing fout and gout to calculate the first time step, Y(2)
        if ConstGFUN
            Y(2) = y0+fout*dt+gout*dW;
        else
            if Derivative
                if Stratonovich
                    dW2 = 0.5*dW^2;
                else
                    dW2 = 0.5*(dW^2-dt);
                end
                Y(2) = y0+fout*dt+gout*(dW+dgout*dW2);
            else
                sdt = sh(1);
                if Stratonovich
                    dW2 = (0.5/sdt)*dW^2;
                else
                    dW2 = (0.5/sdt)*(dW^2-dt);
                end
              	Ybar = y0+gout*dW;
                Y(2) = Ybar+fout*dt+(g(tspan(1),Ybar)-gout)*dW2;
            end
        end
        if NonNegative
            Y(2) = max(Y(2),0);
        end
        
        % Check for and handle zero-crossing events
        if isEvents
            Wi = dW;
            [te,ye,we,ie,EventsValue,IsTerminal] = sdezero(EventsFUN,tspan(2),Y(2),Wi,EventsValue);
            if ~isempty(te)
                if nargout >= 3
                    TE = [TE;te];
                    if nargout >= 4
                        YE = [YE;ye];
                        if nargout >= 5
                            WE = [WE;we];
                            if nargout >= 6
                                IE = [IE;ie];
                            end
                        end
                    end
                end
                if IsTerminal
                    Y = Y(1:2);
                    if nargout >= 2
                        W = W(1:2);
                    end
                    return;
                end
            end
        end
        
        % Integration loop using Weiner increments stored in Y(i+1)
        for i = 2:lt-1
            % Step size and Wiener increments
            if ~ConstStep
                dt = h(i);
                sdt = sh(i);
            end
            dW = Y(i+1);    % Wiener increment

            % Calculate next time step
            if ConstGFUN
              	Y(i+1) = Y(i)+f(tspan(i),Y(i))*dt+gout*dW;
            else
                if ~ConstFFUN
                    fout = f(tspan(i),Y(i))*dt;
                end

                gout = g(tspan(i),Y(i));
                if Derivative   % Use Milstein with derivative if specified
                    if Stratonovich
                        dW2 = 0.5*dW^2;
                    else
                        dW2 = 0.5*(dW^2-dt);
                    end
                    if ~ConstDGFUN
                        dgout = dg(tspan(i),Y(i));
                    end
                    Y(i+1) = Y(i)+fout+gout*(dW+dgout*dW2);
                else            % Otherwise use derivative-free Milstein method
                    if Stratonovich
                        dW2 = (0.5/sdt)*dW^2;
                    else
                        dW2 = (0.5/sdt)*(dW^2-dt);
                    end
                 	Ybar = Y(i)+gout*dW;
                	Y(i+1) = Ybar+fout+(g(tspan(i),Ybar)-gout)*dW2;
                end
            end
            if NonNegative
                Y(i+1) = max(Y(i+1),0);
            end
            
            % Check for and handle zero-crossing events
            if isEvents
                if nargout >= 2
                    Wi = W(i+1);
                else
                    Wi = Wi+dW; % Integrate Wiener increment
                end
                [te,ye,we,ie,EventsValue,IsTerminal] = sdezero(EventsFUN,tspan(i+1),Y(i+1),Wi,EventsValue);
                if ~isempty(te)
                    if nargout >= 3
                        TE = [TE;te];               %#ok<AGROW>
                        if nargout >= 4
                            YE = [YE;ye];           %#ok<AGROW>
                            if nargout >= 5
                                WE = [WE;we];       %#ok<AGROW>
                                if nargout >= 6
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
                        return;
                    end
                end
            end
        end
    else
        % step size and Wiener increment for first time step
        dt = h(1);
        sdt = sh(1);
        if D > N 
            if nargout >= 2         % dW already exists if D > N && nargin == 1
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
                    Yi = Ybar+fout*dt+(g(tspan(1),Ybar)-gout).*dW2;
                else
                    Ybar = y0+gout*dW;
                    Yi = Ybar+fout*dt+(g(tspan(1),Ybar)-gout)*dW2;
                end
            end
        end
        if NonNegative
            Yi(idxNonNegative) = max(Yi(idxNonNegative),0);
        end
        Y(2,:) = Yi;
        
        % Check for and handle zero-crossing events
        if isEvents
            Wi = dW;
            [te,ye,we,ie,EventsValue,IsTerminal] = sdezero(EventsFUN,tspan(2),Yi,Wi,EventsValue);
            if ~isempty(te)
                if nargout >= 3
                    TE = [TE;te];
                    if nargout >= 4
                        YE = [YE;ye];
                        if nargout >= 5
                            WE = [WE;we];
                            if nargout >= 6
                                IE = [IE;ie];
                            end
                        end
                    end
                end
                if IsTerminal
                    Y = Y(1:2,:);
                    if nargout >= 2
                        W = W(1:2,:);
                    end
                    return;
                end
            end
        end

        % Integration loop using cached state, Yi, and increment, dW
        for i = 2:lt-1
            % Step size and Wiener increments
            if ~ConstStep
                dt = h(i);
                sdt = sh(i);
            end
            if D > N
                if nargout >= 2
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
                    Yi = Yi+f(tspan(i),Yi)*dt+gout.*dW;
                else
                    Yi = Yi+f(tspan(i),Yi)*dt+gout*dW;
                end
            else
                if ~ConstFFUN
                    fout = f(tspan(i),Yi)*dt;
                end

                gout = g(tspan(i),Yi);
                if Derivative   % Use Milstein with derivative function if specified
                    if Stratonovich
                        dW2 = 0.5*dW.^2;
                    else
                        dW2 = 0.5*(dW.^2-dt);
                    end
                    if ~ConstDGFUN
                        dgout = dg(tspan(i),Yi);
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
                        Yi = Ybar+fout+(g(tspan(i),Ybar)-gout).*dW2;
                    else
                        Ybar = Yi+gout*dW;
                        Yi = Ybar+fout+(g(tspan(i),Ybar)-gout)*dW2;
                    end
                end
            end
            if NonNegative
                Yi(idxNonNegative) = max(Yi(idxNonNegative),0);
            end
            Y(i+1,:) = Yi;
            
            % Check for and handle zero-crossing events
            if isEvents
                if nargout >= 2
                    Wi = W(i+1,:);
                else
                    Wi = Wi+dW; % Integrate Wiener increments
                end
                [te,ye,we,ie,EventsValue,IsTerminal] = sdezero(EventsFUN,tspan(i+1),Yi,Wi,EventsValue);
                if ~isempty(te)
                    if nargout >= 3
                        TE = [TE;te];               %#ok<AGROW>
                        if nargout >= 4
                            YE = [YE;ye];           %#ok<AGROW>
                            if nargout >= 5
                                WE = [WE;we];       %#ok<AGROW>
                                if nargout >= 6
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
                        return;
                    end
                end
            end
        end
    end
end