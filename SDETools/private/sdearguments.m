function [N,D,D0,tspan,tdir,lt,y0,f0,g0,h,ConstStep,dataType,idxNonNegative,...
          NonNegative,DiagonalNoise,ScalarNoise,idxConstFFUN,ConstFFUN,...
          idxConstGFUN,ConstGFUN,Stratonovich,RandFUN,CustomRandFUN,...
          ResetStream,EventsFUN,EventsValue]...
          = sdearguments(solver,f,g,tspan,y0,options)
%SDEARGUMENTS  Process arguments for all SDE solvers.
%
%   See also:
%       SDE_EULER, SDE_MILSTEIN, SDEARGUMENTS_SPECIAL, SDEGET, FUNCTION_HANDLE,
%       RANDSTREAM
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-12-11
%   Revision: 1.0, 1-12-13

%   SDEARGUMENTS is partially based on an updating of version 1.12.4.15 of
%   Matlab's ODEARGUMENTS.


% Check that tspan is internally consistent
lt = length(tspan);         % Number of time steps
if lt < 2
    error('SDETools:sdearguments:InvalidTSpanSize',...
          'Input vector TSPAN must have length >= 2.  See %s.',solver);
end
if ~isfloat(tspan) || ~isreal(tspan)
    error('SDETools:sdearguments:InvalidTSpanDataType',...
          'TSPAN vector must be a real floating point value.  See %s.',solver);
end
if any(~isfinite(tspan))
    warning('SDETools:sdearguments:TSpanNotFinite',...
            'One or more elements of the input TSPAN are not finite.');
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error('SDETools:sdearguments:TspanNotMonotonic',...
         ['The entries in TSPAN must strictly increase or decrease.'...
          '  See %s.'],solver);
end
dtspan = abs(dtspan);       % Length time steps
htspan = abs(tspan(2)-t0);  % Length of first time step
if all(dtspan == htspan)
    h = htspan;
    ConstStep = true;
else
    h = dtspan;
    ConstStep = false;
end

if isempty(y0) || ~isfloat(y0)
    error('SDETools:sdearguments:Y0EmptyOrNotFloat',...
         ['The initial conditions, Y0, must be non-empty vector of singles '...
          'or doubles.  See %s.'],solver);
end
y0 = y0(:);
N = length(y0);             % Number of state variables

% Check for non-negative components
idxNonNegative = sdeget(options,'NonNegative','no','flag');
if strcmp(idxNonNegative,'yes')
    idxNonNegative = 1:N;
    y0 = max(y0,0);
    NonNegative = true;
elseif ~strcmp(idxNonNegative,'no') && ~isempty(idxNonNegative)
    if ~isnumeric(idxNonNegative) || ~isreal(idxNonNegative) ...
            || ~all(isfinite(idxNonNegative)) || ~isvector(idxNonNegative)
        error('SDETools:sdearguments:InvalidNonNegative',...
             ['NonNegative option must be a finite real numeric vector.'...
              '  See %s.'],solver);
    end
    if any(idxNonNegative < 1) || any(idxNonNegative > N) ...
            || ~all(idxNonNegative-floor(idxNonNegative) == 0)
        error('SDETools:sdearguments:InvalidIndexNonNegative',...
             ['NonNegative option must be a vector of integer indices no '...
              'greater than the length of Y0.  See %s.'],solver);
    end
    if any(diff(sort(idxNonNegative)) == 0)
        error('SDETools:sdearguments:RepeatedIndexNonNegative',...
              'NonNegative vector cannot contain repeated indices.  See %s.',...
              solver);
    end
    y0(idxNonNegative) = max(y0(idxNonNegative),0);
    NonNegative = true;
else
    idxNonNegative = [];
    NonNegative = false;
end

% Ensure first solver input is function handle, or matrix for constant function
if ~isa(f,'function_handle')
    if isempty(f) && isvector(f) && all(size(f) == 0) && isnumeric(f)
        f0 = 0;
    elseif ~isempty(f) && isvector(f) && isfloat(f)
        f0 = f(:);
    else
        error('SDETools:sdearguments:InvalidFFUN',...
             ['The input FFUN must be a function handle, a vector of '...
              'floating point values, or an empty numeric matrix, [], '...
              'denoting no deterministic function.  See %s.'],solver);
    end
    idxConstFFUN = 1:N;
    ConstFFUN = true;
else
    % Check output of FFUN and save it
    try
        f0 = f(t0,y0);
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error('SDETools:sdearguments:FFUNTooFewInputs',...
                      'FFUN must have at least two inputs.  See %s.',solver);
            case 'MATLAB:TooManyOutputs'
                error('SDETools:sdearguments:FFUNNoOutput',...
                     ['The output of FFUN was not specified. FFUN must '...
                      'return a non-empty column vector.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error('SDETools:sdearguments:FFUNUnassignedOutput',...
                      'The first output of FFUN was not assigned.  See %s.',...
                      solver);
            case 'MATLAB:minrhs'
                error('SDETools:sdearguments:FFUNTooManyInputs',...
                     ['FFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
    
    % Check for constant FFUN components
    idxConstFFUN = sdeget(options,'ConstFFUN','no','flag');
    if strcmp(idxConstFFUN,'yes')
        idxConstFFUN = 1:N;
        ConstFFUN = true;
    elseif ~strcmp(idxConstFFUN,'no') && ~isempty(idxConstFFUN)
        if ~isnumeric(idxConstFFUN) || ~isreal(idxConstFFUN) ...
                || ~all(isfinite(idxConstFFUN)) || ~isvector(idxConstFFUN)
            error('SDETools:sdearguments:InvalidConstFFUN',...
                 ['ConstFFUN option must be a finite real numeric vector.'...
                  '  See %s.'],solver);
        end
        if any(idxConstFFUN < 1) || any(idxConstFFUN > N) ...
                || ~all(idxConstFFUN-floor(idxConstFFUN) == 0)
            error('SDETools:sdearguments:InvalidIndexConstFFUN',...
                 ['ConstFFUN option must be a vector of integer indices no '...
                  'greater than the length of Y0.  See %s.'],solver);
        end
        if any(diff(sort(idxConstFFUN)) == 0)
            error('SDETools:sdearguments:RepeatedIndexConstFFUN',...
                 ['ConstFFUN vector cannot contain repeated indices.'...
                  '  See %s.'],solver);
        end
        ConstFFUN = true;
    else
        idxConstFFUN = [];
        ConstFFUN = false;
    end
end

% Ensure that size and type of output of FFUN is consistent
[m,n] = size(f0);
if ~isvector(f0) || n > 1 || isempty(f0) || ~isfloat(f0)
    error('SDETools:sdearguments:FFUNNotColumnVector',...
         ['FFUN must return a non-empty column vector of floating point '...
          'values.  See %s.'],solver);
end
if m ~= N
    error('SDETools:sdearguments:FFUNDimensionMismatch',...
          'FFUN must return a vector the same length as Y0.  See %s.',solver);
end

% Ensure second solver input is function handle, or matrix for constant function
if ~isa(g,'function_handle')
    if isempty(g) && sde_ismatrix(g) && all(size(g) == 0) && isnumeric(g)
        g0 = 0;
    elseif ~isempty(g) && sde_ismatrix(g) && isfloat(g)
        g0 = g;
    else
        error('SDETools:sdearguments:InvalidGFUN',...
             ['The input GFUN must be a function handle, a matrix of '...
              'floating point values, or an empty numeric matrix, [], '...
              'denoting no stochastic function.  See %s.'],solver);
    end
    idxConstGFUN = 1:size(g0,1);
    ConstGFUN = true;
else
    % Check output of GFUN and save it
    try
        g0 = g(t0,y0);
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error('SDETools:sdearguments:GFUNTooFewInputs',...
                      'GFUN must have at least two inputs.  See %s.',solver);
            case 'MATLAB:TooManyOutputs'
                error('SDETools:sdearguments:GFUNNoOutput',...
                     ['The output of GFUN was not specified. GFUN must '...
                      'return a non-empty matrix.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error('SDETools:sdearguments:GFUNUnassignedOutput',...
                      'The first output of GFUN was not assigned.  See %s.',...
                      solver);
            case 'MATLAB:minrhs'
                error('SDETools:sdearguments:GFUNTooManyInputs',...
                     ['GFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
    
    % Check for constant GFUN components
    idxConstGFUN = sdeget(options,'ConstGFUN','no','flag');
    if strcmp(idxConstGFUN,'yes')
        idxConstGFUN = 1:size(g0,1);
        ConstGFUN = true;
    elseif ~strcmp(idxConstGFUN,'no') && ~isempty(idxConstGFUN)
        if ~isnumeric(idxConstGFUN) || ~isreal(idxConstGFUN) ...
                || ~all(isfinite(idxConstGFUN)) || ~isvector(idxConstGFUN)
            error('SDETools:sdearguments:InvalidConstGFUN',...
                 ['ConstGFUN option must be a finite real numeric vector.'...
                  '  See %s.'],solver);
        end
        if any(idxConstGFUN < 1) || any(idxConstGFUN > N) ...
                || ~all(idxConstGFUN-floor(idxConstGFUN) == 0)
            error('SDETools:sdearguments:InvalidIndexConstGFUN',...
                 ['ConstGFUN option must be a vector of integer indices no '...
                  'greater than the length of Y0.  See %s.'],solver);
        end
        if any(diff(sort(idxConstGFUN)) == 0)
            error('SDETools:sdearguments:RepeatedIndexConstGFUN',...
                 ['ConstGFUN vector cannot contain repeated indices.'...
                  '  See %s.'],solver);
        end
        ConstGFUN = true;
    else
        idxConstGFUN = [];
        ConstGFUN = false;
    end
end

% If noise is specified as diagonal, i.e., uncorrelated
DiagonalNoise = strcmp(sdeget(options,'DiagonalNoise','yes','flag'),'yes');

% Ensure that size and type of output of GFUN is consistent
if ~sde_ismatrix(g0) || isempty(g0) || ~isfloat(g0)
    error('SDETools:sdearguments:GFUNNotMatrix',...
         ['GFUN must return a non-empty matrix of floating point values.'...
          '  See %s.'],solver);
end
[m,n] = size(g0);
if DiagonalNoise
    if n ~= 1 || ~(m == N || m == 1)
        error('SDETools:sdearguments:GFUNDimensionMismatchDiagonal',...
             ['For diagonal noise, GFUN must return a scalar value or a '...
              'non-empty column vector the same length as Y0.  See %s.'],...
              solver);
    end
	ScalarNoise = (N == 1);
    D = N;
elseif m == 1 && n == 1	% Scalar noise doesn't depend on N
    ScalarNoise = true;
    D = 1;
else
    if m ~= N
        error('SDETools:sdearguments:GFUNDimensionMismatchNonDiagonal',...
             ['For non-diagonal noise, GFUN must return a scalar value or a '...
              'non-empty matrix with the same number of rows as the length '...
              'as Y0.  See %s.'],solver);
    end
    ScalarNoise = false;
    D = n;
end
if ConstGFUN && ~DiagonalNoise
    D0 = find(g0(idxConstGFUN) ~= 0);
    D = length(D0);
else
    D0 = 1:D;
end

% Determine the dominant data type, single or double
if ~all(strcmp(class(t0),{class(y0),class(f0),class(g0)}))
    warning('SDETools:sdearguments:InconsistentDataType',...
           ['Mixture of single and double data for inputs TSPAN and Y0 and '...
            'outputs of FFUN and GFUN.']);
end
dataType = superiorfloat(t0,y0,f0,g0);

% Create function handle to be used for generating Wiener increments
RandFUN = sdeget(options,'RandFUN',[],'flag');
if ~isempty(RandFUN)	% Use alternative random number generator
    if ~isa(RandFUN,'function_handle')
        error('SDETools:sdearguments:RandFUNNotAFunctionHandle',...
              'RandFUN must be a function handle.  See %s.',solver);
    end
    CustomRandFUN = true;
    ResetStream = [];
else                % Use Matlab's random number generator for normal variates
    Stream = sdeget(options,'RandStream',[],'flag');
    if ~isempty(Stream)
        if ~isa(Stream,'RandStream')
            error('SDETools:sdearguments:InvalidRandStream',...
                  'RandStream must be a RandStream object.  See %s.',solver);
        end
    else
        RandSeed = sdeget(options,'RandSeed',[],'flag');
        if ~isempty(RandSeed)
            if ~isscalar(RandSeed) || ~isnumeric(RandSeed) ...
                    || ~isreal(RandSeed) || ~isfinite(RandSeed) ...
                    || RandSeed >= 2^32 || RandSeed < 0
                error('SDETools:sdearguments:InvalidRandSeed',...
                     ['RandSeed must be a non-negative integer value less '...
                      'than 2^32.  See %s.'],solver);
            end
            % Create new stream based on seed value
            Stream = RandStream.create('mt19937ar','Seed',RandSeed);
        else
            % Use default stream
            try
                Stream = RandStream.getGlobalStream;
            catch                                       %#ok<CTCH>
                Stream = RandStream.getDefaultStream;	%#ok<GETRS>
            end
        end
        
        % Set property if antithetic random variates option is specified
        Antithetic = strcmp(sdeget(options,'Antithetic','no','flag'),'yes');
        if Antithetic ~= Stream.Antithetic
            set(Stream,'Antithetic',Antithetic);
        end
    end
    
    RandFUN = @(M,N)randn(Stream,M,N,dataType);
	CustomRandFUN = false;
    
    % Function to be call on completion or early termination of integration
    ResetStream = onCleanup(@()sdereset_stream(Stream));
end

% Integration method is dependent on if SDE is Stratonovich or Ito form
if ConstGFUN
    % Stochastic function is constant or additive, i.e., not a function of state
    Stratonovich = false;
else
    Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich','flag'),...
        'Stratonovich');
end

% Check for events function
EventsFUN = sdeget(options,'EventsFUN',[],'flag');
if ~isempty(EventsFUN)
    if ~isa(EventsFUN,'function_handle')
        error('SDETools:sdearguments:EventsFUNNotAFunctionHandle',...
              'EventsFUN, if specified, must be a function handle.  See %s.',...
              solver);
    end
    
    % Check output of EventsFUN at initial condition and save value
    try
        [EventsValue,isterminal,direction] = EventsFUN(tspan(1),y0);
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error('SDETools:sdearguments:EventsFUNTooFewInputs',...
                      'EventsFUN must have at least two inputs.  See %s.',...
                      solver);
            case 'MATLAB:TooManyOutputs'
                error('SDETools:sdearguments:EventsFUNNoOutput',...
                     ['The output of EventsFUN was not specified. EventsFUN '...
                      'must return three non-empty vectors.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error('SDETools:sdearguments:EventsFUNUnassignedOutput',...
                     ['The first output of EventsFUN was not assigned.  '...
                      'See %s.'],solver);
            case 'MATLAB:minrhs'
                error('SDETools:sdearguments:EventsFUNTooManyInputs',...
                     ['EventsFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
    if ~isvector(EventsValue) || isempty(EventsValue)...
            || ~all(isfinite(EventsValue)) 
        error('SDETools:sdearguments:InvalidEventsValue',...
             ['The first output of EventsFUN, ''Value'', must be a '...
              'non-empty finite vector.  See %s.'],solver);
    end
    if ~isvector(isterminal)...
            || ~any(length(isterminal) == [length(EventsValue) 1])
        error('SDETools:sdearguments:EventsIsterminalDimensionMismatch',...
             ['The second output of EventsFUN, ''IsTerminal'', must be a '...
              'scalar or a vector the same length as the first output.  '...
              'See %s.'],solver);
    end
    if ~all(isterminal == 0 | isterminal == 1)
        error('SDETools:sdearguments:InvalidEventsIsterminal',...
             ['The elements of the second output of EventsFUN, '...
              '''IsTerminal'', must be equal to 0 (false) or 1 (true).  '...
              'See %s.'],solver);
    end
    if ~isempty(direction)
        if ~isvector(direction)...
                || ~any(length(direction) == [length(EventsValue) 1])
            error('SDETools:sdearguments:EventsDirectionDimensionMismatch',...
                 ['If the third output of EventsFUN, ''Direction'', is not '...
                  'specified as the empty matrix, [], it must be a scalar '...
                  'or a vector the same length as first output.  See %s.'],...
                  solver);
        end
        if ~all(direction == 0 | direction == 1 | direction == -1)
            error('SDETools:sdearguments:InvalidEventsDirection',....
                 ['The third output of EventsFUN, ''Direction'', must be '...
                  'equal to 0 (default), 1, or -1.  See %s.'],solver);
        end
    end
else
    EventsValue = [];
end