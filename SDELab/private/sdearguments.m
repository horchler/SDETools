function [N D tspan tdir lt y0 f0 g0 h ConstStep dataType idxNonNegative ...
          NonNegative DiagonalNoise ScalarNoise ConstFFUN ConstGFUN ...
          Stratonovich RandFUN CustomRandFUN] ...
          = sdearguments(solver,f,g,tspan,y0,options,args)
%SDEARGUMENTS  Helper function that processes arguments for all SDE solvers.
%
%   See also SDE_EULER, SDE_MILSTEIN, SDEGET, FUNCTION_HANDLE, RANDSTREAM.
        
%   Andrew D. Horchler, adh9@case.edu, Created 12-12-11
%   Revision: 1.0, 3-28-12

%   SDEARGUMENTS is partially based on an updating of version 1.12.4.15 of
%   Matlab's ODEARGUMENTS.


% Check that tspan is internally consistent
lt = length(tspan);         % number of time steps
if lt < 2
    error(  'SDELab:sdearguments:InvalidTSpanSize',...
            'Input vector TSPAN must have length >= 2.  See %s.',solver);
end
if ~isfloat(tspan) || ~isreal(tspan)
    error(  'SDELab:sdearguments:InvalidTSpanDataType',...
           ['Datatype of the input vector TSPAN must be real single or '...
            'real double.  See %s.'],solver);
end
if any(~isfinite(tspan))
    warning(    'SDELab:sdearguments:TSpanNotFinite',...
                'One or more elements of the input TSPAN are not finite.');
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error(	'SDELab:sdearguments:TspanNotMonotonic',...
           ['The entries in TSPAN must strictly increase or decrease.'...
            '  See %s.'],solver);
end
dtspan = abs(dtspan);       % length time steps
htspan = abs(tspan(2)-t0);  % length of first time step
if all(dtspan == htspan)
    h = htspan;
    ConstStep = true;
else
    h = dtspan;
    ConstStep = false;
end

if isempty(y0) || ~isfloat(y0)
    error(  'SDELab:sdearguments:Y0EmptyOrNotFloat',...
           ['The initial conditions, Y0, must be non-empty vector of '...
            'singles or doubles.  See %s.'],solver);
end
y0 = y0(:);
N = length(y0);	% number of state variables

% Ensure first solver input is function handle, or matix for constant function
if ~isa(f,'function_handle')
    if isempty(f) && ndims(f) == 2 && all(size(f) == 0) && isnumeric(f)
        f0 = 0;
    elseif ~isempty(f) && ndims(f) == 2 && isfloat(f)
        f0 = f;
    else
        error(  'SDELab:sdearguments:InvalidFFUN',...
               ['The input FFUN must be a function handle, a matrix of '...
                'singles or doubles, or an empty numeric matrix, [], '...
                'denoting no deterministic function.  See %s.'],solver);
    end
    ConstFFUN = true;
else
    % Check output of FFUN and save it
    try
        %f = @(t,y)feval(f,t,y,args{:});
        %f0 = f(t0,y0);
        f0 = feval(f,t0,y0,args{:});
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error(  'SDELab:sdearguments:FFUNTooFewInputs',...
                        'FFUN must have at least two inputs.  See %s.',solver);
            case 'MATLAB:TooManyOutputs'
                error(  'SDELab:sdearguments:FFUNNoOutput',...
                       ['The output of FFUN was not specified. FFUN must '...
                        'return a non-empty column vector.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error(  'SDELab:sdearguments:FFUNUnassignedOutput',...
                       ['The first output of FFUN was not assigned.'...
                        '  See %s.'],solver);
            case 'MATLAB:minrhs'
                error(  'SDELab:sdearguments:FFUNTooManyInputs',...
                       ['FFUN requires one or more input arguments '...
                        '(parameters) that were not supplied.'...
                        '  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
    % If deterministic function specified as constant
    ConstFFUN = strcmp(sdeget(options,'ConstFFUN','no','flag'),'yes');
end

% Ensure that size and type of output of FFUN is consistent
[m n] = size(f0);
if ndims(f0) ~= 2 || n > 1 || isempty(f0) || ~isfloat(f0)
    error(  'SDELab:sdearguments:FFUNNotColumnVector',...
           ['FFUN must return a non-empty column vector of floating point '...
            'values.  See %s.'],solver);
end
if m ~= N
    error(  'SDELab:sdearguments:FFUNDimensionMismatch',...
            'FFUN must return a vector the same length as Y0.  See %s.',solver);
end

% Check for non-negative components
idxNonNegative = sdeget(options,'NonNegative','no','flag');
if strcmp(idxNonNegative,'yes')
    idxNonNegative = 1:N;
    NonNegative = true;
elseif ~strcmp(idxNonNegative,'no') && ~isempty(idxNonNegative)
    if ~isnumeric(idxNonNegative) || ~isreal(idxNonNegative) || ...
            ~all(isfinite(idxNonNegative)) || ~isvector(idxNonNegative)
        error(  'SDELab:sdearguments:InvalidNonNegative',...
               ['NonNegative option must be a finite real numeric vector.'...
                '  See %s.'],solver);
    end
    if any(idxNonNegative < 1) || any(idxNonNegative > N) || ...
            ~all(idxNonNegative-floor(idxNonNegative) == 0)
        error(  'SDELab:sdearguments:InvalidNonNegative',...
               ['NonNegative option must be a vector of integer indices no '...
                'greater than the length of Y0.  See %s.'],solver);
    end
    if any(diff(sort(idxNonNegative)) == 0)
        error(  'SDELab:sdearguments:InvalidNonNegative',...
               ['NonNegative vector cannot contain repeated indices.'...
                '  See %s.'],solver);
    end
    y0(idxNonNegative) = max(y0(idxNonNegative),0);
    NonNegative = true;
else
    idxNonNegative = [];
    NonNegative = false;
end

% Ensure second solver input is function handle, or matix for constant function
if ~isa(g,'function_handle')
    if isempty(g) && ndims(g) == 2 && all(size(g) == 0) && isnumeric(g)
        g0 = 0;
    elseif ~isempty(g) && ndims(g) == 2 && isfloat(g)
        g0 = g;
    else
        error(  'SDELab:sdearguments:InvalidGFUN',...
               ['The input GFUN must be a function handle, a matrix of '...
                'singles or doubles, or an empty numeric matrix, [], '...
                'denoting no stochastic function.  See %s.'],solver);
    end
    ConstGFUN = true;
else
    % Check output of GFUN and save it
    try
        %g = @(t,y)feval(g,t,y,args{:});
        %g0 = g(t0,y0);
        g0 = feval(g,t0,y0,args{:});
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error(  'SDELab:sdearguments:GFUNTooFewInputs',...
                        'GFUN must have at least two inputs.  See %s.',solver);
            case 'MATLAB:TooManyOutputs'
                error(  'SDELab:sdearguments:GFUNNoOutput',...
                       ['The output of GFUN was not specified. GFUN must '...
                        'return a non-empty matrix.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error(  'SDELab:sdearguments:GFUNUnassignedOutput',...
                       ['The first output of GFUN was not assigned.'...
                        '  See %s.'],solver);
            case 'MATLAB:minrhs'
                error(  'SDELab:sdearguments:GFUNTooManyInputs',...
                       ['GFUN requires one or more input arguments '...
                        '(parameters) that were not supplied.'...
                        '  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
    % If stochastic function specified as constant
    ConstGFUN = strcmp(sdeget(options,'ConstGFUN','no','flag'),'yes');
end

% If noise is specified as diagonal, i.e., uncorrelated
DiagonalNoise = strcmp(sdeget(options,'DiagonalNoise','yes','flag'),'yes');

% Ensure that size and type of output of GFUN is consistent
if ndims(g0) ~= 2 || isempty(g0) || ~isfloat(g0)
    error(  'SDELab:sdearguments:GFUNNot2DArray',...
           ['GFUN must return a non-empty matrix of floating point values.'...
            '  See %s.'],solver);
end
[m n] = size(g0);
if DiagonalNoise
    if n ~= 1 || ~(m == N || m == 1)
        error(  'SDELab:sdearguments:GFUNNotColumnVector',...
               ['For diagonal noise, GFUN must return a scalar value or a '...
                'non-empty column vector the same length as Y0.'...
                '  See %s.'],solver);
    end
    ScalarNoise = false;
    D = N;
elseif m == 1 && n == 1 % scalar noise case same for both, doesn't depend on N
    ScalarNoise = true;
    D = 1;
else
    if m ~= N
        error(  'SDELab:sdearguments:GFUNDimensionMismatchNonDiagonal',...
               ['For non-diagonal noise, GFUN must return a scalar value or '...
                'a non-empty matrix with the same number of rows as the '...
                'length as Y0.  See %s.'],solver);
    end
    ScalarNoise = false;
    D = n;
end

% Determine the dominant data type, single or double
dataType = superiorfloat(t0,y0,f0,g0);
if ~all(strcmp(dataType,{class(t0),class(y0),class(f0),class(g0)}))
    warning( 'SDELab:sdearguments:InconsistentDataType',...
            ['Mixture of single and double data for inputs TSPAN and Y0 and '...
             'outputs of FFUN and GFUN.']);
end

% Create function handle to be used for generating Wiener increments
RandFUN = sdeget(options,'RandFUN',[],'flag');
if ~isempty(RandFUN)	% Use alternative random number generator
    if ~isa(RandFUN,'function_handle')
        error(  'SDELab:sdearguments:RandFUNNotAFunctionHandle',...
                'RandFUN must be a function handle.  See %s.',solver);
    end
    CustomRandFUN = true;
else    % Use Matlab's random number generator for normal variates
    RandSeed = sdeget(options,'RandSeed',[],'flag');
    if ~isempty(RandSeed)
        if ~isscalar(RandSeed) || ~isnumeric(RandSeed) || ...
                ~isreal(RandSeed) || ~isfinite(RandSeed) || ...
                RandSeed >= 2^32 || RandSeed < 0
            error(	'SDELab:sdearguments:InvalidRandSeed',...
                   ['RandSeed must be a non-negative integer value less '...
                    'than 2^32.  See %s.'],solver);
        end
        % Create new stream based on seed value
        Stream = RandStream.create('mt19937ar','Seed',RandSeed);
    else
        % Use default stream
        Stream = RandStream.getGlobalStream;
    end
    
    % Set property if antithetic random variates option is specified
    set(Stream,'Antithetic',strcmp(sdeget(options,'Antithetic','no','flag'),'yes'));
    
    RandFUN = @(M,N)randn(Stream,M,N,dataType);
    CustomRandFUN = false;
end

% Integration method is dependent on if SDE is Stratonovich or Ito form
if ConstGFUN || strcmp(sdeget(options,'AdditiveNoise','no','flag'),'yes')
    % Stochastic function is constant or additive, i.e., not a function of state
    Stratonovich = false;
else
    Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich','flag'),'Stratonovich');
end