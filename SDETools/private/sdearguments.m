function [N,D,tspan,tdir,lt,y0,f0,g0,dg0,dg,h,ConstStep,dataType,NonNegative,...
          idxNonNegative,DiagonalNoise,ScalarNoise,OneDNoise,ConstFFUN,...
          ConstGFUN,ConstDGFUN,Stratonovich,RandFUN,ResetStream,EventsFUN,...
          EventsValue,OutputFUN,WSelect] ...
          = sdearguments(solver,f,g,tspan,y0,options)
%SDEARGUMENTS  Process arguments for all SDE solvers.
%
%   See also:
%       SDE_EULER, SDE_MILSTEIN, SDEARGUMENTS_PROCESS, SDEGET, SDEEVENTS,
%       SDEZERO, SDEOUTPUT, SDERESET_STREAM, FUNCTION_HANDLE
        
%   Andrew D. Horchler, horchler @ gmail . com, Created 12-12-11
%   Revision: 1.3, 6-8-16

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
          'TSPAN vector must be a real floating-point value.  See %s.',solver);
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
%{
% Check if max step size specified
MaxStep = sdeget(options,'MaxStep',0,'flag');
if MaxStep(:) ~= 0
    if ~isscalar(MaxStep) || ~isfloat(MaxStep) || ~isreal(MaxStep)
        
    end
    if ~isfinite(MaxStep) || MaxStep < 0
        
    end
    
end
%}
% Check for non-negative components
idxNonNegative = sdeget(options,'NonNegative','no','flag');
if strcmp(idxNonNegative,'yes')
    idxNonNegative = 1:N;
    y0 = abs(y0);
    NonNegative = true;
elseif ~strcmp(idxNonNegative,'no') && ~isempty(idxNonNegative)
    if ~isnumeric(idxNonNegative) || ~isreal(idxNonNegative) ...
            || ~all(isfinite(idxNonNegative)) || ~isvector(idxNonNegative)
        error('SDETools:sdearguments:InvalidNonNegative',...
             ['NonNegative option must be a finite real numeric vector.'...
              '  See %s.'],solver);
    end
    if any(idxNonNegative < 1) || any(idxNonNegative > N) ...
            || ~all(idxNonNegative == floor(idxNonNegative))
        error('SDETools:sdearguments:InvalidIndexNonNegative',...
             ['NonNegative option must be a vector of integer indices no '...
              'greater than the length of Y0.  See %s.'],solver);
    end
    if any(diff(sort(idxNonNegative)) == 0)
        error('SDETools:sdearguments:RepeatedIndexNonNegative',...
              'NonNegative vector cannot contain repeated indices.  See %s.',...
              solver);
    end
    y0(idxNonNegative) = abs(y0(idxNonNegative));
    NonNegative = true;
else
    idxNonNegative = [];
    NonNegative = false;
end

% Check for events function
[EventsFUN,EventsValue] = sdeeventsfun(solver,t0,y0,options);

% Check for output function
[OutputFUN,WSelect] = sdeoutputfun(solver,tspan,y0,N,options);

% Ensure first solver input is function handle, or vector for constant function
if ~isa(f,'function_handle')
    if isempty(f) && isvector(f) && all(size(f) == 0) && isnumeric(f)
        f0 = zeros(1,class(f));   	% Driftless diffusion
    elseif ~isempty(f) && isvector(f) && isfloat(f)
        f0 = f(:);
        if length(f0) ~= N
            error('SDETools:sdearguments:FFUNDimensionMismatchVector',...
                 ['FFUN must return a vector the same length as Y0.  ',...
                  'See %s.'],solver);
        end
        if all(f0 == 0)
            f0 = zeros(1,class(f));	% Driftless diffusion
        end
    else
        error('SDETools:sdearguments:InvalidFFUN',...
             ['The input FFUN must be a function handle, a vector of '...
              'floating-point values, or an empty numeric matrix, [], '...
              'denoting no deterministic function.  See %s.'],solver);
    end
    ConstFFUN = true;   % FFUN is constant
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
    
    % Ensure that size and type of output of FFUN is consistent
    [m,n] = size(f0);
    if ~isvector(f0) || n > 1 || isempty(f0) || ~isfloat(f0)
        error('SDETools:sdearguments:FFUNNotColumnVector',...
             ['FFUN must return a non-empty column vector of floating-point '...
              'values.  See %s.'],solver);
    end
    if m ~= N
        error('SDETools:sdearguments:FFUNDimensionMismatch',...
              'FFUN must return a vector the same length as Y0.  See %s.',...
              solver);
    end
    
    % Check if FFUN specified as constant
    ConstFFUN = strcmp(sdeget(options,'ConstFFUN','no','flag'),'yes');
    if ConstFFUN && all(f0 == 0)
        f0 = zeros(1,class(f0));	% Driftless diffusion
    end
end

% Ensure second solver input is function handle, or matrix for constant function
if ~isa(g,'function_handle')
    if isempty(g) && sde_ismatrix(g) && all(size(g) == 0) && isnumeric(g)
        g0 = zeros(1,class(g));     % Noiseless drift
    elseif ~isempty(g) && sde_ismatrix(g) && isfloat(g)
        g0 = g;
        if all(g0(:) == 0)
            g0 = zeros(1,class(g));	% Noiseless drift
        end
    else
        error('SDETools:sdearguments:InvalidGFUN',...
             ['The input GFUN must be a function handle, a matrix of '...
              'floating-point values, or an empty numeric matrix, [], '...
              'denoting no stochastic function.  See %s.'],solver);
    end
    ConstGFUN = true;   % GFUN is constant
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
    if ~sde_ismatrix(g0) || isempty(g0) || ~isfloat(g0)
        error('SDETools:sdearguments:GFUNNotMatrix',...
             ['GFUN must return a non-empty matrix of floating-point values.'...
              '  See %s.'],solver);
    end
    
    % Check if GFUN specified as constant
    ConstGFUN = strcmp(sdeget(options,'ConstGFUN','no','flag'),'yes');
    if ConstGFUN && all(g0(:) == 0)
        g0 = zeros(1,class(g0));	% Noiseless drift
    end
end

% If diffusion function exists (not constant and all zero)
isDiffusion = ~(ConstGFUN && isscalar(g0) && g0 == 0);
if isDiffusion
    [d,D] = size(g0);
    
    if ConstGFUN
        dg0 = zeros(0,class(g0));
        dg = [];
        ConstDGFUN = false;
        Derivative = false;
    else
        % If GFUN not constant and optional derivative function property set
        dg = sdeget(options,'DGFUN',[],'fast');

        % Ensure DGFUN is function handle, or matrix for constant function
        if ~isa(dg,'function_handle')
            if isempty(dg) && sde_ismatrix(dg) && all(size(dg) == 0) ...
                    && isnumeric(dg)
                dg0 = zeros(0,class(dg));
                ConstDGFUN = false;
                Derivative = false;
            elseif ~isempty(dg) && sde_ismatrix(dg) && isfloat(dg)
                dg0 = dg;
                ConstDGFUN = true;  % DGFUN is constant
                Derivative = true;
            else
                error('SDETools:sde_milstein:InvalidDGFUN',...
                     ['The input DGFUN must be a function handle, a matrix '...
                      'of floating-point values, or an empty numeric '...
                      'matrix, [], denoting no stochastic derivative '...
                      'function.  See %s.'],solver);
            end
        else
            % If stochastic derivative function specified as constant
            ConstDGFUN = strcmp(sdeget(options,'ConstDGFUN','no','flag'),'yes');

            % Check output of DGFUN and save it
            if ConstDGFUN || ~strcmp(solver,'SDE_EULER')
                try
                    dg0 = dg(t0,y0);
                catch err
                    switch err.identifier
                        case 'MATLAB:TooManyInputs'
                            error('SDETools:sde_milstein:DGFUNTooFewInputs',...
                                 ['DGFUN must have at least two inputs.  '...
                                  'See %s.'],solver);
                        case 'MATLAB:TooManyOutputs'
                            error('SDETools:sde_milstein:DGFUNNoOutput',...
                                 ['The output of DGFUN was not specified. '...
                                  'DGFUN must return a non-empty matrix.  '...
                                  'See %s.'],solver);
                        case 'MATLAB:unassignedOutputs'
                            error('SDETools:sde_milstein:DGFUNUnassignedOutput',...
                                 ['The first output of DGFUN was not '...
                                  'assigned.  See %s.'],solver);
                        case 'MATLAB:minrhs'
                            error('SDETools:sde_milstein:DGFUNTooManyInputs',...
                                 ['DGFUN requires one or more input '...
                                  'arguments (parameters) that were not '...
                                  'supplied.  See %s.'],solver);
                        otherwise
                            rethrow(err);
                    end
                end
                if ~sde_ismatrix(dg0) || isempty(dg0) || ~isfloat(dg0)
                    error('SDETools:sde_milstein:DGFUNNot2DArray',...
                         ['DGFUN must return a non-empty matrix of '...
                          'floating-point values.  See %s.'],solver);
                end
                Derivative = true;
            else
                dg0 = zeros(1,class(g0));
                Derivative = false;
            end
        end
        
        % Ensure that size of output of GFUN and DGFUN is consistent
        if Derivative && any([d D] ~= size(dg0))
            error('SDETools:sdearguments:GFUNDimensionMismatchDGFUN',...
                 ['The output of DGFUN must be the same dimension as the '...
                  'output of GFUN.  See %s.'],solver);
        end
        if ConstDGFUN && all(dg0(:) == 0)
            ConstGFUN = true;   % Derivative is always zero, so GFUN is constant
            if all(g0(:) == 0)
                g0 = zeros(1,class(g0));	% Noiseless drift
            end
        end
    end
else
    % No diffusion function (or all zero), use defaults
    D = 1;
    dg0 = zeros(0,class(g0));
    dg = [];
    ConstDGFUN = false;
    Derivative = false;
end

% Determine the dominant data type, single or double
if Derivative
    if ~all(strcmp(class(t0),{class(y0),class(f0),class(g0),class(dg0)}))
        warning('SDETools:sdearguments:InconsistentDataTypeDerivative',...
               ['Mixture of single and double data for inputs TSPAN and Y0 '...
                'and outputs of FFUN, GFUN, and DGFUN.']);
    end
    dataType = superiorfloat(t0,y0,f0,g0,dg0);
else
    % Determine the dominant data type, single or double
    if ~all(strcmp(class(t0),{class(y0),class(f0),class(g0)}))
        warning('SDETools:sdearguments:InconsistentDataType',...
               ['Mixture of single and double data for inputs TSPAN and Y0 '...
                'and outputs of FFUN and GFUN.']);
    end
    dataType = superiorfloat(t0,y0,f0,g0);
end

% Check again if diffusion function exists, ConstGFUN and g0 updated
if isDiffusion && ~(ConstGFUN && isscalar(g0) && g0 == 0)
    % If noise is specified as diagonal, i.e., uncorrelated
    DiagonalNoise = strcmp(sdeget(options,'DiagonalNoise','yes','flag'),'yes');
    
    if DiagonalNoise
        if D ~= 1 || (d ~= N && d ~= 1)
            error('SDETools:sdearguments:GFUNDimensionMismatchDiagonal',...
                 ['For diagonal noise, GFUN must return a non-empty scalar '...
                  'or a column vector the same length as Y0. Set the '...
                  'DiagonalNoise option to ''no'' with SDESET for general '...
                  'noise case.  See %s.'],solver);
        end
        D = N;
        ScalarNoise = (N == 1);
    else
        if d ~= N && d ~= 1
            error('SDETools:sdearguments:GFUNDimensionMismatchNonDiagonal',...
                 ['For non-diagonal noise, GFUN must return a non-empty '...
                  'scalar or matrix with the same number of rows as the '...
                  'length as Y0.  See %s.'],solver);
        end
        if D == 1
            ScalarNoise = true;
            DiagonalNoise = (N == 1);
        elseif ConstGFUN && sde_isdiag(g0)
            D = N;
            g0 = g0(1:N+1:end).';
            if Derivative
                dg0 = dg0(1:N+1:end).';
            end
            ScalarNoise = false;
            DiagonalNoise = true;   % Not specified, noise is actually diagonal
        else
            ScalarNoise = false;
            DiagonalNoise = false;
        end
    end
    
    % GFUN is 1-by-1 scalar and diagonal, i.e., N = 1, D = 1
    OneDNoise = (DiagonalNoise && ScalarNoise);
    
    % Create function handle to be used for generating Wiener increments
    [RandFUN,ResetStream] = sderandfun(solver,dataType,options);

    % Integration method is dependent on if SDE is Stratonovich or Ito form
    if ConstGFUN
        % Stochastic function constant or additive, i.e., not function of state
        Stratonovich = false;
    else
        Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich',...
            'flag'),'Stratonovich');
    end
else
    % No diffusion function (or all zero), use defaults
    DiagonalNoise = true;
    ScalarNoise = true;
    OneDNoise = (N==1);
    Stratonovich = false;
    RandFUN = [];
    ResetStream = [];
end