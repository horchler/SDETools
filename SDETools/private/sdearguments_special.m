function [N,tspan,tdir,lt,y0,h,ConstStep,Stratonovich,RandFUN,CustomRandFUN,...
          ResetStream,EventsFUN,EventsValue]...
          = sdearguments_special(func,tspan,y0,dataType,options,args)
%SDEARGUMENTS_SPECIAL  Process arguments for all SDE special functions.
%
%   See also:
%       SDE_GBM, SDE_OU, SDEARGUMENTS, SDEGET, FUNCTION_HANDLE, RANDSTREAM
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 4-4-12
%   Revision: 1.0, 12-31-12

%   SDEARGUMENTS_SPECIAL is partially based on an updating of version 1.12.4.15
%   of Matlab's ODEARGUMENTS.


% Check that tspan is internally consistent
lt = length(tspan);         % number of time steps
if lt < 2
    error('SDETools:sdearguments_special:InvalidTSpanSize',...
          'Input vector TSPAN must have length >= 2.  See %s.',func);
end
if ~isfloat(tspan) || ~isreal(tspan)
    error('SDETools:sdearguments_special:InvalidTSpanDataType',...
         ['Datatype of the input vector TSPAN must be real single or real '...
          'double.  See %s.'],func);
end
if any(~isfinite(tspan))
    warning('SDETools:sdearguments_special:TSpanNotFinite',...
           ['One or more elements of the input TSPAN are not finite'...
            '  See %s.'],func);
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error('SDETools:sdearguments_special:TspanNotMonotonic',...
         ['The entries in TSPAN must strictly increase or decrease.'...
          '  See %s.'],func);
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

% Check y0
if isempty(y0) || ~isfloat(y0)
    error('SDETools:sdearguments_special:Y0EmptyOrNotFloat',...
         ['The initial conditions, Y0, must be non-empty vector of singles '...
          'or doubles.  See %s.'],func);
end
y0 = y0(:);
N = length(y0);	% number of state variables

% Create function handle to be used for generating Wiener increments
RandFUN = sdeget(options,'RandFUN',[],'flag');
if ~isempty(RandFUN)	% Use alternative random number generator
    if ~isa(RandFUN,'function_handle')
        error('SDETools:sdearguments_special:RandFUNNotAFunctionHandle',...
              'RandFUN must be a function handle.  See %s.',func);
    end
    CustomRandFUN = true;
    ResetStream = [];
else    % Use Matlab's random number generator for normal variates
    Stream = sdeget(options,'RandStream',[],'flag');
    if ~isempty(Stream)
        if ~isa(Stream,'RandStream')
            error('SDETools:sdearguments_special:InvalidRandStream',...
                  'RandStream must be a RandStream object.  See %s.',solver);
        end
    else
        RandSeed = sdeget(options,'RandSeed',[],'flag');
        if ~isempty(RandSeed)
            if ~isscalar(RandSeed) || ~isnumeric(RandSeed) ...
                    || ~isreal(RandSeed) || ~isfinite(RandSeed) ...
                    || RandSeed >= 2^32 || RandSeed < 0
                error('SDETools:sdearguments_special:InvalidRandSeed',...
                     ['RandSeed must be a non-negative integer value less '...
                      'than 2^32.  See %s.'],func);
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
    
    % Function to be called on completion or early termination of integration
    ResetStream = onCleanup(@()reset_stream(Stream));
end

% Solution method is dependent on if SDE is Stratonovich or Ito form
Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich','flag'),...
    'Stratonovich');

% Check for events function
EventsFUN = sdeget(options,'Events',[],'flag');
if ~isempty(EventsFUN)
    if ~isa(EventsFUN,'function_handle')
        error('SDETools:sdearguments_special:EventsFUNNotAFunctionHandle',...
              'EventsFUN, if specified, must be a function handle.  See %s.',...
              solver);
    end
    
    % Check output of EventsFUN at initial condition and save value
    try
        [EventsValue,isterminal,direction] = feval(EventsFUN,tspan(1),y0,args{:});
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                error('SDETools:sdearguments_special:EventsFUNTooFewInputs',...
                      'EventsFUN must have at least two inputs.  See %s.',...
                      solver);
            case 'MATLAB:TooManyOutputs'
                error('SDETools:sdearguments_special:EventsFUNNoOutput',...
                     ['The output of EventsFUN was not specified. EventsFUN '...
                      'must return three non-empty vectors.  See %s.'],solver);
            case 'MATLAB:unassignedOutputs'
                error('SDETools:sdearguments_special:EventsFUNUnassignedOutput',...
                     ['The first output of EventsFUN was not assigned.  '...
                      'See %s.'],solver);
            case 'MATLAB:minrhs'
                error('SDETools:sdearguments_special:EventsFUNTooManyInputs',...
                     ['EventsFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.  See %s.'],solver);
            otherwise
                rethrow(err);
        end
    end
    if ~isvector(EventsValue) || isempty(EventsValue)...
            || ~all(isfinite(EventsValue)) 
        error('SDETools:sdearguments_special:InvalidEventsValue',...
             ['The first output of EventsFUN, ''Value'', must be a '...
              'non-empty finite vector.  See %s.'],solver);
    end
    if ~isvector(isterminal)...
            || ~any(length(isterminal) == [length(EventsValue) 1])
        error('SDETools:sdearguments_special:EventsIsterminalDimensionMismatch',...
             ['The second output of EventsFUN, ''IsTerminal'', must be a '...
              'scalar or a vector the same length as the first output.  '...
              'See %s.'],solver);
    end
    if ~all(isterminal == 0 | isterminal == 1)
        error('SDETools:sdearguments_special:InvalidEventsIsterminal',...
             ['The elements of the second output of EventsFUN, '...
              '''IsTerminal'', must be equal to 0 (false) or 1 (true).  '...
              'See %s.'],solver);
    end
    if ~isempty(direction)
        if ~isvector(direction)...
                || ~any(length(direction) == [length(EventsValue) 1])
            error('SDETools:sdearguments_special:EventsDirectionDimensionMismatch',...
                 ['If the third output of EventsFUN, ''Direction'', is not '...
                  'specified as the empty matrix, [], it must be a scalar '...
                  'or a vector the same length as first output.  See %s.'],...
                  solver);
        end
        if ~all(direction == 0 | direction == 1 | direction == -1)
            error('SDETools:sdearguments_special:InvalidEventsDirection',....
                 ['The third output of EventsFUN, ''Direction'', must be '...
                  'equal to 0 (default), 1, or -1.  See %s.'],solver);
        end
    end
else
    EventsValue = [];
end