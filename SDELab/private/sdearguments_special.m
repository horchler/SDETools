function [N tspan tdir lt y0 h ConstStep Stratonovich RandFUN CustomRandFUN] ...
          = sdearguments_special(func,tspan,y0,options,dataType)
%SDEARGUMENTS_SPECIAL  Processes arguments for all SDE special functions.
%
%   See also:
%       SDE_GBM, SDE_OU, SDEARGUMENTS, SDEGET, FUNCTION_HANDLE, RANDSTREAM
        
%   Andrew D. Horchler, adh9@case.edu, Created 4-4-12
%   Revision: 1.0, 4-8-12

%   SDEARGUMENTS_SPECIAL is partially based on an updating of version 1.12.4.15
%   of Matlab's ODEARGUMENTS.


% Check that tspan is internally consistent
lt = length(tspan);         % number of time steps
if lt < 2
    error(  'SDELab:sdearguments_special:InvalidTSpanSize',...
            'Input vector TSPAN must have length >= 2.  See %s.',func);
end
if ~isfloat(tspan) || ~isreal(tspan)
    error(  'SDELab:sdearguments_special:InvalidTSpanDataType',...
           ['Datatype of the input vector TSPAN must be real single or '...
            'real double.  See %s.'],func);
end
if any(~isfinite(tspan))
    warning(    'SDELab:sdearguments_special:TSpanNotFinite',...
               ['One or more elements of the input TSPAN are not finite'...
                '  See %s.'],func);
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error(	'SDELab:sdearguments_special:TspanNotMonotonic',...
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
    error(  'SDELab:sdearguments_special:Y0EmptyOrNotFloat',...
           ['The initial conditions, Y0, must be non-empty vector of '...
            'singles or doubles.  See %s.'],func);
end
y0 = y0(:);
N = length(y0);	% number of state variables

% Create function handle to be used for generating Wiener increments
RandFUN = sdeget(options,'RandFUN',[],'flag');
if ~isempty(RandFUN)	% Use alternative random number generator
    if ~isa(RandFUN,'function_handle')
        error(  'SDELab:sdearguments_special:RandFUNNotAFunctionHandle',...
                'RandFUN must be a function handle.  See %s.',func);
    end
    CustomRandFUN = true;
else    % Use Matlab's random number generator for normal variates
    RandSeed = sdeget(options,'RandSeed',[],'flag');
    if ~isempty(RandSeed)
        if ~isscalar(RandSeed) || ~isnumeric(RandSeed) || ...
                ~isreal(RandSeed) || ~isfinite(RandSeed) || ...
                RandSeed >= 2^32 || RandSeed < 0
            error(	'SDELab:sdearguments_special:InvalidRandSeed',...
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
    set(Stream,'Antithetic',strcmp(sdeget(options,'Antithetic','no','flag'),...
        'yes'));
    
    RandFUN = @(M,N)randn(Stream,M,N,dataType);
    CustomRandFUN = false;
end

% Solution method is dependent on if SDE is Stratonovich or Ito form
Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich','flag'),...
    'Stratonovich');