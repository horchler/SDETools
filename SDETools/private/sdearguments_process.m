function [N,tspan,tdir,lt,y0,h,ConstStep,Stratonovich,RandFUN,ResetStream,...
          EventsFUN,EventsValue,OutputFUN,WSelect] ...
          = sdearguments_process(func,tspan,y0,dataType,options)
%SDEARGUMENTS_PROCESS  Process arguments for all SDE process functions.
%
%   See also:
%       SDE_BM, SDE_GBM, SDE_OU, SDEARGUMENTS, SDEGET, SDEEVENTS, SDEZERO,
%       SDEOUTPUT, SDERESET_STREAM, FUNCTION_HANDLE
        
%   Andrew D. Horchler, horchler @ gmail . com, Created 4-4-12
%   Revision: 1.2, 5-3-13

%   sdearguments_process is partially based on an updating of version 1.12.4.15
%   of Matlab's ODEARGUMENTS.


% Check that tspan is internally consistent
lt = length(tspan);         % Number of time steps
if lt < 2
    error('SDETools:sdearguments_process:InvalidTSpanSize',...
          'Input vector TSPAN must have length >= 2.  See %s.',func);
end
if ~isfloat(tspan) || ~isreal(tspan)
    error('SDETools:sdearguments_process:InvalidTSpanDataType',...
         ['Datatype of the input vector TSPAN must be real single or real '...
          'double.  See %s.'],func);
end
if any(~isfinite(tspan))
    warning('SDETools:sdearguments_process:TSpanNotFinite',...
           ['One or more elements of the input TSPAN are not finite'...
            '  See %s.'],func);
end
tspan = tspan(:);
t0 = tspan(1);
tdir = sign(tspan(end)-t0);
dtspan = diff(tspan);
if tdir == 0 || (tdir > 0 && any(dtspan <= 0)) || (tdir < 0 && any(dtspan >= 0))
	error('SDETools:sdearguments_process:TspanNotMonotonic',...
         ['The entries in TSPAN must strictly increase or decrease.'...
          '  See %s.'],func);
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

% Check y0
if isempty(y0) || ~isfloat(y0)
    error('SDETools:sdearguments_process:Y0EmptyOrNotFloat',...
         ['The initial conditions, Y0, must be non-empty vector of singles '...
          'or doubles.  See %s.'],func);
end
y0 = y0(:);
N = length(y0);             % Number of state variables

% Check for events function
[EventsFUN,EventsValue] = sdeeventsfun(func,t0,y0,options);

% Check for output function
[OutputFUN,WSelect] = sdeoutputfun(func,tspan,y0,N,options);

% Create function handle to be used for generating Wiener increments
[RandFUN,ResetStream] = sderandfun(func,dataType,options);

% Check if noise is not specified as diagonal, i.e., uncorrelated
if strcmp(sdeget(options,'DiagonalNoise','yes','flag'),'no');
    error('SDETools:sdearguments_process:NonDiagonalNoiseUnsupported',...
         ['The DiagonalNoise property is set to ''no'', but this function '...
          'only supports diagonal noise, not the general correlated noise '...
          'case.  See %s.'],func);
end

% Check if non-negative property is specified
if ~strcmp(sdeget(options,'NonNegative','no','flag'),'no')
    error('SDETools:sdearguments_process:NonNegativeUnsupported',...
         ['The NonNegative property is set to a value other than ''no'', '...
          'but this function does not support this option.  See %s.'],func);
end

% Solution method is dependent on if SDE is Stratonovich or Ito form
Stratonovich = strcmp(sdeget(options,'SDEType','Stratonovich','flag'),...
    'Stratonovich');