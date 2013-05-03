function [RandFUN,CustomRandFUN,ResetStream]=sderandfun(solver,dataType,options)
%SDERANDFUN  Process random function arguments for SDE solvers and functions.
%
%   See also:
%       SDEARGUMENTS, SDEARGUMENTS_PROCESS, SDEEVENTSFUN, SDEOUTPUTFUN, SDEGET,
%       FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 5-2-13
%   Revision: 1.2, 5-2-13


% Create function handle to be used for generating Wiener increments
RandFUN = sdeget(options,'RandFUN',[],'flag');
if ~isempty(RandFUN)	% Use alternative random number generator
    if ~isa(RandFUN,'function_handle')
        error('SDETools:sderandfun:RandFUNNotAFunctionHandle',...
              'RandFUN must be a function handle.  See %s.',solver);
    end
    CustomRandFUN = true;
    ResetStream = [];
else    % Use Matlab's random number generator for normal variates
    Stream = sdeget(options,'RandStream',[],'flag');
    if ~isempty(Stream)
        if ~isa(Stream,'RandStream')
            error('SDETools:sderandfun:InvalidRandStream',...
                  'RandStream must be a RandStream object.  See %s.',solver);
        end
    else
        RandSeed = sdeget(options,'RandSeed',[],'flag');
        if ~isempty(RandSeed)
            if ~isscalar(RandSeed) || ~isnumeric(RandSeed) ...
                    || ~isreal(RandSeed) || ~isfinite(RandSeed) ...
                    || RandSeed >= 2^32 || RandSeed < 0
                error('SDETools:sderandfun:InvalidRandSeed',...
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