function sdereset_stream(Stream)
%SDERESET_STREAM  Reset global stream or delete custom random stream.
%
%   See also:
%       SDEARGUMENTS, SDEARGUMENTS_PROCESS, ONCLEANUP, FUNCTION_HANDLE,
%       RANDSTREAM
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-30-12
%   Revision: 1.2, 5-2-13


% Called by onCleanup to reset antihetic property for global stream
try
    isGlobal = isequal(Stream,RandStream.getGlobalStream);
catch                                                       %#ok<CTCH>
    isGlobal = isequal(Stream,RandStream.getDefaultStream);	%#ok<GETRS>
end
if isGlobal
    set(Stream,'Antithetic',~Stream.Antithetic);
else
    delete(Stream);
end