function [OutputFUN,WSelect]=sdeoutputfun(solver,tspan,y0,N,options)
%SDEOUTPUTFUN  Process output function arguments for SDE solvers and functions.
%
%   See also:
%       SDEARGUMENTS, SDEARGUMENTS_PROCESS, SDEEVENTSFUN, SDERANDFUN, SDEGET,
%       FUNCTION_HANDLE
        
%   Andrew D. Horchler, horchler @ gmail . com, Created 5-2-13
%   Revision: 1.2, 5-6-13


% Check for output function
OutputFUN = sdeget(options,'OutputFUN',[],'flag');
if ~isempty(OutputFUN)
    if ~isa(OutputFUN,'function_handle')
        error('SDETools:sdeoutputfun:OutputFUNNotAFunctionHandle',...
              'OutputFUN, if specified, must be a function handle.  See %s.',...
              solver);
    end
    
    % Check for selected Y components
    idxY = sdeget(options,'OutputYSelect','yes','flag');
    if strcmp(idxY,'yes')
        YSelect = true;
        YAll = true;
    elseif ~strcmp(idxY,'no') && ~isempty(idxY)
        if ~isnumeric(idxY) || ~isreal(idxY) || ~all(isfinite(idxY)) ...
                || ~isvector(idxY)
            error('SDETools:sdeoutputfun:InvalidOutputYSelect',...
                 ['OutputYSelect option must be a finite real numeric '...
                  'vector.  See %s.'],solver);
        end
        if any(idxY < 1) || any(idxY > N) || ~all(idxY == floor(idxY))
            error('SDETools:sdeoutputfun:InvalidIndexOutputYSelect',...
                 ['OutputYSelect option must be a vector of integer indices '...
                  'no greater than the length of Y0.  See %s.'],solver);
        end
        if any(diff(sort(idxY)) == 0)
            error('SDETools:sdeoutputfun:RepeatedIndexOutputYSelect',...
                 ['OutputYSelect vector cannot contain repeated indices.  '...
                  'See %s.'],solver);
        end
        idxY = idxY(:);
        YSelect = true;
        YAll = false;
    else
        YSelect = false;
        YAll = false;
    end
    
    % Check for selected W components, create function handles
    idxW = sdeget(options,'OutputWSelect','no','flag');
    if strcmp(idxW,'yes')
        if YAll
            OutputFUN = @(t,y,flag,w)OutputFUN(t,y,flag,w);
        elseif YSelect
            OutputFUN = @(t,y,flag,w)OutputFUN(t,y(idxY(:,~isempty(y))),flag,w);
        else
            OutputFUN = @(t,y,flag,w)OutputFUN(t,[],flag,w);
        end
        WSelect = true;
    elseif ~strcmp(idxW,'no') && ~isempty(idxW)
        if ~isnumeric(idxW) || ~isreal(idxW) || ~all(isfinite(idxW)) ...
                || ~isvector(idxW)
            error('SDETools:sdeoutputfun:InvalidOutputWSelect',...
                 ['OutputWSelect option must be a finite real numeric '...
                  'vector.  See %s.'],solver);
        end
        if any(idxW < 1) || any(idxW > N) || ~all(idxW == floor(idxW))
            error('SDETools:sdeoutputfun:InvalidIndexOutputWSelect',...
                 ['OutputWSelect option must be a vector of integer indices '...
                  'no greater than the length of Y0.  See %s.'],solver);
        end
        if any(diff(sort(idxW)) == 0)
            error('SDETools:sdeoutputfun:RepeatedIndexOutputWSelect',...
                 ['OutputWSelect vector cannot contain repeated indices.  '...
                  'See %s.'],solver);
        end
        idxW = idxW(:);
        if YAll
            OutputFUN = @(t,y,flag,w)OutputFUN(t,y,flag,w(idxW(:,~isempty(w))));
        elseif YSelect
            OutputFUN = @(t,y,flag,w)OutputFUN(t,y(idxY(:,~isempty(y))),flag,...
                w(idxW(:,~isempty(w))));
        else
            OutputFUN = @(t,y,flag,w)OutputFUN(t,[],flag,...
                w(idxW(:,~isempty(w))));
        end
        WSelect = true;
    else
        if YSelect && ~YAll
            OutputFUN = @(t,y,flag,w)OutputFUN(t,y(idxY(:,~isempty(y))),flag);
        elseif ~YSelect
            OutputFUN = @(t,y,flag,w)OutputFUN(t,[],flag);
        else
            OutputFUN = @(t,y,flag,w)OutputFUN(t,y,flag);
        end
        WSelect = false;
    end
    
    % Initialize and check output function
    try
        if WSelect
            OutputFUN(tspan,y0,'init',zeros(N,1));
        else
            OutputFUN(tspan,y0,'init',[]);
        end
    catch err
        switch err.identifier
            case 'MATLAB:TooManyInputs'
                if WSelect
                    error('SDETools:sdeoutputfun:OutputFUNTooFewInputs4',...
                         ['OutputFUN must have at least four inputs.  See '...
                          'SDEPLOT.']);
                else
                    error('SDETools:sdeoutputfun:OutputFUNTooFewInputs3',...
                         ['OutputFUN must have at least three inputs.  See '...
                          'SDEPLOT.']);
                end
            case 'MATLAB:unassignedOutputs'
                error('SDETools:sdeoutputfun:OutputFUNUnassignedOutput',...
                     ['The first output of OutputFUN was not assigned.  See '...
                      'SDEPLOT.']);
            case 'MATLAB:minrhs'
                error('SDETools:sdeoutputfun:OutputFUNTooManyInputs',...
                     ['OutputFUN requires one or more input arguments '...
                      '(parameters) that were not supplied.  See SDEPLOT.']);
            otherwise
                rethrow(err);
        end
    end
else
    WSelect = false;
end