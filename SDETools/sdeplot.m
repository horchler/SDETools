function status=sdeplot(t,y,flag,w)
%SDEPLOT  Time series SDE output function.
%   When the function SDEPLOT is passed to an SDE solver as the 'OutputFUN'
%   property, i.e., OPTIONS = SDESET('OutputFUN',@SDEPLOT), the solver calls 
%   SDEPLOT(T,Y,'') after every timestep. The SDEPLOT function plots all
%   components of the solution, Y, as they are computed, adapting the axis
%   limits of the plot dynamically. To plot only particular components of the
%   solution, specify their indices via the 'OutputYSelect' SDE options
%   property. SDEPLOT is the default output function of the solvers when they
%   are called with no output arguments.
%
%   At the start of integration, the solver calls SDEPLOT(TSPAN,Y0,'init') to
%   initialize the output function. After each integration step to new time, T,
%   and solution vector, Y, the solver calls STATUS = SDEPLOT(T,Y,''). If the
%   solver's 'Refine' property is greater than one (see SDESET), T is a column
%   vector of new output times and Y is an array of corresponding column
%   vectors. The STATUS return value is 1 if the figure window and plot axis are
%   still open and 0 otherwise. When the integration is complete, the solver
%   calls SDEPLOT([],[],'done').
%
%   Set the 'OutputWSelect' SDE options property to 'yes' or to a vector of
%   indices to output the integrated Wiener increments, W. The integrated Wiener
%   increments are passed to SDEPLOT as a fourth argument, SDEPLOT(T,Y,'',W). If
%   the 'OutputWSelect' property is set to 'yes' or a non-empty vector of
%   indices the solver calls SDEPLOT(TSPAN,Y0,'init',W0) to initialize the
%   output function at the start of integration and SDEPLOT([],[],'done',[])
%   when the integration is complete.
%   
%   See also:
%       SDEPHASEPLOT2, SDEPHASEPLOT3, SDEIMGPLOT, SDEPRINT, SDESET, SDEGET,
%       SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU, ODEPLOT

%   SDEPLOT is based on an updating of Matlab's ODEPLOT, version 1.25.4.9

%   Andrew D. Horchler, adh9 @ case . edu, 4-29-13
%   Revision: 1.2, 7-16-13


persistent FIG_HANDLE AX_HANDLE LEN_TSPAN;
status = 1;         % Figure widow still open and and plot axis still present
isW = (nargin > 3); % Have integrated Wiener increments been passed
if nargin < 3
    flag = '';
end

switch flag
    case 'init'
        % Initialize persitent handle variables
        FIG_HANDLE = figure(gcf);
        AX_HANDLE = gca;
        
        % Set units to pixels to get width of axis, set back to default
        units = get(AX_HANDLE,'Units');
        set(AX_HANDLE,'Units','Pixels');
        pos = get(AX_HANDLE,'OuterPosition');
        set(AX_HANDLE,'Units',units);
        
        % Number of time samples to expect
        LEN_TSPAN = length(t);
        
        % Use figure axis width and TSPAN length determine redraw chunk
        chunk = min(ceil(LEN_TSPAN/pos(3)),LEN_TSPAN);
        
        % Number of solution variables, Y (cannot change)
        N = length(y);
        
        % Initialize UserData, T and Y
        ud = [];
        ud.t(1,chunk) = 0;
        ud.y(N,chunk) = 0;
        ud.i = 1;
        ud.t(1) = t(1);
        ud.y(:,1) = y;
        
        % Plot initial data
        if isW
            % Number of integrated Wiener increment variables, W (cannot change)
            D = length(w);
            
            % Initialize UserData, W
            ud.w(D,chunk) = 0;
            ud.w(:,1) = w;
            
            ud.lines = plot(ud.t(1),ud.y(:,1),'-',ud.t(1),ud.w(:,1),'--');
        else
            ud.lines = plot(ud.t(1),ud.y(:,1),'-');
        end
        
        % Set x-axis limits
        if ~ishold(AX_HANDLE)
            set(AX_HANDLE,'XLim',[min(t) max(t)]);
        end
        
        % Store UserData and draw
        set(FIG_HANDLE,'UserData',ud);
        drawnow;
    case ''
        if isempty(FIG_HANDLE) || isempty(AX_HANDLE)
            if isW
                error('SDETools:sdeplot:NotCalledWithInitW',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'',W0).']);
            else
                error('SDETools:sdeplot:NotCalledWithInit',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'').']);

            end
        end
        
        % If figure is open
        if ishghandle(FIG_HANDLE) && ishghandle(AX_HANDLE)
            % Get UserData
            ud = get(FIG_HANDLE,'UserData');
            lt = length(t);
            [N,lY] = size(ud.y);
            
            % Update UserData
            oldi = ud.i;
            newi = oldi+lt;
            if newi > lY
                % Set line data
                XYData = get(ud.lines,{'XData','YData'});
                for j = 1:N
                    XYData{j,1} = [XYData{j,1} ud.t];
                    XYData{j,2} = [XYData{j,2} ud.y(j,:)];
                end
                if isW
                    for j = N+1:N+size(ud.w,1)
                        XYData{j,1} = [XYData{j,1} ud.t];
                        XYData{j,2} = [XYData{j,2} ud.w(j-N,:)];
                    end
                end
                set(ud.lines,{'XData','YData'},XYData);
                
                % Set x-axis limits to auto if exceeded
                if ~ishold(AX_HANDLE)
                    if strcmp(get(AX_HANDLE,'XLimMode'),'manual')
                        XLim = get(AX_HANDLE,'XLim');
                        if min(ud.t) < XLim(1) || max(ud.t) > XLim(2)
                            set(AX_HANDLE,'XLimMode','auto');
                            LEN_TSPAN = max(length(XYData{1,1}),LEN_TSPAN);
                        end
                    else
                        LEN_TSPAN = max(length(XYData{1,1}),LEN_TSPAN);
                    end
                end
                
                % Check if figure width has changed
                units = get(AX_HANDLE,'Units');
                set(AX_HANDLE,'Units','Pixels');
                pos = get(AX_HANDLE,'OuterPosition');
                set(AX_HANDLE,'Units',units);
                
                % Use figure axis width and TSPAN length determine redraw chunk
                chunk = min(ceil(LEN_TSPAN/pos(3)),LEN_TSPAN);
                
                % Reset UserData
                ud.t(1,chunk) = 0;
                ud.y(:,chunk) = 0;
                if isW
                    ud.w(:,chunk) = 0;
                end
                oldi = 0;
                newi = lt;
            end
            
            % Append new data to UserData
            ud.t(oldi+1:newi) = t;
            ud.y(:,oldi+1:newi) = y;
            if isW
                ud.w(:,oldi+1:newi) = w;
            end
            ud.i = newi;
            
            % Store updated UserData and draw if redraw chunk was full
            set(FIG_HANDLE,'UserData',ud);
            if oldi == 0
                drawnow;
            end
        else
            status = 0;
        end
    case 'done'
        % Rename and delete persistent handles
        hf = FIG_HANDLE;
        FIG_HANDLE = [];
        ha = AX_HANDLE;
        AX_HANDLE = [];
        
        % If figure is open
        if ~isempty(hf) && ishghandle(hf) && ~isempty(ha) && ishghandle(ha)
            % Get non-zero UserData
            ud = get(hf,'UserData');
            ud.t = ud.t(1:ud.i);
            ud.y = ud.y(:,1:ud.i);
            N = size(ud.y,1);
            
            % Set any remaining line data
            XYData = get(ud.lines,{'XData','YData'});
            for j = 1:N
                XYData{j,1} = [XYData{j,1} ud.t];
                XYData{j,2} = [XYData{j,2} ud.y(j,:)];
            end
            if isW
                ud.w = ud.w(:,1:ud.i);
                for j = N+1:N+size(ud.w,1)
                    XYData{j,1} = [XYData{j,1} ud.t];
                    XYData{j,2} = [XYData{j,2} ud.w(j-N,:)];
                end
            end
            set(ud.lines,{'XData','YData'},XYData);
            
            % Delete UserData
            set(hf,'UserData',[]);
            
            % Refresh or draw
            if ishold(ha)
                drawnow;
            else
                % Set x-axis limits
                XLim = get(ha,'XLim');
                if min(XYData{1,1}) ~= XLim(1)
                    set(ha,'XLim',[min(XYData{1,1}) XLim(2)]);
                elseif max(XYData{1,1}) ~= XLim(2)
                    set(ha,'XLim',[XLim(1) max(XYData{1,1})]);
                end
                refresh;
            end
        else
            status = 0;
        end
    otherwise
        error('SDETools:sdeplot:InvalidFlag',...
              'Invalid status flag passed to output function.');
end