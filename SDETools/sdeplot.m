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
%   At the start of integration, the SDE solver calls SDEPLOT(TSPAN,Y0,'init')
%   to initialize the output function. After each integration step to a new time
%   point T and solution vector Y, the solver calls STATUS = SDEPLOT(T,Y,'').
%   If the solver's 'Refine' property is greater than one (see SDESET), T is a 
%   column vector of new output times and Y is an array of corresponding column
%   vectors. The STATUS return value is 1 if the figure window and plot axis are
%   still open and 0 otherwise. When the integration is complete, the solver
%   calls SDEPLOT([],[],'done').
%
%   Set the 'OutputWSelect' SDE options property to 'yes' or to a vector of
%   indices to output the integrated Wiener increments, W. The Wiener increments
%   are passed to SDEPLOT as a fourth argument, SDEPLOT(T,Y,'',W). If the
%   'OutputWSelect' property is set to 'yes' or a non-empty vector of indices
%   the SDE solver calls SDEPLOT(TSPAN,Y0,'init',W0) to initialize the output
%   function at the start of integration and SDEPLOT([],[],'done',[]) when the
%   integration is complete.
%   
%   See also:
%       SDESET, SDEGET, SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU,
%       ODEPLOT

%   SDEPLOT is based on an updating of version 1.25.4.9 of Matlab's ODEPLOT

%   Andrew D. Horchler, adh9 @ case . edu, 4-29-13
%   Revision: 1.2, 5-8-13


persistent FIG_HANDLE AX_HANDLE LEN_TSPAN;
status = 1;
isW = (nargin > 3);

switch flag
    case 'init'
        % Initialize persitent handle variables
        FIG_HANDLE = figure;
        AX_HANDLE = gca;
        
        % Set units to pixels to get width of axis, set back to default
        units = get(AX_HANDLE,'Units');
        set(AX_HANDLE,'Units','Pixels');
        pos = get(AX_HANDLE,'OuterPosition');
        set(AX_HANDLE,'Units',units);
        
        % Use figure axis width and TSPAN length determine redraw chunk
        LEN_TSPAN = 0.3*length(t);
        chunk = min(ceil(LEN_TSPAN/pos(3)),LEN_TSPAN);
        N = length(y);
        
        % Initialize UserData
        ud = [];
        ud.t(1,chunk) = 0;
        ud.y(N,chunk) = 0;
        ud.i = 1;
        ud.t(1) = t(1);
        ud.y(:,1) = y;
        
        % Plot initial data and set axis limits
        if isW
            D = length(w);
            ud.w(D,chunk) = 0;
            ud.w(:,1) = w;
            
            if ishold
                ud.lines = plot(ud.t(1),ud.y(:,1),'-',ud.t(1),ud.w(:,1),'--',...
                    'EraseMode','none');
            else
                ud.lines = plot(ud.t(1),ud.y(:,1),'-',ud.t(1),ud.y(:,1),'--');
                set(AX_HANDLE,'XLim',[min(t) max(t)]);  
            end
        else
            if ishold
                ud.lines = plot(ud.t(1),ud.y(:,1),'-','EraseMode','none');
            else
                ud.lines = plot(ud.t(1),ud.y(:,1),'-');
                set(AX_HANDLE,'XLim',[min(t) 0.3*max(t)]);  
            end
        end
        
        % Store UserData and draw
        set(FIG_HANDLE,'UserData',ud);
        drawnow;
    case ''
        if isempty(FIG_HANDLE)
            if nargin > 3
                error('SDETools:sdeplot:NotCalledWithInitW',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(tspan,y0,''init'',w).']);
            else
                error('SDETools:sdeplot:NotCalledWithInit',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(tspan,y0,''init'').']);

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
            if newi <= lY
                % Append new data to UserData
                ud.t(oldi+1:newi) = t;
                ud.y(:,oldi+1:newi) = y;
                if isW
                    ud.w(:,oldi+1:newi) = w;
                end
                ud.i = newi;
            else
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
                if ~ishold
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
                
                % Reset UserData and append new data
                ud.t(1,chunk) = 0;
                ud.y(:,chunk) = 0;
                ud.t(1:lt) = t;
                ud.y(:,1:lt) = y;
                if isW
                    ud.w(:,chunk) = 0;
                    ud.w(:,1:lt) = w;
                end
                ud.i = 1;
            end
            
            % Store updated UserData and draw if redraw chunk was full
            set(FIG_HANDLE,'UserData',ud);
            if newi > lY
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
        if ishghandle(hf) && ishghandle(ha)
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
            if ishold
                drawnow;
            else
                if ishghandle(ha)
                    set(ha,'XLimMode','auto');
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