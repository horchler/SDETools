function status=sdephaseplot2(t,y,flag,w)
%SDEPHASEPLOT2  2-D phase plane SDE output function.
%   When the function SDEPHASEPLOT2 is passed to an SDE solver as the
%   'OutputFUN' property, i.e., OPTIONS = SDESET('OutputFUN',@SDEPHASEPLOT2),
%   the solver calls SDEPHASEPLOT2(T,Y,'') after every timestep. The
%   SDEPHASEPLOT2 function plots the first two components of the solution, Y, as
%   they are computed, adapting the axis limits of the plot dynamically. To plot
%   particular components of the solution, specify their indices via the
%   'OutputYSelect' SDE options property. The first specified component is
%   plotted with respect to the X-axis and the second component with respect to
%   the Y-axis.
%
%   At the start of integration, the solver calls SDEPHASEPLOT2(TSPAN,Y0,'init')
%   to initialize the output function. After each integration step to new time,
%   T, and solution vector, Y, the solver calls STATUS = SDEPHASEPLOT2(T,Y,'').
%   If the solver's 'Refine' property is greater than one (see SDESET), T is a
%   column vector of new output times and Y is an array of corresponding column
%   vectors. The STATUS return value is 1 if the figure window and plot axis are
%   still open and 0 otherwise. When the integration is complete, the solver
%   calls SDEPHASEPLOT2([],[],'done').
%
%   If fewer than two components of the solution, Y, are specified, a component
%   of the integrated Wiener increments, W, can be plotted versus a Y solution
%   component or two integrated Wiener increment components can be plotted. Set
%   the 'OutputWSelect' SDE options property to 'yes' or to a vector of indices
%   to pass the integrated Wiener increments to SDEPHASEPLOT2 as a fourth
%   argument, SDEPHASEPLOT2(T,Y,'',W). If the 'OutputWSelect' property is set to
%   'yes' or a non-empty vector of indices the solver calls
%   SDEPHASEPLOT2(TSPAN,Y0,'init',W0) to initialize the output function at the
%   start of integration and SDEPHASEPLOT2([],[],'done',[]) when the integration
%   is complete.
%   
%   See also:
%       SDEPLOT, SDEPHASEPLOT3, SDEIMGPLOT, SDEPRINT, SDESET, SDEGET, SDE_EULER,
%       SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU, ODEPHAS2

%   SDEPHASEPLOT2 is based on an updating of Matlab's ODEPHAS2, version 1.26.4.9

%   Andrew D. Horchler, horchler @ gmail . com, 5-11-13
%   Revision: 1.2, 5-13-13


persistent FIG_HANDLE AX_HANDLE;
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
        
        % Initialize UserData and Y
        ud = [];
        ud.y(2,chunk) = 0;
        ud.i = 1;
        
        if N >= 2
            ud.y(:,1) = y(1:2);
        elseif isW
            % Number of integrated Wiener increment variables, W (cannot change)
            D = length(w);
            
            if D+N >= 2
                if N == 1
                    ud.y(:,1) = [y(1);w(1)];
                else
                    ud.y(:,1) = w(1:2);
                end
            else
                error('SDETools:sdephaseplot2:TooFewInputsW',...
            	     ['Output function requires at least two solution '...
                      'components from Y or W.']);
            end
        else
            error('SDETools:sdephaseplot2:TooFewInputs',...
            	 ['Output function requires at least two solution '...
                  'components from Y.']);
        end
        
        % Plot initial data
        ud.lines = plot(ud.y(1),ud.y(2),'-');
        
        % Plot start and end markers
        nextplot = get(AX_HANDLE,'NextPlot');
        set(AX_HANDLE,'NextPlot','add');
        plot(AX_HANDLE,ud.y(1),ud.y(2),'g.');
        ud.marker = plot(AX_HANDLE,ud.y(1),ud.y(2),'r.');
        set(AX_HANDLE,'NextPlot',nextplot);
        
        % Store UserData and draw
        set(FIG_HANDLE,'UserData',ud);
        drawnow;
    case ''
        if isempty(FIG_HANDLE) || isempty(AX_HANDLE)
            if isW
                error('SDETools:sdephaseplot2:NotCalledWithInitW',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'',W0).']);
            else
                error('SDETools:sdephaseplot2:NotCalledWithInit',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'').']);

            end
        end
        
        % If figure is open
        if ishghandle(FIG_HANDLE) && ishghandle(AX_HANDLE)
            % Get UserData
            ud = get(FIG_HANDLE,'UserData');
            lt = length(t);
            N = size(y,1);
            lY = size(ud.y,2);
            
            % Update UserData
            oldi = ud.i;
            newi = oldi+lt;
            if newi > lY
                % Set line data
                XYData = get(ud.lines,{'XData','YData'});
                set(ud.lines,{'XData','YData'},{[XYData{1} ud.y(1,:)],...
                                                [XYData{2} ud.y(2,:)]});
                
                % Set marker point
                set(ud.marker,{'XData','YData'},{ud.y(1,end),ud.y(2,end)});
                
                % Reset UserData
                ud.y(:,lY) = 0;
                oldi = 0;
                newi = lt;
            end
            
            % Append new data to UserData
            if N >= 2
                ud.y(:,oldi+1:newi) = y(1:2,:);
            elseif N == 1
                ud.y(:,oldi+1:newi) = [y(1,:);w(1,:)];
            else
                ud.y(:,oldi+1:newi) = w(1:2,:);
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
            ud.y = ud.y(:,1:ud.i);
            
            % Set any remaining line data
            XYData = get(ud.lines,{'XData','YData'});
         	set(ud.lines,{'XData','YData'},{[XYData{1} ud.y(1,:)],...
                                            [XYData{2} ud.y(2,:)]});
            
            % Set final marker point
          	set(ud.marker,{'XData','YData'},{ud.y(1,end),ud.y(2,end)});
            
            % Delete UserData
            set(hf,'UserData',[]);
            
            % Refresh or draw
            if ishold(ha)
                drawnow;
            else
                refresh;
            end
        else
            status = 0;
        end
    otherwise
        error('SDETools:sdephaseplot2:InvalidFlag',...
              'Invalid status flag passed to output function.');
end