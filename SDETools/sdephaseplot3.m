function status=sdephaseplot3(t,y,flag,w)
%SDEPHASEPLOT3  3-D phase space SDE output function.
%   When the function SDEPHASEPLOT3 is passed to an SDE solver as the
%   'OutputFUN' property, i.e., OPTIONS = SDESET('OutputFUN',@SDEPHASEPLOT3),
%   the solver calls SDEPHASEPLOT3(T,Y,'') after every timestep. The
%   SDEPHASEPLOT3 function plots the first three components of the solution, Y,
%   as they are computed, adapting the axis limits of the plot dynamically. To
%   plot particular components of the solution, specify their indices via the
%   'OutputYSelect' SDE options property. The first specified component is
%   plotted with respect to the X-axis, the second component with respect to the
%   Y-axis, and the third component with respect to the Z-axis.
%
%   At the start of integration, the solver calls SDEPHASEPLOT3(TSPAN,Y0,'init')
%   to initialize the output function. After each integration step to new time,
%   T, and solution vector, Y, the solver calls STATUS = SDEPHASEPLOT3(T,Y,'').
%   If the solver's 'Refine' property is greater than one (see SDESET), T is a
%   column vector of new output times and Y is an array of corresponding column
%   vectors. The STATUS return value is 1 if the figure window and plot axis are
%   still open and 0 otherwise. When the integration is complete, the solver
%   calls SDEPHASEPLOT3([],[],'done').
%
%   If fewer than three components of the solution, Y, are specified, one or two
%   components of the integrated Wiener increments, W, can be plotted versus Y
%   solution components or three integrated Wiener increment components can be
%   plotted. Set the 'OutputWSelect' SDE options property to 'yes' or to a
%   vector of indices to pass the integrated Wiener increments to SDEPHASEPLOT3
%   as a fourth argument, SDEPHASEPLOT3(T,Y,'',W). If the 'OutputWSelect'
%   property is set to 'yes' or a non-empty vector of indices the solver calls
%   SDEPHASEPLOT3(TSPAN,Y0,'init',W0) to initialize the output function at the
%   start of integration and SDEPHASEPLOT3([],[],'done',[]) when the integration
%   is complete.
%   
%   See also:
%       SDEPLOT, SDEPHASEPLOT2, SDEIMGPLOT, SDEPRINT, SDESET, SDEGET, SDE_EULER,
%       SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU, ODEPHAS3

%   SDEPHASEPLOT3 is based on an updating of Matlab's ODEPHAS3,
%   version 1.27.4.10

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
        ud.y(3,chunk) = 0;
        ud.i = 1;
        
        if N >= 3
            ud.y(:,1) = y(1:3);
        elseif isW
            % Number of integrated Wiener increment variables, W (cannot change)
            D = length(w);
            
            if D+N >= 3
                if N == 2
                    ud.y(:,1) = [y(1:2);w(1)];
                elseif N == 1
                    ud.y(:,1) = [y(1);w(1:2)];
                else
                    ud.y(:,1) = w(1:3);
                end
            else
                error('SDETools:sdephaseplot3:TooFewInputsW',...
            	     ['Output function requires at least three solution '...
                      'components from Y or W.']);
            end
        else
            error('SDETools:sdephaseplot3:TooFewInputs',...
            	 ['Output function requires at least three solution '...
                  'components from Y.']);
        end
        
        % Plot initial data
        ud.lines = plot3(ud.y(1),ud.y(2),ud.y(3),'-');
        
        % Plot start and end markers, turn on grid
        nextplot = get(AX_HANDLE,'NextPlot');
        set(AX_HANDLE,'NextPlot','add');
        plot3(ud.y(1),ud.y(2),ud.y(3),'g.');
        ud.marker = plot3(ud.y(1),ud.y(2),ud.y(3),'r.');
        set(AX_HANDLE,'NextPlot',nextplot,'DrawMode','fast');
        grid on;
        
        % Store UserData and draw
        set(FIG_HANDLE,'UserData',ud);
        drawnow;
    case ''
        if isempty(FIG_HANDLE) || isempty(AX_HANDLE)
            if isW
                error('SDETools:sdephaseplot3:NotCalledWithInitW',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'',W0).']);
            else
                error('SDETools:sdephaseplot3:NotCalledWithInit',...
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
                XYZData = get(ud.lines,{'XData','YData','ZData'});
                set(ud.lines,{'XData','YData','ZData'},...
                             {[XYZData{1} ud.y(1,:)],...
                              [XYZData{2} ud.y(2,:)],...
                              [XYZData{3} ud.y(3,:)]});
                
                % Set marker point
                set(ud.marker,{'XData','YData','ZData'},...
                              {ud.y(1,end),ud.y(2,end),ud.y(3,end)});
                
                % Reset UserData
                ud.y(:,lY) = 0;
                oldi = 0;
                newi = lt;
            end
            
            % Append new data to UserData
            if N >= 3
                ud.y(:,oldi+1:newi) = y(1:3,:);
            elseif N == 2
                ud.y(:,oldi+1:newi) = [y(1:2,:);w(1,:)];
            elseif N == 1
                ud.y(:,oldi+1:newi) = [y(1,:);w(1:2,:)];
            else
                ud.y(:,oldi+1:newi) = w(1:3,:);
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
            XYZData = get(ud.lines,{'XData','YData','ZData'});
            set(ud.lines,{'XData','YData','ZData'},...
                         {[XYZData{1} ud.y(1,:)],...
                          [XYZData{2} ud.y(2,:)],...
                          [XYZData{3} ud.y(3,:)]});
            
            % Set final marker point
            set(ud.marker,{'XData','YData','ZData'},...
                          {ud.y(1,end),ud.y(2,end),ud.y(3,end)});
            
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
        error('SDETools:sdephaseplot3:InvalidFlag',...
              'Invalid status flag passed to output function.');
end