function status=sdeplot(t,y,flag,w)
%SDEPLOT  Time series SDE output function.
%
%   
%   See also:
%       SDESET, SDEGET, SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU

%   SDEPLOT is based on an updating of version 1.25.4.9 of Matlab's ODEPLOT

%   Andrew D. Horchler, adh9 @ case . edu, 4-29-13
%   Revision: 1.2, 5-3-13


persistent fig_handle ax_handle;
status = 0;
chunk = 128;
isW = false;    % Disabled for now

switch flag
    case 'init'
        % Initialize persitent handle variables
        fig_handle = figure;
        ax_handle = gca;
        
        N = length(y);
        
        % Initialize UserData
        ud = [];
        ud.t(chunk,1) = 0;
        ud.y(chunk,N) = 0;
        ud.i = 1;
        ud.t(1) = t(1);
        ud.y(1,:) = y;
        
        if isW
            D = length(w);
            ud.w(chunk,D) = 0;
            ud.w(1,:) = w;
        end
        
        if ishold
            if isW
                ud.lines = plot(ud.t(1),ud.y(1,:),'-',ud.t(1),ud.w(1,:),'--',...
                    'EraseMode','none');
            else
                ud.lines = plot(ud.t(1),ud.y(1,:),'-','EraseMode','none');
            end
        else
            if isW
                ud.lines = plot(ud.t(1),ud.y(1,:),'-',ud.t(1),ud.w(1,:),'--');
            else
                ud.lines = plot(ud.t(1),ud.y(1,:),'-');
            end
            set(ax_handle,'XLim',[min(t) max(t)]);  
        end
        
        set(fig_handle,'UserData',ud);
        drawnow;
    case ''
        if isempty(fig_handle)
            if nargin > 3
                error('SHCTools:sdeplot:NotCalledWithInitW',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(tspan,y0,''init'',w).']);
            else
                error('SHCTools:sdeplot:NotCalledWithInit',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(tspan,y0,''init'').']);

            end
        end
        
        % If figure is open
        if ishghandle(fig_handle) && ishghandle(ax_handle)
            ud = get(fig_handle,'UserData');
            nt = length(t);
            chunk = max(chunk,nt);
            
            [ltY,N] = size(ud.y);
            
            oldi = ud.i;
            newi = oldi+nt;
            if newi > ltY
                ud.t = [ud.t;zeros(chunk,1)];
                ud.y = [ud.y;zeros(chunk,N)];
            end
            ud.t(oldi+1:newi) = t;
            ud.y(oldi+1:newi,:) = y;
            
            if isW
                [ltW,D] = size(ud.w);
                if newi > ltW
                    ud.w = [ud.w;zeros(chunk,D)];
                end
                ud.w(oldi+1:newi,:) = w;
            end
            
            ud.i = newi;
            
            set(fig_handle,'UserData',ud);
            
            for j = 1:N
                set(ud.lines(j),'Xdata',ud.t(1:newi),'Ydata',ud.y(1:newi,j));
            end
            drawnow;
        end
    case 'done'
        % Rename and delete persistent variables
        hf = fig_handle;
        fig_handle = [];
        ha = ax_handle;
        ax_handle = [];
        
        % If figure is open
        if ishghandle(hf) && ishghandle(ha)
            ud = get(hf,'UserData');
            ud.t = ud.t(1:ud.i);
            ud.y = ud.y(1:ud.i,:);
            set(hf,'UserData',ud);
            for j = 1:size(ud.y,2)
                set(ud.lines(j),'Xdata',ud.t,'Ydata',ud.y(:,j));
            end
            if ~ishold
                if ishghandle(ha)
                    set(ha,'XLimMode','auto');
                end
                refresh;
            end
        end
    otherwise
        error('SHCTools:sdeplot:InvalidFlag',...
              'Invalid status flag passed to output function.');
end