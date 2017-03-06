function status=sdeprint(t,y,flag,w)
%SDEPRINT  Command window printing SDE output function.
%   When the function SDEPRINT is passed to an SDE solver as the 'OutputFUN'
%   property, i.e., OPTIONS = SDESET('OutputFUN',@SDEPRINT), the solver calls 
%   SDEPRINT(T,Y,'') after every timestep. The SDEPRINT function prints the
%   components of the solution, Y, as they are computed. If more than 16
%   components are specified, an abreviated output is given. To print only
%   particular components of the solution, specify their indices via the
%   'OutputYSelect' SDE options property.
%
%   At the start of integration, the solver calls SDEPRINT(TSPAN,Y0,'init') to
%   initialize the output function. After each integration step to new time, T,
%   and solution vector, Y, the solver calls STATUS = SDEPRINT(T,Y,''). If the
%   solver's 'Refine' property is greater than one (see SDESET), T is the most
%   recent output time and Y is the corresponding column vector. The STATUS
%   return value is 1. When the integration is complete, the solver calls
%   SDEPRINT([],[],'done').
%
%   Set the 'OutputWSelect' SDE options property to 'yes' or to a vector of
%   indices to output the integrated Wiener increments, W. The integrated Wiener
%   increments are passed to SDEPRINT as a fourth argument, SDEPRINT(T,Y,'',W).
%   If the 'OutputWSelect' property is set to 'yes' or a non-empty vector of
%   indices the solver calls SDEPRINT(TSPAN,Y0,'init',W0) to initialize the
%   output function at the start of integration and SDEPRINT([],[],'done',[])
%   when the integration is complete.
%
%   The display format of the printed output can be changed via FORMAT.
%   
%   See also:
%       SDEPLOT, SDEPHASEPLOT2, SDEPHASEPLOT3, SDEIMGPLOT, SDESET, SDEGET,
%       SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU, ODEPRINT, FORMAT

%   SDEPRINT is based on an updating of Matlab's ODEPRINT, version 1.17.4.4

%   Andrew D. Horchler, horchler @ gmail . com, 5-11-13
%   Revision: 1.2, 5-13-13


persistent LEN_TSPAN ITERATION STR BSSTR FSTR C;	%#ok<PUSE>
status = 1;         % Figure widow still open and and plot axis still present
isW = (nargin > 3); % Have integrated Wiener increments been passed
if nargin < 3
    flag = '';
end
bs = char(8);       % '\b', backspace character to overwrite previous line
id = '|';           % Progress bar indicator

switch flag
    case 'init'
        % Number of time samples to expect
        LEN_TSPAN = length(t);
        
        % Number of solution variables, Y (cannot change)
        N = length(y);
        
        % Allocate blank array for output text, char(10) is '\n'
        sp = ' ';
        ln = [sp(ones(60,1));10];
        if isW
            % Number of integrated Wiener increment variables, W (cannot change)
            D = length(w);
            
            STR = ln(:,ones(1,min(max(N,D),17)+5));
        else
        	STR = ln(:,ones(1,min(N,17)+5));
        end
        
        % Get output format
        isDouble = isa(y,'double');
        switch get(0,'format')
            case 'short'
                FSTR = '% -.5g';
            case 'long'
                if isDouble
                    FSTR = '% -.15g';
                else
                    FSTR = '% -.7g';
                end
            case 'shorte'
                FSTR = '% -.5d';
            case 'longe'
                if isDouble
                    FSTR = '% -.15d';
                else
                    FSTR = '% -.7d';
                end
            case 'shortg'
                FSTR = '% -.5g';
            case 'longg'
                if isDouble
                    FSTR = '% -.15g';
                else
                    FSTR = '% -.7g';
                end
            case 'shorteng'
                FSTR = '% -.5e';
            case 'longeng'
                if isDouble
                    FSTR = '% -.15e';
                else
                    FSTR = '% -.7e';
                end
            otherwise
                FSTR = '% -.15g';
        end
        
        strout = STR;
        
        % Progress bar
        strout(1,2) = '[';
        strout(end-9:end-2,2) = ']   0.0%';
        
        % Time
        tout = sprintf([' T = ' FSTR],t(1));
        strout(1:length(tout),4) = tout;
        
        % Step number
        ITERATION = 0;
        tout = sprintf(' Step = %d',ITERATION);
        strout(31:30+length(tout),4) = tout;
        
        % Divider
        strout(1:end-2,5) = '-';
        
        % Y
        if N >= 1
            yout = sprintf([' Y = ' FSTR],y(1));
            strout(1:length(yout),6) = yout;

            if N > 1
                ystr = ['     ' FSTR];
                if N > 16
                    for i = 2:15
                        yout = sprintf(ystr,y(i));
                        strout(1:length(yout),5+i) = yout;
                    end
                    strout(1:9,21) = '      ...';
                    yout = sprintf(ystr,y(end));
                 	strout(1:length(yout),22) = yout;
                else
                    for i = 2:N
                        yout = sprintf(ystr,y(i));
                        strout(1:length(yout),5+i) = yout;
                    end
                end
            end
            
            offset = 30;
        else
            offset = 0;
        end
        
        % W
        if isW && D >= 1
            wout = sprintf([' W = ' FSTR],w(1));
            strout(offset+1:offset+length(wout),6) = wout;

            if D > 1
                wstr = ['     ' FSTR];
                if D > 16
                    for i = 2:15
                        wout = sprintf(wstr,w(i));
                        strout(offset+1:offset+length(wout),5+i) = wout;
                    end
                    strout(offset+1:offset+9,21) = '      ...';
                    wout = sprintf(wstr,w(end));
                  	strout(offset+1:offset+length(wout),22) = wout;
                else
                    for i = 2:D
                        wout = sprintf(wstr,w(i));
                        strout(offset+1:offset+length(wout),5+i) = wout;
                    end
                end
            end
        end
        
        % Print text, pass as string to handle '%' and '\n' characters
        nchars = fprintf(1,'%s',strout(:));
        
        % Allocate array of backspace characters to be used to overwrite text
        BSSTR = bs(ones(nchars,1));
        
        % Get and set state of lastwarn so warning can be caught
        [lastmsg,lastid] = lastwarn('','SDETools:sdeprint:CatchWarning');
        
        % Persistent cleanup function to ensure lastwarn reset to original state
        C = onCleanup(@()lastwarn(lastmsg,lastid));
    case ''
        if isempty(LEN_TSPAN) || isempty(STR)
            if isW
                error('SDETools:sdeprint:NotCalledWithInitW',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'',W0).']);
            else
                error('SDETools:sdeprint:NotCalledWithInit',...
                     ['Output function has not been initialized. Use syntax '...
                      'OutputFUN(TSPAN,Y0,''init'').']);

            end
        end
        
        % Blank out previous data with spaces
        strout = STR;
        
        % Progress bar
        ITERATION = ITERATION+1;
        pct = ITERATION/LEN_TSPAN;
        idout = ['[';id(ones(min(round(51*pct),51),1))];
        strout(1:length(idout),2) = idout;
        strout(end-9:end-2,2) = ']      %';
        pctout = sprintf('%.1f',100*pct);
        strout(end-2-length(pctout):end-3,2) = pctout;
        
        % Time
        tout = sprintf([' T = ' FSTR],t(end));
        strout(1:length(tout),4) = tout;
        
        % Step number
        tout = sprintf(' Step = %d',ITERATION);
        strout(31:30+length(tout),4) = tout;
        
        % Divider
        strout(1:end-2,5) = '-';
        
        % Y
        N = length(y);
        if N >= 1
            yout = sprintf([' Y = ' FSTR],y(1,end));
            strout(1:length(yout),6) = yout;
            
            if N > 1
                ystr = ['     ' FSTR];
                if N > 16
                    for i = 2:15
                        yout = sprintf(ystr,y(i,end));
                        strout(1:length(yout),5+i) = yout;
                    end
                    strout(1:9,21) = '      ...';
                    yout = sprintf(ystr,y(end,end));
                 	strout(1:length(yout),22) = yout;
                else
                    for i = 2:N
                        yout = sprintf(ystr,y(i,end));
                        strout(1:length(yout),5+i) = yout;
                    end
                end
            end
            
            offset = 30;
        else
            offset = 0;
        end
        
        % W
        if isW
            D = length(w);
            if D >= 1
                wout = sprintf([' W = ' FSTR],w(1,end));
                strout(offset+1:offset+length(wout),6) = wout;

                if D > 1
                    wstr = ['     ' FSTR];
                    if D > 16
                        for i = 2:15
                            wout = sprintf(wstr,w(i,end));
                            strout(offset+1:offset+length(wout),5+i) = wout;
                        end
                        strout(offset+1:offset+9,21) = '      ...';
                        wout = sprintf(wstr,w(end,end));
                     	strout(offset+1:offset+length(wout),22) = wout;
                    else
                        for i = 2:D
                            wout = sprintf(wstr,w(i,end));
                            strout(offset+1:offset+length(wout),5+i) = wout;
                        end
                    end
                end
            end
        end
        
        % Catch warnings so they are not overwritten
        [msg,msgid] = lastwarn;	%#ok<*ASGLU>
        if ~strcmp(msgid,'SDETools:sdeprint:CatchWarning')
            % Reset state of lastwarn
            lastwarn('','SDETools:sdeprint:CatchWarning');
            
            % Empty backspace array so warning(s) not overwritten
            bsout = [];
        else
            % Backspace array to overwrite previous text
            bsout = BSSTR;
        end
        
        % Simultaneously overwrite (backspace) previous text and print new text
        fprintf(1,'%s',[bsout;strout(:)]);
    case 'done'
        if ~isempty(LEN_TSPAN) && ~isempty(STR)
            fprintf(1,'\n');
        end
    otherwise
        error('SDETools:sdeprint:InvalidFlag',...
              'Invalid status flag passed to output function.');
end