function [te,ye,we,ie,vnew,stop] = sdezero(EventsFUN,t,y,w,value)
%SDEZERO  Locate any zero-crossings of event functions in a time step.
%
%   See also:
%       SDE_EULER, SDE_MILSTEIN, SDE_BM, SDE_GBM, SDE_OU, SDEARGUMENTS,
%       SDEARGUMENTS_PROCESS, SDEGET, SDESET, FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-30-11
%   Revision: 1.2, 5-4-13

%   SDEZERO is loosely based on Matlab's ODEZERO helper function.


[vnew,isterminal,direction] = EventsFUN(t,y(:));
if isempty(direction)
    direction = 0;
end
z = sign(vnew).*sign(value) < 0 & (direction == 0 | direction == 1 ...
    & vnew > value | direction == -1 & vnew < value);
if any(z)
    ie = find(z(:));
    q = ones(length(ie),1);
    te = t+q-1;
    y = y(:).';
    ye = y(q,:);
    w = w(:).';
    we = w(q,:);
    stop = any(isterminal & z);
else
    te = [];
    ye = [];
    we = [];
    ie = [];
    stop = false;
end