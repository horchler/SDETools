function [te,ye,ie,vnew,stop] = sdezero(EventsFUN,t,y,value,args)
%SDEZERO  Locate any zero-crossings of event functions in a time step.
%
%   See also:
%       SDE_EULER, SDE_MILSTEIN, SDEARGUMENTS, SDEGET, SDESET, FUNCTION_HANDLE
        
%   Andrew D. Horchler, adh9 @ case . edu, Created 12-30-11
%   Revision: 1.0, 12-31-12

%   SDEZERO is loosely based on Matlab's ODEZERO helper function.


[vnew,isterminal,direction] = feval(EventsFUN,t,y,args{:});
if isempty(direction)
    direction = 0;
end
z = sign(vnew).*sign(value) <= 0 & (direction == 0 | direction == 1 ...
    & vnew > value | direction == -1 & vnew < value);
if any(z)
    ie = find(z(:));
    q = zeros(length(ie),1);
    te = t+q;
    y = y';
    ye = y(1+q,:);
    stop = any(isterminal == 1);
else
    te = [];
    ye = [];
    ie = [];
    stop = false;
end