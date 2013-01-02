function xi=sde_interpq(t,x,ti,eta)
%SDE_INTERPQ  Quick Brownian Bridge interpolation.
%   XI = SDE_INTERPQ(T,X,TI) interpolates to find XI, the values of the
%   underlying 1-D function X at the points of column vector TI. The column
%   vector T specifies the independent coordinates of the underlying interval on
%   which the rows of X are defined.
%
%   If X is a column vector, then XI is a column vector of LENGTH(TI).
%   
%   If X is a matrix, then the interpolation is performed for each column of X,
%   in which case XI is LENGTH(TI)-by-SIZE(X,2).
%
%   NaN is returned for any value of TI outside the range of T as well as for
%   any non-finite value even if within T.
%
%   SDE_INTERPQ is quicker than SDE_INTERP for non-uniformly spaced data because
%   it does no input checking and does not handle N-D array X inputs, matrix TI
%   inputs, or optional parameters. For SDE_INTERPQ to work properly:
%
%       T must be a monotonically increasing column vector. T must have at least
%           two elements and all values must be distinct, finite, and real.
%       X must be a real column vector or matrix with LENGTH(T) rows.
%       TI must be a real column vector and can be specified in any order.
%
%   Class support for inputs T, X, TI, ETA:
%       float: double, single
%
%   See also SDE_INTERP, SDE_INTERPNQ, INTERP1Q, INTERP1

% Some code based on version 1.15.4.4 of Matlab's INTERP1Q

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-28-12
%   Revision: 1.0, 1-1-13


if nargin == 3
    eta = 1;
    dtype=superiorfloat(t,x,ti);
else
    dtype=superiorfloat(t,x,ti,eta);
end
if isempty(ti)
	xi=zeros(0,size(x,2),dtype);
else
    [M,N]=size(x);
    D=length(ti);
    if D == 1
        % Find index of TI on subintervals of T
        i=find(ti >= t,1,'last');
        i(ti == t(end))=M-1;
        
        % Check if TI within bounds, it either is or is not
        if isempty(i) || i >= M
            xi=NaN(1,N,dtype);
        else
            % Mean of bridge values
            dti=ti-t(i);
            mu=dti/(t(i+1)-t(i));

            % Calculate XI, output will be 1-by-N row vector
            dx=diff(x);
            xi=x(i,:)+mu*dx(i,:)+eta(:).*sqrt((1-mu)*dti)*randn(1,N,dtype);
        end
    else
        % Allocate as NaN to set out of bounds values of TI
        xi=NaN(D,N,dtype);
        
        % Only two T values (one subinterval), TI values either in bounds or not
        if M == 2
            % Indices of TI on subintervals of T
            i=(ti >= t(1) & ti <= t(2));
            
            if any(i)
                % Mean and standard deviation of bridge values
                dti=ti(i)-t(1);
                ni=length(dti);
                mu=dti/(t(2)-t(1));
                sig=eta(:).*sqrt((1-mu).*dti);

                % Calculate XI, output is D-by-N matrix or D-by-1 column vector
                xi(i,:)=x(ones(ni,1),:)+mu*diff(x)+bsxfun(@times,randn(ni,N,dtype),sig);
            end
        else
            % Indices of TI on subintervals of T, sorted/bounded indices of T, X
            [tii,p]=sort(ti);
            i=(tii >= t(1) & tii <= t(end));
            if any(i)
                q=p(i);
                ti=ti(q);
                [~,ii]=histc(ti,t);
                ii(ti == t(end))=M-1;
                dt=diff(t);
                dt=dt(ii);
                dti=ti-t(ii);

                % Mean and standard deviation of bridge values
                mu=dti./dt;
                sig=eta(:).*sqrt((1-mu).*dti);

                % Expand and re-index X, dX, to match size of re-indexed TI, XI
                dx=diff(x);
                x=x(bsxfun(@plus,ii,0:M:M*(N-1)));
                dx=dx(bsxfun(@plus,ii,0:(M-1):(M-1)*(N-1)));
                
                % Calculate XI, use sort indices for TI to reduce dimensions
                xi(q,:)=x+bsxfun(@times,dx,mu)+bsxfun(@times,randn(length(dti),N,dtype),sig);
            end
        end
    end
end