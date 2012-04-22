function xi=sde_interpqn(t,x,nt)
%SDE_INTERPQN  Quick Linearly-spaced Brownian Bridge interpolation.
%   XI = SDE_INTERPQN(T,X,NT) interpolates to find XI, the values of the
%   underlying 1-D function X at the points of the vector T and NT
%   linearly-spaced points in each subinterval of T. The column vector T
%   specifies the independent coordinates of the underlying interval on which
%   the rows of X are defined.
%
%   If X is an M element column vector, then XI is a column vector of length 
%   NT*(M-1)+M.
%   
%   If X is an M-by-N matrix, then the interpolation is performed for each of
%   the columns of X, in which case XI is (NT*(M-1)+M)-by-N.
%
%   SDE_INTERPQN is quicker than SDE_INTERP for non-uniformly spaced data
%   because it does no input checking and does not handle N-D array X inputs or
%   optional parameters. For SDE_INTERPQN to work properly:
%
%       T must be a monotonically increasing column vector. T must have at least
%           two elements and all values must be distinct.
%       X must be a real column vector or matrix with M = LENGTH(T) rows.
%       NT must be a finite real positive scalar.
%
%   Class support for inputs T, X, NT:
%       float: double, single
%
%   See also SDE_INTERP, SDE_INTERPQ, INTERP1Q, INTERP1

% Some code based on version 1.15.4.4 of Matlab's INTERP1Q

%   Andrew D. Horchler, adh9@case.edu, Created 2-29-12
%   Revision: 1.0, 4-21-12


dtype=superiorfloat(t,x,nt);
if isempty(nt) || nt < 0
	xi=zeros(0,size(x,2),dtype);
else
    %j=~isinf(t);
    %x=x(j,:);
    [M N]=size(x);
    %{
    if M == 0
        % Calculate XI, output is a 1-by-N vector or a scalar
        xi=NaN(1,N,dtype);
    elseif M == 1
        % Calculate XI, output is a 1-by-N vector or a scalar
        xi=cast(x,dtype);
    else
    %}
    j=isinf(t);
    if any(j)
        D=M*nt-nt+M;
        xi=NaN(D,N,dtype);
        j=~j;
        x=x(j,:);
        M=size(x,1);
        t=t(j);
    else
        j=~j;
    end
    
    % Mean of bridge values
    mu=0:1/(nt+1):1;
    mu=mu(:);

    if M == 1
        % Calculate XI, output is a 1-by-N vector or a scalar
        xi(j,:)=x;
    else
        if M == 2
            % Expanded indices
            D=nt+2;
            i=ones(D,1);
            mu=mu(i(j));
            % Standard deviation of bridge values
            sig=sqrt((1-mu).*mu*(t(2)-t(1)));
        else
            % Expanded indices
            D=M*nt-nt+M;
            i=1:M-1;
            i=i(ones(1,nt+1),:);
            i=[i(:);M-1];

            % Expanded mean of bridge values
            q=1:nt+1;
            q=q(:);
            mu=[mu(q(:,ones(1,M-1)),1);1];

            % Standard deviation of bridge values
            dt=diff(t);
            sig=sqrt((1-mu).*mu.*dt(i));
        end
        i=i(j);
        dx=diff(x);
        
        xi(j,:)
        x(i,:)
        bsxfun(@times,dx(i,:),mu)
        bsxfun(@times,randn(D,N,dtype),sig)
        % Calculate XI, output is a D-by-N matrix or D-by-1 column vector
        xi(j,:)=x(i,:)+bsxfun(@times,dx(i,:),mu)+bsxfun(@times,randn(D,N,dtype),sig);
    end
end