function xi=sde_interp(varargin)
%SDE_INTERP  Brownian Bridge interpolation.
%   XI = SDE_INTERP(T,X,TI) interpolates to find XI, the values of the 
%   underlying function X at the points in the array TI. T must be a vector of
%   length M.
%
%   If X is a vector, then it must also have length M, and XI will be the same
%   size as TI. If X is an array of size [M,D1,D2,...,Dk], then the
%   interpolation is performed for each D1-by-D2-by-...-Dk value in
%   X(i,:,:,...,:).
%
%   If TI is a vector of length N, then XI has size [N,D1,D2,...,Dk].
%
%   If TI is an array of size [N1,N2,...,Nj], then XI is of size
%   [N1,N2,...,Nj,D1,D2,...,Dk].
%
%   XI = SDE_INTERP(X,TI) assumes T = 0:M-1, where M = LENGTH(X) for vector X
%   or M = SIZE(X,1) for array X.
%
%   XI = SDE_INTERP(T,X,'linear',NT) creates TI internally ...... returns NT
%   linearly-spaced interpolated values betwen each subinterval of T. NT must be
%   a finite real positive scalar.
%
%   XI = SDE_INTERP(T,X,TI,'extrap',EXTRAPVAL) replaces values outside of the
%   interval spanned by T with the scalar EXTRAPVAL instead of the default NaN.
%
%   Class support for inputs T, X, TI, NT, ETA:
%       float: double, single
%
%   See also SDE_INTERPQ, SDE_INTERPNQ, INTERP1, INTERP1Q

%	Some code based on version 5.41.4.18 of Matlab's INTERP1

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-26-12
%   Revision: 1.0, 1-1-13


%XI = SDE_INTERP(T,X,TI)
%XI = SDE_INTERP(T,X,TI)
%XI = SDE_INTERP(T,X,TI,'extrap',EXTRAPVAL)
%XI = SDE_INTERP(T,X,TI,'scaling',ETA)
%XI = SDE_INTERP(T,X,TI,'extrap',EXTRAPVAL,'scaling',ETA)
%XI = SDE_INTERP(T,X,TI,'scaling',ETA,'extrap',EXTRAPVAL)
%XI = SDE_INTERP(X,TI)
%XI = SDE_INTERP(X,TI,'extrap',EXTRAPVAL)
%XI = SDE_INTERP(X,TI,'scaling',ETA)
%XI = SDE_INTERP(X,TI,'extrap',EXTRAPVAL,'scaling',ETA)
%XI = SDE_INTERP(X,TI,'scaling',ETA,'extrap',EXTRAPVAL)

%XI = SDE_INTERP(T,X,'linear',NT)
%XI = SDE_INTERP(T,X,'linear',NT,'scaling',ETA)
%XI = SDE_INTERP(X,'linear',NT)
%XI = SDE_INTERP(X,'linear',NT,'scaling',ETA)

if nargin < 2
    error('SDETools:sde_interp:TooFewInputs','Not enough input arguments.');
end
if nargin > 7
    error('SDETools:sde_interp:TooManyInputs','Too many input arguments.');
end

% Read variable inputs
eqsp=false;
extrap=false;
linspc=false;
scaling=false;
if nargin == 2
    x=varargin{1};
	ti=varargin{2};
    eqsp=true;
else
    for offset=2:nargin
        if ischar(varargin{offset})
            break;
        end
    end
    if offset == 2
        x=varargin{1};
        if ~strcmp(varargin{2},'linear')
            error('SDETools:sde_interp:InvalidLinearString',...
                 ['The parameter string ''linear'' must be used to specify '...
                  'the number of linearly-spaced interpolation points, NT, '...
                  'in the subsequent input argument.']);
        end
        nt=varargin{3};
        eqsp=true;
        linspc=true;
    elseif offset == 3
        if nargin == 3
            t=varargin{1};
            x=varargin{2};
            ti=varargin{3};
        else
            if nargin < 4
                error('SDETools:sde_interp:PropertNameValueMismatch4',...
              	     ['Property name string inputs must be followed by '...
                      'property value inputs.']);
            end
            str={'extrap','scaling','linear'};
            s=strcmp(varargin{3},str);
            if s(3)
                t=varargin{1};
                x=varargin{2};
                nt=varargin{4};
                linspc=true;
                if nargin > 4
                    if nargin ~= 6
                        error('SDETools:sde_interp:PropertNameValueMismatch46a',...
                             ['Property name string inputs must be followed '...
                              'by property value inputs.']);
                    end
                    if ~strcmp(varargin{5},str{2})
                        error('SDETools:sde_interp:DoublePropertyNameScaling',...
                             ['The last property name has already been '...
                              'specified.']);
                    end
                    eta=varargin{6};
                    scaling=true;
                end
            elseif s(1) || s(2)
                x=varargin{1};
                ti=varargin{2};
                if s(1)
                    extrapval=varargin{4};
                    extrap=true;
                else
                    eta=varargin{4};
                    scaling=true;
                end
                eqsp=true;
                if nargin > 4
                    if nargin ~= 6
                        error('SDETools:sde_interp:PropertNameValueMismatch46b',...
                             ['Property name string inputs must be followed '...
                              'by property value inputs.']);
                    end
                    if ~strcmp(varargin{5},str{~s(1:2)})
                        error('SDETools:sde_interp:DoublePropertyName5',...
                             ['The last property name has already been '...
                              'specified.']);
                    end
                    if s(2)
                        extrapval=varargin{6};
                        extrap=true;
                    else
                        eta=varargin{6};
                        scaling=true;
                    end
                end
            else 
                error('SDETools:sde_interp:a',...
                      '');
            end
        end
    elseif offset == 4
        if nargin < 5
            error('SDETools:sde_interp:PropertNameValueMismatch45',...
              	 ['Property name string inputs must be followed by property '...
                   'value inputs.']);
        end
        t=varargin{1};
        x=varargin{2};
        ti=varargin{3};
        str={'extrap','scaling'};
        s=strcmp(varargin{4},str);
        if s(1)
            extrapval=varargin{5};
            extrap=true;
        else
            eta=varargin{5};
            scaling=true;
        end
        if nargin > 5
            if nargin ~= 7
                error('SDETools:sde_interp:PropertNameValueMismatch57',...
              	     ['Property name string inputs must be followed by '...
                      'property value inputs.']);
            end
            if ~strcmp(varargin{6},str{~s})
                error('SDETools:sde_interp:DoublePropertyName6',...
                      'The last property name has already been specified.');
            end
            if s(2)
                extrapval=varargin{7};
                extrap=true;
            else
                eta=varargin{7};
                scaling=true;
            end
        end
    else
        error('SDETools:sde_interp:InvalidArguments',...
              'Unknown input arguments.');
    end
end

% Check X
if isempty(x) || ~isfloat(x)
    error('SDETools:sde_interp:XEmptyOrNotFloat',...
          'The input X must be a non-empty floating-point array.');
end
szx=size(x);
isVecX=(isvector(x) && ~isempty(x));
if isVecX
    M=length(x);
    N=1;
    NN=1;
    x=x(:);
else
    M=szx(1);
    N=szx(2:end);
    NN=prod(N);
    x=reshape(x,[M NN]);
end
if M < 2
    error('SDETools:sde_interp:XNotLongEnough',...
	      'The input X should contain at least two points.');
end

% Construct effective T vector if not specified, otherwise check T
if eqsp
    t=[0;M-1];
    dt=1;
else
    if ~isvector(t)
        error('SDETools:sde_interp:TNotVector','The input T must be a vector.');
    end
    if length(t) ~= M
        if N == 1
            error('SDETools:sde_interp:XVectorDimensionMismatch',...
                 ['If the input X is a vector, it must have the same length '...
                  'as the input T.']);
        else
            error('SDETools:sde_interp:XArrayDimensionMismatch',...
                 ['If the input X is an array, size(X,1) must equal the '...
                  'length of the input T.']);
        end
    end
    if ~isfloat(t) || ~isreal(t)
        error('SDETools:sde_interp:InvalidTDataType',...
              'Datatype of input vector T must be real single or real double.');
    end
    if any(~isfinite(t))
        error('SDETools:sde_interp:TNotNaN',...
              'One or more elements of the input T are not finite.');
    end
    t=t(:);
    dt=diff(t);
    
    % Sort T and X if T is not strictly ascending
    if all(dt < 0)
        p=M:-1:1;
        t=t(p);
        x=x(p,:);
        dt=-dt;
    elseif any(dt < 0)
        [t,p]=sort(t);
        x=x(p,:);
        dt=diff(t);
    end
    
    if any(dt == 0)
        error('SDETools:sde_interp:TNotDistinct',...
              'The entries in the input vector T must be distinct.');
    end
    
    % Determine if T is effectively equispaced in numerical sense
    eqsp=(M == 2 || norm(diff(dt),Inf) <= eps(max(abs(t([1 end])))));
    if eqsp
        dt=(t(end)-t(1))/(M-1);
    end
end

% Check extrapolation value if present
if extrap && (~isscalar(extrapval) || ~isfloat(extrapval))
    error('SDETools:sde_interp:InvalidExtrapval',...
         ['The optional extrapolation value parameter must be a '...
          'floating-point scalar.']);
end

% Check scaling parameter if present
if scaling
    if linspc
        if ~isvector(eta) || ~isfloat(eta)
            error('SDETools:sde_interp:InvalidScaling',...
                 ['The optional scaling parameter must be a non-empty '...
                  'floating-point scalar or vector equal in length to TI.']);
        end
    else
        if ~isvector(eta) || ~isfloat(eta) || ~any(length(eta) == [N 1])
            error('SDETools:sde_interp:InvalidScaling',...
                 ['The optional scaling parameter must be a non-empty '...
                  'floating-point scalar or vector equal in length to TI.']);
        end
    end
    eta = eta(:);
else
    eta=1;
end

% Construct linearly-spaced TI from NT, otherwise check TI
if linspc
    dtype=superiorfloat(t,x,nt,eta);
    if isempty(nt)
        xi=zeros(0,NN,dtype);
        D=0;
    else
        % Mean of bridge values
        mu=0:1/(nt+1):1;
        mu=mu(:);
        if M == 2
            % Standard deviation of bridge values
            sig=sqrt((1-mu).*mu*dt);
            
            % Expanded indices
            D=nt+2;
            ii=ones(D,1);
        else
            % Expanded indices
            D=M*nt-nt+M;
            ii=1:M-1;
            ii=ii(ones(1,nt+1),:);
            ii=[ii(:);M-1];
            
            % Expanded mean of bridge values
            q=1:nt+1;
            q=q(:);
            mu=[mu(q(:,ones(1,M-1)),1);1];
            
            % Standard deviation of bridge values
            if eqsp
                sig=sqrt((1-mu).*mu*dt);
            else
                sig=sqrt((1-mu).*mu.*dt(ii));
            end
        end
        
        dx=diff(x);
        if isVecX
            % Calculate XI, output is a D-by-1 column vector
            xi=x(ii)+mu.*dx(ii)+eta.*sig.*randn(D,1,dtype);
        else
            % Expand and re-index X, and dX, to match size of re-indexed XI
            x=x(bsxfun(@plus,ii,0:M:M*(NN-1)));
            dx=dx(bsxfun(@plus,ii,0:(M-1):(M-1)*(NN-1)));

            % Calculate XI, output is a D-by-N matrix
            if isscalar(eta)
                xi=x+bsxfun(@times,dx,mu)+bsxfun(@times,randn(D,NN,dtype),eta.*sig);
            else
                xi=x+bsxfun(@times,dx,mu)+(sig*eta').*randn(D,NN,dtype);
            end
        end
    end
    szti=[D 1];
else
    if ~isfloat(ti) || ~isreal(ti)
        error('SDETools:sde_interp:InvalidTiDataType',...
              'Datatype of input Ti must be real single or real double.');
    end
    szti=size(ti);
    ti=ti(:);
    D=length(ti);
    
    % Allocate out of bounds values as NaN, or override default extrapolation
    dtype=superiorfloat(t,x,ti,eta);
    if extrap && ~isnan(extrapval)
        switch extrapval
            case -Inf
                xi=-Inf(D,NN,dtype);
            case -1
                xi=-ones(D,NN,dtype);
            case 0
                if isempty(D)
                    xi=zeros(D,NN,dtype);
                elseif strcmp(dtype,'double');
                    xi(D,NN)=0;
                else
                    xi(D,NN)=single(0);
                end
            case 1
                xi=ones(D,NN,dtype);
            case Inf
                xi=Inf(D,NN,dtype);
            otherwise
                xi=zeros(D,NN,dtype)+extrapval;
        end
    else
        xi=NaN(D,NN,dtype);
    end
    
    % Sort TI if not in monotonically increasing order, find in-bounds indices
    if ~eqsp && ~issorted(ti)
        [~,p]=sort(ti);
        p=p(~(ti(p) < t(1) | ti(p) > t(end)));
    else
        p=1:D;
        p=p(~(ti < t(1) | ti > t(end)));
    end
    
    lp=length(p);
    if lp ~= 0
        % Re-index TI and find indices on subintervals of T
        ti=ti(p);
        if offset == 0
            ii=min(max(1+floor(ti),1),M-1);
            dti=ti-ii+1;
        elseif eqsp
            ii=min(max(1+floor((ti-t(1))/dt),1),M-1);
            dti=ti-t(ii);
        else
            [~,ii]=histc(ti,t);
            ii(end)=M-1;
            dt=dt(ii);
            dti=ti-t(ii);
        end
        
        % Mean and standard deviation of bridge values
        mu=dti./dt;
        sig=sqrt((1-mu).*dti);
        
        dx=diff(x);
        if isVecX
            % Calculate XI, use sort indices for TI to reduce dimensions
            xi(p)=x(ii)+mu.*dx(ii)+eta.*sig.*randn(lp,1,dtype);
        else
            % Expand and re-index X, and dX, to match size of re-indexed TI, XI
            x=x(bsxfun(@plus,ii,0:M:M*(NN-1)));
            dx=dx(bsxfun(@plus,ii,0:(M-1):(M-1)*(NN-1)));

            % Calculate XI, use sort indices for TI to reduce dimensions
            if isscalar(eta)
                xi(p,:)=x+bsxfun(@times,dx,mu)+bsxfun(@times,randn(lp,NN,dtype),eta.*sig);
            else
                xi(p,:)=x+bsxfun(@times,dx,mu)+sig*eta'.*randn(lp,NN,dtype);
            end
        end
    end
end

% Reshape output to produce N-D array or to match TI orientation
if ~isVecX
    if linspc || (length(szti) == 2 && any(szti == 1))
        szti=[D N];
    else
        szti=[szti N];
    end
end
if length(szti) ~= 2 || any(szti ~= [D NN])
    xi=reshape(xi,szti);
end