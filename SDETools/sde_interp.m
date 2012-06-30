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
%   Class support for inputs T, X, TI, NT:
%       float: double, single
%
%   See also SDE_INTERPQ, SDE_INTERPNQ, INTERP1, INTERP1Q

% Some code based on version 5.41.4.18 of Matlab's INTERP1

%   Andrew D. Horchler, adh9 @ case . edu, Created 2-26-12
%   Revision: 1.0, 6-30-12


%XI = SDE_INTERP(T,X,TI) nargin=3, offset=1
%XI = SDE_INTERP(T,X,'linear',NT) nargin=4, offset=1
%XI = SDE_INTERP(T,X,TI,'extrap',EXTRAPVAL) nargin=5, offset=1

%XI = SDE_INTERP(X,TI) nargin=2, offset=0
%XI = SDE_INTERP(X,'linear',NT) nargin=3, offset=0
%XI = SDE_INTERP(X,TI,'extrap',EXTRAPVAL) nargin=4, offset=0

if nargin < 2
    error('SDETools:sde_interp:TooFewInputs','Not enough input arguments.');
end
if nargin > 5
    error('SDETools:sde_interp:TooManyInputs','Too many input arguments.');
end

% Handle linearly-spaced parameter and optional extrapolation value
extrap=false;
linspc=false;
offset=1;
if nargin == 5
    str=varargin{4};
    if ~ischar(str) || ~strcmp(str,'extrap')
        error('SDETools:sde_interp:InvalidExtrapString',...
             ['The parameter string ''extrap'' must be used to specify an '...
              'optional EXTRAPVAL in the subsequent input argument.']);
    end
    extrapval=varargin{5};
    extrap=true;
elseif nargin == 4
    str=varargin{3};
    if ~ischar(str) || ~any(strcmp(str,{'linear','extrap'}))
        error('SDETools:sde_interp:InvalidParameterString',...
             ['The parameter strings ''linear'' or ''extrap'' must be used '...
              'to specify the the number of linearly-spaced interpolation '...
              'points, NT, or an optional EXTRAPVAL in the subsequent input '...
              'argument.']);
    end
    if strcmp(str,'extrap')
        extrapval=varargin{4};
        extrap=true;
        offset=0;
    else
        nt=varargin{4};
        linspc=true;
    end
elseif nargin == 3 && ischar(varargin{2})
    if ~strcmp(varargin{2},'linear')
        error('SDETools:sde_interp:InvalidLinearString',...
             ['The parameter string ''linear'' must be used to specify the '...
              'the number of linearly-spaced interpolation points, NT, in '...
              'the subsequent input argument.']);
    end
    nt=varargin{3};
    linspc=true;
    offset=0;
elseif nargin == 2
    offset=0;
end

if extrap
    if numel(extrapval) ~= 1
        error('SDETools:sde_interp:NonScalarExtrapVal',...
              'The input EXTRAPVAL must be a scalar value.');
    end
elseif linspc
    if ~iscalar(nt) || ~isnumeric(nt) || ~isreal(nt) || ~isempty(nt) ...
            && (~isfinite(nt) || (~isinteger(nt) && nt-floor(nt) ~= 0) ...
            || nt < 0)
        error('SDETools:sde_interp:InvalidNT',...
              'The input NT must be finite real positive scalar value.');
    end
end

% Check X
x=varargin{offset+1};
if isempty(x) || ~isfloat(x)
    error('SDETools:sde_interp:XEmptyOrNotFloat',...
          'The input X must be a non-empty array of singles or doubles.');
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
if offset == 0
    t=[0;M-1];
    dt=1;
    eqsp=true;
else
    t=varargin{offset};
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
        [t p]=sort(t);
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

% Check TI or construct linearly-spaced TI from NT
if linspc
    dtype=superiorfloat(t,x,nt);
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
        
        if isVecX
            % Calculate XI, output is a D-by-1 column vector
            xi=x+mu.*diff(x)+sig.*randn(D,1,dtype);
        else
            % Expand and re-index X, and dX, to match size of re-indexed XI
            dx=diff(x);
            x=x(bsxfun(@plus,ii,0:M:M*(NN-1)));
            dx=dx(bsxfun(@plus,ii,0:(M-1):(M-1)*(NN-1)));

            % Calculate XI, output is a D-by-N matrix
            xi=x+bsxfun(@times,dx,mu)+bsxfun(@times,randn(D,NN,dtype),sig);
        end
    end
    szti=[D 1];
else
    ti=varargin{offset+2};

    if ~isfloat(ti) || ~isreal(ti)
        error('SDETools:sde_interp:InvalidTiDataType',...
              'Datatype of input Ti must be real single or real double.');
    end
    szti=size(ti);
    ti=ti(:);
    D=length(ti);
    
    % Allocate out of bounds values as NaN, or override default extrapolation
    dtype=superiorfloat(t,x,ti);
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
        sig=0*sqrt((1-mu).*dti);
        
        dx=diff(x);
        if isVecX
            % Calculate XI, use sort indices for TI to reduce dimensions
            xi(p)=x(ii)+mu.*dx(ii)+sig.*randn(lp,1,dtype);
        else
            % Expand and re-index X, and dX, to match size of re-indexed TI, XI
            x=x(bsxfun(@plus,ii,0:M:M*(NN-1)));
            dx=dx(bsxfun(@plus,ii,0:(M-1):(M-1)*(NN-1)));

            % Calculate XI, use sort indices for TI to reduce dimensions
            xi(p,:)=x+bsxfun(@times,dx,mu)+bsxfun(@times,randn(lp,NN,dtype),...
                sig);
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