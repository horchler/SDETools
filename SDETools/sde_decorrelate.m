function varargout=sde_decorrelate(d,r2)
%SDE_DECORRELATE  Decorrelate correlated values.
%   R1 = SDE_DECORRELATE(D, R2) returns the matrix R1 of M decorrelated values
%   given the N-by-N diffusion matrix D and the M-by-N matrix of (correlated)
%   values R2. The first column of R1 is always SQRT(D(1,1)) divided by the
%   first column of R2.
%
%   [R1, C] = SDE_DECORRELATE(D, R2) also returns the N-by-N correlation matrix,
%   C, decorrelated from the diffusion matrix, D.
%
%   C = SDE_DECORRELATE(D) without a second input argument returns the N-by-N
%   correlation matrix, C, decorrelated from the diffusion matrix, D.
%
%   Note:
%       R1 = SDE_DECORRELATE(D, SDE_CORRELATE(C, R1)) is only accurate to within
%       the floating point machine precision, EPS.
%
%   Example:
%       % Correlate and then decorrelate normally-distributed samples
%       r1 = randn(5,3)
%       c = [1 0.8 0.2;0.8 1 0.5;0.2 0.5 1];
%       [r2, d] = sde_correlate(c,r1)
%       r3 = sde_decorrelate(d,r2)
%
%   See also: SDE_CORRELATE, RAND, RANDN, RANDSTREAM, RANDSTREAM/RANDN, EPS

%   Andrew D. Horchler, adh9 @ case . edu, Created 5-20-13
%   Revision: 1.2, 7-12-13


if nargout > min(nargin,2)
    error('SDETools:sde_decorrelate:TooManyOutputs',...
          'Too many output arguments for number of supplied inputs.');
end

if isempty(d) || isempty(r2)
    if nargin == 1
        varargout{1} = [];
    else
        varargout{2} = [];
    end
else
    [m,n] = size(d);
    if ndims(d) ~= 2 || m ~= n              %#ok<ISMAT>
        error('SDETools:sde_decorrelate:NonSquareMatrix',...
              'The diffusion matrix must be square.');
    end
    
    if isscalar(d)
        c = d^2;
        isDiag = false;
	elseif sde_isdiag(d)
        d = diag(d);
        c = diag(d.^2);
        isDiag = true;
    else
        c = d^2;
        isDiag = false;
    end
    
    if nargin == 1
        varargout{1} = c;
    else
        if ndims(r2) ~= 2 || size(r2,2) ~= m	%#ok<ISMAT>
            error('SDETools:sde_decorrelate:DimensionMismatch',...
                 ['The number of columns in the matrix of normally '...
                  'distributed values must equal the dimension of the '...
                  'correlation matrix.']);
        end
        
        if isDiag
            varargout{1} = bsxfun(@mrdivide,r2,d);
        else
            varargout{1} = r2/d;
        end
        if nargout == 2
            varargout{2} = c;
        end
    end
end