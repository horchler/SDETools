function varargout=sde_correlate(c,r1)
%SDE_CORRELATE  Correlated values.
%   R2 = SDE_CORRELATE(C, R1) returns the matrix R2 of M correlated values given
%   the N-by-N correlation matrix C and the M-by-N matrix of (uncorrelated)
%   values R1. The first column of R2 is always SQRT(C(1,1)) times the first
%   column of R1.
%
%   [R2, D] = SDE_CORRELATE(C, R1) also returns the N-by-N diffusion matrix, D,
%   created from the correlation matrix, C.
%
%   D = SDE_CORRELATE(C) without a second input argument returns only the N-by-N
%   diffusion matrix, D, created from the correlation matrix, C.
%
%   Note:
%       C may either be a symmetric positive semidefinite correlation matrix or
%       covariance matrix.
%
%   Example:
%       % Plot uncorrelated and correlated normally-distributed points
%       r1 = randn(1e4,2); corr(r1)
%       r2 = sde_correlate([1 -0.8;-0.8 1],r1); corr(r2)
%       figure; plot(r1(:,1),r1(:,2),'b.',r2(:,1),r2(:,2),'r*'); axis equal;
%
%   See also: SDE_DECORRELATE, RAND, RANDN, RANDSTREAM, RANDSTREAM/RANDN

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-20-13
%   Revision: 1.2, 7-17-13


if nargout > min(nargin,2)
    error('SDETools:sde_correlate:TooManyOutputs',...
          'Too many output arguments for number of supplied inputs.');
end

if isempty(c) || isempty(r1)
    if nargin == 1
        varargout{1} = [];
    else
        varargout{2} = [];
    end
else
    [m,n] = size(c);
    if ndims(c) ~= 2 || m ~= n              %#ok<ISMAT>
        error('SDETools:sde_correlate:NonSquareMatrix',...
              'The correlation matrix must be square.');
    end
    
    if sde_isdiag(c)
        num = rank(c);
        c = diag(c);
        if any(c) < 0
            error('SDETools:sde_correlate:NegativeDiagonal',...
                 ['All elements of the diagonal of the correlation matrix '...
                  'must be non-negative.']);
        end
        
        d = sqrt(c);
        isDiag = isscalar(c);
    else
        [d,num] = cholcov(c);
        if isempty(d)
            error('SDETools:sde_correlate:NonSymmetricSemiDefinite',...
                 ['The correlation matrix must be a symmetric positive '...
                  'semidefinite matrix.']);
        end
        isDiag = false;
    end
    
    if nargin == 1
        if num ~= m
            warning('SDETools:sde_correlate:NonFullRank1',...
                    'The correlation matrix is not full rank.');
        end
        
        if isDiag
            varargout{1} = diag(d);
        else
            varargout{1} = d;
        end
    else
        if ndims(r1) ~= 2 || size(r1,2) ~= m	%#ok<ISMAT>
            error('SDETools:sde_correlate:DimensionMismatch',...
                 ['The number of columns in the matrix of normally '...
                  'distributed values must equal the dimension of the '...
                  'correlation matrix.']);
        end
        if num ~= m
            warning('SDETools:sde_correlate:NonFullRank2',...
                   ['The correlation matrix is not full rank. Only the '...
                    'first %d columns of the input values will be '...
                    'correlated.'],num);
        end
        
        if isDiag
            varargout{1} = bsxfun(@mtimes,r1,d);
            if nargout == 2
                varargout{2} = diag(d);
            end
        else
            varargout{1} = r1*d;
            if nargout == 2
                varargout{2} = d;
            end
        end
    end
end