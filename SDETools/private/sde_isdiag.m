function bool=sde_isdiag(V,k)
%SDE_ISDIAG  True for diagonal matrices.
%   ISDIAG(V) returns logical 1 (true) if V is a diagonal matrix and logical 0
%   (false) otherwise. A diagonal matrix is a non-empty square matrix with
%   arbitrary entries along the main diagonal and is zero elsewhere.
%
%   ISDIAG(V,K) returns logical 1 (true) if V is diagonal along the K-th
%   diagonal and logical 0 (false) otherwise. K > 0 and K < 0 are above and
%   below the main diagonal, respectively. ISDIAG(V,0) is the same as ISDIAG(V).
%
%   See also: SDE_ISMATRIX, SDE_ISSQUARE, DIAG

%   http://mathworks.com/matlabcentral/newsreader/view_thread/299383#806642

%   Andrew D. Horchler, horchler @ gmail . com, Created 4-30-13
%   Revision: 1.2, 7-12-13


bool = false;
if ~isempty(V)
    if nargin == 1
        if isscalar(V)
            bool = true;
        else
            [m,n] = size(V);
            if ndims(V) == 2 && m == n	%#ok<ISMAT>
                bool = (nnz(V) == nnz(diag(V)));
            end
        end
    else
        if ~isscalar(k) || ~isnumeric(k) || any(k ~= floor(k))
            error('SDETools:sde_isdiag:kthDiagInputNotInteger',...
                  'K-th diagonal input must be a scalar integer.')
        end
        
        if isscalar(V)
            bool = (k == 0);
        else
            [m,n] = size(V);
            if ndims(V) == 2 && m == n && abs(k) < n	%#ok<ISMAT>
                bool = (nnz(V) == nnz(diag(V,k)));
            end
        end
    end
end