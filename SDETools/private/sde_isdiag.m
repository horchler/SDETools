function bool=sde_isdiag(V,k)
%SDE_ISDIAG  True for diagonal matrices.
%
%   See also: SDE_ISMATRIX, SDE_ISSQUARE

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-30-13
%   Revision: 1.0, 4-30-13


if ~isempty(V) && ndims(V) == 2	%#ok<ISMAT>
    [m,n] = size(V);
    if m == n
        if nargin == 2
            if ~isscalar(k) || k ~= floor(k)
                error('SDETools:sde_isdiag:kthDiagInputNotInteger',...
                      'K-th diagonal input must be an integer scalar.')
            end
            if abs(k) < n
                if k < 0
                    bool = all(V(mod(1-k*n:end,n+1) ~= 1) == 0);
                else
                    bool = all(V(mod(1+k:end,n+1) ~= 1) == 0);
                end
                if issparse(V)
                    bool = full(bool);
                end
            else
                bool = false;
            end
        else
            bool = all(V(mod(1:end,n+1) ~= 1) == 0);
            if issparse(V)
                bool = full(bool);
            end
        end
    else
        bool = false;
    end
else
    bool = false;
end