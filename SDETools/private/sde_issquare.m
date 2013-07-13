function bool=sde_issquare(V)
%SDE_ISSQUARE  True for square matrices.
%
%   See also: SDE_ISMATRIX, SDE_ISDIAG

%   Andrew D. Horchler, adh9 @ case . edu, Created 4-30-13
%   Revision: 1.2, 5-4-13


bool = (ndims(V) == 2 && size(V,1) == size(V,2));	%#ok<ISMAT>