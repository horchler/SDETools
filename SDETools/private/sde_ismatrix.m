function bool=sde_ismatrix(V)
%SDE_ISMATRIX  Replicate functionality of builtin ismatrix for pre-R2010b Matlab
%
%   See also: ISMATRIX, SDE_ISSQUARE, SDE_ISDIAG

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-9-12
%   Revision: 1.2, 5-4-13


try
    bool = ismatrix(V);
catch err
    if strcmp(err.identifier,'MATLAB:UndefinedFunction')
        bool = ndims(V) == 2 || size(V(:,:,:),3) == 1;	%#ok<ISMAT>
    else
        rethrow(err);
    end
end