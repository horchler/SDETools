function sde_install(opt)
%SDE_INSTALL  Add or remove SDETools from Matlab search path and save path.
%   SDE_INSTALL adds the SDETools directory (where this function is located) to
%   the Matlab search path, saves the path, and prints the help for the toolbox.
%
%   SDE_INSTALL('remove') uninstalls SDETools by removing the SDETools directory
%   from the Matlab search path and saving the path. If SDETools is not
%   installed (on the path), a warning is issued.
%
%   See also: PATH, ADDPATH, RMPATH, SAVEPATH

%   Andrew D. Horchler, horchler @ gmail . com, Created 8-12-13
%   Revision: 1.2, 11-16-13


if nargin == 0 || any(strcmp(opt,{'add','install','addpath'}))
    addpath(fileparts(mfilename('fullpath')));
    status = true;
elseif any(strcmp(opt,{'remove','uninstall','rmpath'}))
    rmpath(fileparts(mfilename('fullpath')));
    status = false;
else
    error('SDETools:sde_install:UnknownOption',...
         ['Input argument must be the string ''install'' to install or '...
          '''remove'' to uninstall.']);
end

if savepath
    error('SDETools:sde_install:SavePathError',...
          'Unable to save pathdef.m file.');
end
rehash('toolbox');
clear('sde_install');

if status
    fprintf(1,'\n SDETools installed.\n\n');
    help('SDETools');
else
    fprintf(1,'\n SDETools uninstalled.\n');
end