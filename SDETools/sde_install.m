function sde_install(opt)
%SDE_INSTALL  Add or remove SDETools from Matlab search path and save path.
%   SDE_INSTALL adds the SDETools directory (where this function is located) to
%   the Matlab search path, saves the path, and prints the help for the toolbox.
%
%   SDE_INSTALL('remove') uninstalls SDETools by removing the SDETools directory
%   from the Matlab search path and saving the path. If SDETools is not
%   installed (on the path), a warning is issued.
%
%   See also path, addpath, rmpath, savepath

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-12-13
%   Revision: 1.2, 8-30-13


[success,msg] = fileattrib;
if success
    if nargin == 0 || any(strcmp(opt,{'add','install','addpath'}))
        addpath(msg.Name);
        status = true;
    elseif any(strcmp(opt,{'remove','uninstall','rmpath','delete'}))
        rmpath(msg.Name);
        status = false;
    else
        error('SDETools:sde_install:UnknownOption',...
              'Input argument must be the string ''remove'' to uninstall.');
    end
else
    error('SDETools:sde_install:AddPathError',...
          'Unable to find absolute path of this directory.');
end

if savepath
    error('SDETools:sde_install:SavePathError',...
          'Unable to save pathdef.m file.');
end
rehash('toolbox');

if status
    fprintf(1,'\n SDETools installed.\n\n');
    help('SDETools');
else
    fprintf(1,'\n SDETools uninstalled.\n');
end