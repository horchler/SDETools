function install(opt)
%INSTALL  Add or remove SDETools from Matlab search path and save path.
%   INSTALL adds the SDETools directory (where this function is located) to the
%   Matlab search path, saves the path, and prints the help for the toolbox.
%
%   INSTALL('remove') uninstalls SDETools by removing the SDETools directory
%   from the Matlab search path and saving the path. If SDETools is not
%   installed (on the path), a warning is issued.
%
%   See also path, addpath, rmpath, savepath

%   Andrew D. Horchler, adh9 @ case . edu, Created 8-12-13
%   Revision: 1.2, 8-12-13


[success,msg] = fileattrib;
if success
    if nargin == 0 || any(strcmp(opt,{'add','install','addpath'}))
        addpath(msg.Name);
        status = true;
    elseif any(strcmp(opt,{'remove','uninstall','rmpath','delete'}))
        rmpath(msg.Name);
        status = false;
    else
        error('SDETools:install:UnknownOption',...
              'Input argument must be the string ''remove'' to uninstall.');
    end
else
    error('SDETools:install:AddPathError',...
          'Unable to find absolute path of this directory.');
end

if savepath
    error('SDETools:install:SavePathError','Unable to save pathdef.m file.');
end
rehash('toolbox');

if status
    fprintf(1,'\n SDETools installed.\n\n');
    help('SDETools');
else
    fprintf(1,'\n SDETools uninstalled.\n');
end