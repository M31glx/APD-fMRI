function [status,output] = call_fsl_edit(cmd)
% [status, output] = call_fsl(cmd)
% 
% Wrapper around calls to FSL binaries
% clears LD_LIBRARY_PATH and ensures
% the FSL envrionment variables have been
% set up
% Debian/Ubuntu users should uncomment as
% indicated

if strcmp(cmd(1),'!')
    cmd = cmd(2:end);
end

if isempty(getenv('FSLDIR'))
    setenv('FSLDIR','/usr/local/fsl')
end
path_default=getenv('PATH');
path =[path_default,':',getenv('FSLDIR'),'/bin'];
setenv('PATH',path)
setenv('FSLOUTPUTTYPE','NIFTI_GZ')

fsldir=getenv('FSLDIR');

% Debian/Ubuntu - uncomment the following
%fsllibdir=sprintf('%s/%s', fsldir, 'bin');

if ismac
  dylibpath=getenv('DYLD_LIBRARY_PATH');
  setenv('DYLD_LIBRARY_PATH');
else
  ldlibpath=getenv('LD_LIBRARY_PATH');
  setenv('LD_LIBRARY_PATH');
  % Debian/Ubuntu - uncomment the following
  % setenv('LD_LIBRARY_PATH',fsllibdir);
end

command = sprintf('/bin/sh -c ''. %s/etc/fslconf/fsl.sh; %s''', fsldir, cmd);
[status,output] = system(command);

if ismac
  setenv('DYLD_LIBRARY_PATH', dylibpath);
else
    setenv('LD_LIBRARY_PATH', ldlibpath);
end

if status
    error('FSL call (%s) failed, %s', command, output)
end