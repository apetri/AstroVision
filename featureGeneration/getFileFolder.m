%% Returns: a cell array of all the file names in a folder
%
function files = getFileFolder(folder)
 
files = dir(folder);
% listing also includes two special dirs, './' and '../'
files = {files(~[files.isdir]).name};
