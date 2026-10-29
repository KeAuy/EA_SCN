%OUTPUTDIR Folder that collects every reproduced figure (created on demand).

function folder = outputDir()
rootDir = fileparts(fileparts(mfilename('fullpath')));
folder  = fullfile(rootDir, 'figures');
if ~isfolder(folder)
    mkdir(folder);
end
end
