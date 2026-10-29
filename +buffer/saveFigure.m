%SAVEFIGURE Export a figure to the shared figures folder as a PNG.

function saveFigure(handle, name)
arguments
    handle
    name (1,:) char
end

exportgraphics(handle, fullfile(buffer.outputDir(), name), 'Resolution', 160);
end
