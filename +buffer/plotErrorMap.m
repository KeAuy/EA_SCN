%PLOTERRORMAP Filled contour map of an error field on a lambda-nu grid.
%
%   OPTIONS: Levels, CLim, ColorTicks, Title, XLabel, YLabel, ColorMap and
%   OutputFile.  Both axes are logarithmic.

function ax = plotErrorMap(nuValues, lambdaValues, errorMap, options)

arguments
    nuValues double
    lambdaValues double
    errorMap double
    options.Levels double = []
    options.CLim double = []
    options.ColorTicks double = []
    options.Title (1,:) char = ''
    options.XLabel (1,:) char = '\nu'
    options.YLabel (1,:) char = '\lambda'
    options.ColorMap (1,:) char = 'bone'
    options.OutputFile (1,:) char = ''
end

nuValues     = nuValues(:);
lambdaValues = lambdaValues(:);

if ~isequal(size(errorMap), [numel(lambdaValues), numel(nuValues)])
    error('buffer:sizeMismatch', ...
        'ERRORMAP must be NLAMBDA-by-NNU.');
end

figure('Color', 'w');
ax = axes();
if isempty(options.Levels)
    contourf(ax, nuValues, lambdaValues, errorMap);
else
    contourf(ax, nuValues, lambdaValues, errorMap, options.Levels);
end

set(ax, 'XScale', 'log', 'YScale', 'log', ...
    'XLim', [1e-3 1e2], 'YLim', [1e-3 1e2], ...
    'XTick', 10.^(-3:2), 'YTick', 10.^(-3:2), ...
    'XTickLabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '10^1', '10^2'}, ...
    'YTickLabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '10^1', '10^2'});

colormap(ax, options.ColorMap);
if isempty(options.CLim)
    clim(ax, 'auto');
else
    clim(ax, options.CLim);
end

cb = colorbar(ax);
if ~isempty(options.ColorTicks)
    cb.Ticks = options.ColorTicks;
end

xlabel(ax, options.XLabel);
ylabel(ax, options.YLabel);
if ~isempty(options.Title)
    title(ax, options.Title);
end

if ~isempty(options.OutputFile)
    exportgraphics(ax.Parent, options.OutputFile, 'Resolution', 160);
end
end
