%PLOTWINNERMAP Method-region map on a lambda-nu grid.
%
%   Draws WINNER with one colour per method, adds an optional set of zero
%   contours (the alpha real/complex boundaries) and exports the panel.

function ax = plotWinnerMap(nuValues, lambdaValues, winner, methodNames, colors, options)

arguments
    nuValues double
    lambdaValues double
    winner double
    methodNames cell
    colors double
    options.Title (1,:) char = ''
    options.Contours struct = struct([])
    options.OutputFile (1,:) char = ''
end

nuValues     = nuValues(:);
lambdaValues = lambdaValues(:);
nMethod      = numel(methodNames);

if ~isequal(size(winner), [numel(lambdaValues), numel(nuValues)])
    error('buffer:sizeMismatch', ...
        'WINNER must be NLAMBDA-by-NNU.');
end
if ~isequal(size(colors), [nMethod, 3])
    error('buffer:sizeMismatch', ...
        'COLORS must be NMETHOD-by-3.');
end

nanColor = [0.97 0.97 0.97];

figure('Color', 'w', 'Position', [100 100 720 560]);
ax = axes();
imagesc(ax, nuValues, lambdaValues, winner);
set(ax, 'XScale', 'log', 'YScale', 'log', ...
    'XLim', [1e-3 1e2], 'YLim', [1e-3 1e2], ...
    'XTick', 10.^(-3:2), 'YTick', 10.^(-3:2), ...
    'XTickLabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '10^1', '10^2'}, ...
    'YTickLabel', {'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '10^1', '10^2'});
axis(ax, 'xy');
clim(ax, [0.5, nMethod + 1.5]);
colormap(ax, [colors; nanColor]);

cb = colorbar(ax);
cb.Ticks      = 1:nMethod + 1;
cb.TickLabels = [methodNames, {'none'}];
cb.FontSize   = 9;

xlabel(ax, '\nu');
ylabel(ax, '\lambda');
if ~isempty(options.Title)
    title(ax, options.Title);
end
box(ax, 'on');

axes(ax);
hold(ax, 'on');

for k = 1:numel(options.Contours)
    contourData = options.Contours(k).Z;
    if any(contourData(:) < 0) && any(contourData(:) > 0)
        contour(nuValues, lambdaValues, contourData, [0 0], ...
            'LineColor', options.Contours(k).Color, ...
            'LineWidth', options.Contours(k).LineWidth);
    end
end

if ~isempty(options.OutputFile)
    exportgraphics(ax.Parent, options.OutputFile, 'Resolution', 160);
end
end
