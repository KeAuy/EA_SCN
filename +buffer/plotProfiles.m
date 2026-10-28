function ax = plotProfiles(x, values, labels, options)
%PLOTPROFILES Plot one or more b(r) profiles on a logarithmic r axis.
%
%   AX = BUFFER.PLOTPROFILES(X, VALUES, LABELS, OPTIONS) creates a figure
%   and plots REAL(VALUES{K}) against X for every entry of the cell array
%   VALUES.  OPTIONS controls line styles, colors, widths, axis limits,
%   labels, title, subtitle and legend location.
%
%   See also BUFFER.PLOTERRORS, BUFFER.PALETTE.

arguments
    x double {mustBePositive}
    values cell
    labels cell
    options.Styles cell = {}
    options.Colors cell = {}
    options.LineWidths double = []
    options.Title = ''
    options.Subtitle (1,:) char = ''
    options.XLabel (1,:) char = 'r (distance)'
    options.YLabel (1,:) char = 'b(r)'
    options.XLim (1,2) double = [1e-3 1e2]
    options.XTicks double = [1e-3 1e-2 1e-1 1 10 100]
    options.XTickLabels cell = {}
    options.YLim double = []
    options.YTicks double = []
    options.YTickLabels cell = {}
    options.LegendLocation (1,:) char = 'southeast'
    options.Box (1,:) char {mustBeMember(options.Box, {'on', 'off'})} = 'off'
end

x = x(:);
n = numel(values);
if n ~= numel(labels)
    error('buffer:sizeMismatch', ...
        'VALUES and LABELS must have the same number of entries.');
end
for k = 1:n
    if numel(values{k}) ~= numel(x)
        error('buffer:sizeMismatch', ...
            'Every profile must have the same length as X.');
    end
end

styles = options.Styles;
if isempty(styles), styles = repmat({'-'}, 1, n); end
colors = options.Colors;
if isempty(colors), colors = repmat({'k'}, 1, n); end
widths = options.LineWidths;
if isempty(widths), widths = repmat(2, 1, n); end
if isscalar(widths), widths = repmat(widths, 1, n); end

figure;
ax = gca;
hold(ax, 'on');
box(ax, options.Box);

for k = 1:n
    plot(x, real(values{k}), ...
        'LineStyle', styles{k}, ...
        'Color', colors{k}, ...
        'LineWidth', widths(k));
end

set(ax, 'XScale', 'log', 'XLim', options.XLim, 'XTick', options.XTicks);
if ~isempty(options.XTickLabels)
    set(ax, 'XTickLabel', options.XTickLabels);
end

xlabel(options.XLabel);
ylabel(options.YLabel);
if ~isempty(options.Title)
    title(options.Title);
end
if ~isempty(options.Subtitle)
    subtitle(options.Subtitle);
end
if ~isempty(options.YLim)
    ylim(options.YLim);
end
if ~isempty(options.YTicks)
    yticks(options.YTicks);
end
if ~isempty(options.YTickLabels)
    yticklabels(options.YTickLabels);
end

legend(labels, 'Location', options.LegendLocation);
end
