function ax = plotErrors(x, errors, labels, options)
%PLOTERRORS Plot one or more error curves against nu on a logarithmic x axis.
%
%   AX = BUFFER.PLOTERRORS(X, ERRORS, LABELS, OPTIONS) creates a figure and
%   plots every entry of the cell array ERRORS against X.  OPTIONS.RANGES
%   can restrict individual curves to a subset of indices (for example the
%   first 60 nu values of the Exp-Pade curve).  OPTIONS.MARKMINIMUM marks the
%   minimum of the selected curve.
%
%   See also BUFFER.PLOTPROFILES, BUFFER.MEANERROR.

arguments
    x double {mustBePositive}
    errors cell
    labels cell
    options.Styles cell = {}
    options.Colors cell = {}
    options.LineWidths double = []
    options.Ranges cell = {}
    options.MarkMinimum double = []
    options.Title = ''
    options.Subtitle (1,:) char = ''
    options.XLabel (1,:) char = '\nu'
    options.YLabel (1,:) char = 'log_{10}||b_{num} - b_{approx}||'
    options.XLim (1,2) double = [1e-3 1e2]
    options.XTicks double = [1e-3 1e-2 1e-1 1 10 100]
    options.XTickLabels cell = ...
        {'10^{-3}', '10^{-2}', '10^{-1}', '10^0', '10^1', '10^2'}
    options.YLim (1,2) double = [-5.5 -0.5]
    options.YTicks double = [-5.5 -5 -4 -3 -2 -1 -0.5]
    options.YTickLabels cell = {' ', '-5', '-4', '-3', '-2', '-1', ' '}
    options.LegendLocation (1,:) char = 'southeast'
    options.Box (1,:) char {mustBeMember(options.Box, {'on', 'off'})} = 'on'
end

x = x(:);
n = numel(errors);
if n ~= numel(labels)
    error('buffer:sizeMismatch', ...
        'ERRORS and LABELS must have the same number of entries.');
end
for k = 1:n
    if numel(errors{k}) ~= numel(x)
        error('buffer:sizeMismatch', ...
            'Every error vector must have the same length as X.');
    end
end

styles = options.Styles;
if isempty(styles), styles = repmat({'-'}, 1, n); end
colors = options.Colors;
if isempty(colors), colors = repmat({'k'}, 1, n); end
widths = options.LineWidths;
if isempty(widths), widths = repmat(2, 1, n); end
if isscalar(widths), widths = repmat(widths, 1, n); end

ranges = options.Ranges;
if ~isempty(ranges) && numel(ranges) ~= n
    error('buffer:sizeMismatch', ...
        'RANGES must be empty or have one entry per error curve.');
end

figure;
ax = gca;
hold(ax, 'on');
box(ax, options.Box);

for k = 1:n
    if isempty(ranges) || isempty(ranges{k})
        idx = 1:numel(x);
    else
        idx = ranges{k};
    end

    lineStyle = styles{k};
    marker    = '';
    if numel(lineStyle) > 1 && lineStyle(end) == '.' ...
            && any(strcmp(lineStyle(1:end-1), {'-', '--', ':', '-.'}))
        marker    = '.';
        lineStyle = lineStyle(1:end-1);
    end

    if isempty(marker)
        plot(x(idx), errors{k}(idx), ...
            'LineStyle', lineStyle, ...
            'Color', colors{k}, ...
            'LineWidth', widths(k));
    else
        plot(x(idx), errors{k}(idx), ...
            'LineStyle', lineStyle, ...
            'Marker', marker, ...
            'Color', colors{k}, ...
            'LineWidth', widths(k));
    end
end

if ~isempty(options.MarkMinimum)
    k = options.MarkMinimum;
    if k < 1 || k > n
        error('buffer:invalidIndex', ...
            'MarkMinimum must be a valid index into ERRORS.');
    end
    [minValue, minIndex] = min(errors{k});
    plot(x(minIndex), minValue, 'o', ...
        'MarkerSize', 8, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colors{k});
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

ylim(options.YLim);
yticks(options.YTicks);
if ~isempty(options.YTickLabels)
    yticklabels(options.YTickLabels);
end

legend(labels, 'Location', options.LegendLocation);
end
