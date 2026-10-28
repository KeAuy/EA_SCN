function colors = palette()
%PALETTE Shared line colors for the figure panels.
%
%   COLORS = BUFFER.PALETTE returns a struct with the colors used by the
%   original figure scripts.
%
%   See also BUFFER.PLOTPROFILES, BUFFER.PLOTERRORS.

colors = struct();
colors.numerical = [0.75 0.75 0.75];
colors.pade      = [46 139 87] / 256;
colors.blue      = [0 54 253] / 256;
colors.red       = 'r';
colors.green     = 'g';
colors.black     = 'k';
colors.magenta   = 'm';
end
