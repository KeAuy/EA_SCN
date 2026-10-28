function grid = geometricGrid(rMin, rMax, nGrid)
%GEOMETRICGRID Geometric grid and local spacings used by the figure panels.
%
%   GRID = BUFFER.GEOMETRICGRID(RMIN, RMAX, NGRID) returns a struct with
%   fields R, HLEFT and HRIGHT for the grid
%
%       r(k) = 10^(log10(rMin) + (log10(rMax) - log10(rMin)) k / NGRID),
%
%   which satisfies r(k+1)/r(k) = constant.  HLEFT(k) = r(k) - r(k-1) and
%   HRIGHT(k) = r(k+1) - r(k).
%
%   See also BUFFER.SOLVEGEOMETRIC.

arguments
    rMin (1,1) double {mustBePositive}
    rMax (1,1) double {mustBePositive}
    nGrid (1,1) double {mustBeInteger, mustBeGreaterThan(nGrid, 2)}
end

if rMax <= rMin
    error('buffer:invalidGrid', 'rMax must be larger than rMin.');
end

k     = (1:nGrid)';
r     = 10.^(log10(rMin) + (log10(rMax) - log10(rMin)) * k / nGrid);
alpha = r(2) / r(1);

grid = struct( ...
    'r',      r, ...
    'hLeft',  r * (1 - 1 / alpha), ...
    'hRight', r * (alpha - 1));
end
