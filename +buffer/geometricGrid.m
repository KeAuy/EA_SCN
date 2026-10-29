%GEOMETRICGRID Geometric radial grid r(k+1)/r(k) = const with local spacings.
%
%   GRID = BUFFER.GEOMETRICGRID(RMIN, RMAX, NGRID) returns R, HLEFT and
%   HRIGHT, where HLEFT(k) = r(k) - r(k-1) and HRIGHT(k) = r(k+1) - r(k).

function grid = geometricGrid(rMin, rMax, nGrid)

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
