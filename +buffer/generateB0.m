%GENERATEB0 Inner boundary value b0(nu) computed with BVP5C.
%
%   B0 = BUFFER.GENERATEB0(LAMBDA, ETA, NUVALUES, OPTIONS) solves the buffer
%   BVP on the geometric grid for every entry of NUVALUES, using the previous
%   solution as the initial guess for the next nu (continuation).
%
%   OPTIONS: GridSize, Tolerance, NMax, InitialB0, Verbose.  The result is
%   returned in memory only: no data file is written, so the repository stays
%   code-only.

function b0 = generateB0(lambda, eta, nuValues, options)

arguments
    lambda (1,1) double {mustBePositive}
    eta (1,1) double
    nuValues double
    options.GridSize = []
    options.Tolerance (1,1) double {mustBePositive} = 1e-8
    options.NMax (1,1) double {mustBeInteger, mustBePositive} = 50000
    options.InitialB0 (1,1) double = 0.5
    options.Verbose (1,1) logical = true
end

nuValues = nuValues(:);
if isempty(options.GridSize)
    options.GridSize = numel(nuValues);
end
if options.GridSize < 3
    error('buffer:invalidGrid', 'GridSize must be at least three.');
end

grid = buffer.geometricGrid(1e-3, 1e2, options.GridSize);
r    = grid.r;
nNu  = numel(nuValues);
b0   = zeros(1, nNu);

bvpOptions = bvpset('NMax', options.NMax, ...
    'RelTol', options.Tolerance, 'AbsTol', options.Tolerance);

nu = nuValues(1);
sol = bvp5c( ...
    @(rr, y) buffer.ode(rr, y, lambda, nu, eta), ...
    @(ya, yb) buffer.bc(ya, yb, lambda), ...
    bvpinit(r, @(rr) buffer.guess(rr, options.InitialB0, lambda, nu, eta)), ...
    bvpOptions);

y      = deval(sol, r);
b0(1)  = y(1, 1);
if options.Verbose
    fprintf('nu = %.3e completed\n', nu);
end

for k = 2:nNu
    nu  = nuValues(k);
    sol = bvp5c( ...
        @(rr, y) buffer.ode(rr, y, lambda, nu, eta), ...
        @(ya, yb) buffer.bc(ya, yb, lambda), ...
        sol, bvpOptions);

    y       = deval(sol, r);
    b0(k)   = y(1, 1);
    if options.Verbose
        fprintf('nu = %.3e completed\n', nu);
    end
end
end
