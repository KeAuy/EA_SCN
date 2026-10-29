%GENERATEB0ROW One row of the b0(lambda, nu) grid.
%
%   B0 = BUFFER.GENERATEB0ROW(LAMBDA, NUVALUES, ETA, R, TOLERANCE, NMAX,
%   INITIALB0) solves the buffer BVP for a single lambda, continuing along the
%   entries of NUVALUES.  R is the radial grid shared by all rows.
%
%   This is the body of the PARFOR loop in BUFFER.GENERATEB0GRID, kept in its
%   own file so that every worker can call it directly.

function b0 = generateB0Row(lambda, nuValues, eta, r, tolerance, nMax, initialB0)

nuValues = nuValues(:);
nNu      = numel(nuValues);
b0       = zeros(1, nNu);

bvpOptions = bvpset('NMax', nMax, ...
    'RelTol', tolerance, 'AbsTol', tolerance);

nu  = nuValues(1);
sol = bvp5c( ...
    @(rr, y) buffer.ode(rr, y, lambda, nu, eta), ...
    @(ya, yb) buffer.bc(ya, yb, lambda), ...
    bvpinit(r, @(rr) buffer.guess(rr, initialB0, lambda, nu, eta)), ...
    bvpOptions);

y     = deval(sol, r);
b0(1) = y(1, 1);

for k = 2:nNu
    nu  = nuValues(k);
    sol = bvp5c( ...
        @(rr, y) buffer.ode(rr, y, lambda, nu, eta), ...
        @(ya, yb) buffer.bc(ya, yb, lambda), ...
        sol, bvpOptions);

    y     = deval(sol, r);
    b0(k) = y(1, 1);
end
end
