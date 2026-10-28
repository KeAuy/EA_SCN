function b0 = generateB0(lambda, eta, nuValues, outputFile, options)
%GENERATEB0 Generate the inner boundary data b0(nu) with BVP5C.
%
%   B0 = BUFFER.GENERATEB0(LAMBDA, ETA, NUVALUES, OUTPUTFILE, OPTIONS)
%   solves the buffer BVP with BVP5C for every entry of NUVALUES on the
%   geometric grid used by the figure panels, continuing from the previous
%   solution.  The row vector B0 contains b(r_min) for each nu.
%
%   NUVALUES should be a reasonably resolved log-spaced grid; the
%   continuation assumes that neighbouring values are close enough for
%   BVP5C to converge.
%
%   OUTPUTFILE may be empty to skip saving.  OPTIONS:
%       GridSize   number of radial grid points (default NUMEL(NUVALUES))
%       Tolerance  RelTol and AbsTol for BVP5C (default 1e-8)
%       NMax       maximum number of mesh points (default 50000)
%       InitialB0  parameter of the initial BVP5C guess (default 0.5)
%       Verbose    print progress (default true)
%
%   See also BUFFER.ODE, BUFFER.BC, BUFFER.GUESS.

arguments
    lambda (1,1) double {mustBePositive}
    eta (1,1) double
    nuValues double
    outputFile = ''
    options.GridSize (1,1) double {mustBeInteger, mustBeGreaterThan(options.GridSize, 2)} = []
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

if ~isempty(outputFile)
    outputFile = char(outputFile);
    outDir     = fileparts(outputFile);
    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end
    save(outputFile, 'b0');
end
end
