rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);

fprintf('Test: approximations ... ');
r  = (1e-3:1e-2:1e2)';
nu = [1e-3; 1; 100];
names = {'Pade', 'Pade2', 'RBA2', 'Linear', 'ExpSer', 'ExpVar', ...
         'ExpGlobal', 'ExpPade', 'DblExpSer', 'DblExpVar', 'DblExpGlobal'};
approx = buffer.approximations(r, 0.1, nu, 1, names, ...
    'PadeRefine', true, 'PadeRealCoefficients', true, 'RealBranchIndex', 2);
for k = 1:numel(names)
    field = lower(names{k});
    switch field
        case 'pade',            value = approx.pade;
        case 'pade2',           value = approx.pade2;
        case 'rba2',            value = approx.rba2;
        case 'linear',          value = approx.linear;
        case 'expser',          value = approx.expSeries;
        case 'expvar',          value = approx.expVar;
        case 'expglobal',       value = approx.expGlobal;
        case 'exppade',         value = approx.expPade;
        case 'dblexpser',       value = approx.doubleExpSeries;
        case 'dblexpvar',       value = approx.doubleExpVar;
        case 'dblexpglobal',    value = approx.doubleExpGlobal;
    end
    assert(isequal(size(value), [numel(r), numel(nu)]), ...
        'Unexpected size for %s.', names{k});
    assert(all(isfinite(real(value)), 'all'), ...
        'Non-finite real part for %s.', names{k});
end

fprintf('ok\nTest: Pade2 coefficient modes ... ');
a1 = buffer.pade2(r, 0.1, 0.1, 1, 'Refine', false, ...
    'RealCoefficients', false);
a2 = buffer.pade2(r, 0.1, 0.1, 1, 'Refine', true, ...
    'RealCoefficients', true);
assert(all(isfinite(real(a1))) && all(isfinite(real(a2))));
fprintf('ok\nTest: uniform solver ... ');

rUniform = linspace(1e-3, 1e2, 2000)';
b = buffer.solveUniform(0.131, 0.1, 0.1, 1, rUniform, ...
    'SwitchRadius', 3, 'Tolerance', 1e-8, 'MaxIterations', 50, ...
    'Verbose', false);
[F, ~] = buffer.residualJacobianUniform(b, rUniform, 0.131, 0.1, 0.1, 1);
assert(norm(F, inf) < 1e-6, 'Uniform solver residual is too large.');

fprintf('ok\nTest: geometric solver ... ');
nuSmall = logspace(-3, 2, 100)';
grid = buffer.geometricGrid(1e-3, 1e2, 100);
b0Small = buffer.generateB0(0.02, 1, nuSmall, ...
    'Tolerance', 1e-6, 'Verbose', false);
cols = 1:3;
bGeo = buffer.solveGeometric(b0Small(cols), 0.02, nuSmall(cols), 1, grid, ...
    'Tolerance', 1e-4, 'MaxIterations', 500, 'Verbose', false);
assert(all(isfinite(bGeo), 'all'), 'Geometric solver returned non-finite values.');
for k = 1:numel(cols)
    F = buffer.residualGeometric(bGeo(:, k), grid.r, grid.hLeft, ...
        grid.hRight, b0Small(k), 0.02, nuSmall(k), 1);
    assert(norm(F, inf) < 1e-3, 'Geometric solver residual is too large.');
end

fprintf('ok\nTest: charge density and errors ... ');
cNum = buffer.chargeDensity(bGeo, grid.r, nuSmall(cols));
err  = buffer.meanError(cNum, cNum, 'log');
assert(numel(err) == numel(cols) && all(isfinite(err)));

fprintf('ok\nTest: lambda-batch solver ... ');
lambdaBatch = logspace(-3, 2, 100)';
nuBatch     = logspace(-3, 2, 100)';
gridBatch   = buffer.geometricGrid(1e-3, 1e2, 100);
b0Batch = buffer.generateB0Grid(lambdaBatch, nuBatch, 1, ...
    'GridSize', 60, 'Tolerance', 1e-6, 'Verbose', false);
bBatch = buffer.solveGeometricLambdaGrid( ...
    b0Batch, lambdaBatch, nuBatch, 1, gridBatch, ...
    'Tolerance', 1e-3, 'MaxIterations', 200, 'Verbose', false);
assert(isequal(size(bBatch), [100, 100, 100]));
assert(all(isfinite(bBatch), 'all'));
F = buffer.residualGeometricLambda(bBatch(:, :, 1), gridBatch.r, ...
    gridBatch.hLeft, gridBatch.hRight, b0Batch(:, 1).', ...
    lambdaBatch, nuBatch(1), 1);
assert(max(abs(F), [], 'all') < 1e-4);

approxBatch = buffer.approximationsGrid(gridBatch.r, lambdaBatch, nuBatch, 1, ...
    {'Pade2', 'RBA2', 'DblExpGlobal'}, ...
    'PadeRefine', true, 'PadeRefineMethod', 'newton', ...
    'PadeRealCoefficients', true);
assert(isequal(size(approxBatch.pade2), [100,100,100]));
fprintf('ok\n');
fprintf('Test: all checks passed.\n');
