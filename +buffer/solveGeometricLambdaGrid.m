%SOLVEGEOMETRICLAMBDA GRID Batched Newton solver on a lambda-nu grid.
%
%   [BNUM, INFO] = BUFFER.SOLVEGEOMETRICLAMBDAGRID(B0, LAMBDAVALUES, NUVALUES,
%   ETA, GRID, OPTIONS) iterates all lambda values of one nu simultaneously:
%   residual and Jacobian are N-by-NLAMBDA arrays and the linear systems are
%   solved column-wise by BUFFER.SOLVETRIDIAGONALBATCH, with a per-column
%   backtracking line search.
%
%   OPTIONS: Tolerance, StepTolerance, MaxIterations, Verbose.

function [bNum, info] = solveGeometricLambdaGrid(b0, lambdaValues, nuValues, eta, grid, options)

arguments
    b0 double
    lambdaValues double
    nuValues double
    eta (1,1) double
    grid struct
    options.Tolerance (1,1) double {mustBePositive} = 1e-5
    options.StepTolerance (1,1) double {mustBePositive} = 1e-12
    options.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 500
    options.Verbose (1,1) logical = true
end

lambdaValues = lambdaValues(:);
nuValues     = nuValues(:);
nLambda      = numel(lambdaValues);
nNu          = numel(nuValues);

if ~isequal(size(b0), [nLambda, nNu])
    error('buffer:sizeMismatch', ...
        'B0 must be NLAMBDA-by-NNU.');
end

r  = grid.r;
hL = grid.hLeft;
hR = grid.hRight;

bNum = zeros(numel(r), nLambda, nNu);
info = struct( ...
    'iterations', zeros(1, nNu), ...
    'residual',   zeros(nLambda, nNu), ...
    'step',       zeros(1, nNu), ...
    'converged',  false(1, nNu));

for j = 1:nNu
    b0j = b0(:, j).';
    b   = buffer.initialProfileLambda(r, b0j, lambdaValues, ...
        nuValues(j), eta);

    step      = inf(1, nLambda);
    converged = false;
    iteration = 0;

    for iter = 1:options.MaxIterations
        iteration = iter;
        F = buffer.residualGeometricLambda(b, r, hL, hR, b0j, ...
            lambdaValues, nuValues(j), eta);
        residual = max(abs(F), [], 1);
        if max(residual) < options.Tolerance
            converged = true;
            break;
        end

        [lower, main, upper] = buffer.jacobianGeometricLambda(b, r, hL, hR, ...
            lambdaValues, nuValues(j), eta);
        delta = buffer.solveTridiagonalBatch(lower, main, upper, -F);

        step = max(abs(delta), [], 1) ./ (1 + max(abs(b), [], 1));
        if all(step < options.StepTolerance)
            converged = true;
            break;
        end

        factors = ones(1, nLambda);
        for ls = 1:options.MaxIterations
            bTrial  = b + factors .* delta;
            FTrial  = buffer.residualGeometricLambda(bTrial, r, hL, hR, ...
                b0j, lambdaValues, nuValues(j), eta);
            rTrial  = max(abs(FTrial), [], 1);
            bad     = (step >= options.StepTolerance) ...
                      & (rTrial > (1 - 1e-4 * factors) .* residual);
            if ~any(bad)
                break;
            end
            factors(bad) = 0.5 * factors(bad);
        end

        b = b + factors .* delta;
    end

    bNum(:, :, j) = b;

    F = buffer.residualGeometricLambda(b, r, hL, hR, b0j, ...
        lambdaValues, nuValues(j), eta);
    residual = max(abs(F), [], 1);

    info.iterations(j) = iteration;
    info.residual(:, j) = residual;
    info.step(j)        = max(step);
    info.converged(j)   = converged;

    if options.Verbose
        fprintf('nu = %.3e completed: Newton %2d, max|F| = %.2e, max step %.1e\n', ...
            nuValues(j), iteration, max(residual), max(step));
    end
end
end
