%SOLVEGEOMETRIC Newton solver for the geometric-grid discretization.
%
%   [BNUM, INFO] = BUFFER.SOLVEGEOMETRIC(B0VALUES, LAMBDA, NUVALUES, ETA, GRID,
%   OPTIONS) solves one column per nu with the exact tridiagonal Jacobian.  If
%   the plain Newton iteration stalls, a backtracking relaxation over
%   OPTIONS.RELAXATIONGRID is tried from the same initial profile.
%
%   OPTIONS: Tolerance, MaxIterations, RelaxationGrid, Verbose.

function [bNum, info] = solveGeometric(b0Values, lambda, nuValues, eta, grid, options)

arguments
    b0Values double
    lambda (1,1) double {mustBePositive}
    nuValues double
    eta (1,1) double
    grid struct
    options.Tolerance (1,1) double {mustBePositive} = 1e-8
    options.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 50
    options.RelaxationGrid (1,:) double = 0.1:0.05:1.0
    options.Verbose (1,1) logical = true
end

b0Values = b0Values(:);
nuValues = nuValues(:);
if numel(b0Values) ~= numel(nuValues)
    error('buffer:sizeMismatch', ...
        'b0Values and nuValues must have the same number of entries.');
end

r  = grid.r;
hL = grid.hLeft;
hR = grid.hRight;

bNum = buffer.initialProfile(r, b0Values, lambda, nuValues, eta, []);

nNu  = numel(nuValues);
info = struct( ...
    'iterations',     zeros(nNu, 1), ...
    'residual',       nan(nNu, 1), ...
    'converged',      false(nNu, 1), ...
    'usedRelaxation', false(nNu, 1));

for j = 1:nNu
    initial    = bNum(:, j);
    converged  = false;
    residualNorm = inf;
    iteration  = 0;

    for iter = 1:options.MaxIterations
        iteration = iter;
        F = buffer.residualGeometric(bNum(:, j), r, hL, hR, ...
            b0Values(j), lambda, nuValues(j), eta);
        residualNorm = norm(F, inf);
        if residualNorm < options.Tolerance
            converged = true;
            break;
        end

        J = buffer.jacobianGeometric(bNum(:, j), r, hL, hR, ...
            lambda, nuValues(j), eta);
        bNum(:, j) = bNum(:, j) - J \ F;
    end

    if converged
        if options.Verbose
            fprintf('  Newton converged, j = %d\n', j);
        end
    else
        if options.Verbose
            fprintf('  Newton failed at j = %d, switching to relaxation.\n', j);
        end

        [bNum(:, j), converged, iteration, residualNorm] = relaxColumn( ...
            initial, r, hL, hR, b0Values(j), lambda, nuValues(j), eta, ...
            options, j);

        info.usedRelaxation(j) = true;
    end

    info.iterations(j) = iteration;
    info.residual(j)   = residualNorm;
    info.converged(j)  = converged;

    if ~converged
        warning('buffer:notConverged', ...
            'Column %d (nu = %.3e) did not converge.', j, nuValues(j));
    end
end
end

function [b, converged, iteration, residualNorm] = relaxColumn( ...
        b, r, hL, hR, b0, lambda, nu, eta, options, column)

converged    = false;
residualNorm = inf;
iteration    = 0;

for iter = 1:options.MaxIterations
    iteration = iter;
    F = buffer.residualGeometric(b, r, hL, hR, b0, lambda, nu, eta);
    resOld = norm(F, inf);
    if resOld < options.Tolerance
        converged    = true;
        residualNorm = resOld;
        if options.Verbose
            fprintf('  Relaxation converged, j = %d\n', column);
        end
        return;
    end

    J     = buffer.jacobianGeometric(b, r, hL, hR, lambda, nu, eta);
    delta = -J \ F;

    bestRes = inf;
    bestB   = b;
    for omega = options.RelaxationGrid
        bTrial   = b + omega * delta;
        FTrial   = buffer.residualGeometric(bTrial, r, hL, hR, ...
            b0, lambda, nu, eta);
        resTrial = norm(FTrial, inf);
        if resTrial < bestRes
            bestRes = resTrial;
            bestB   = bTrial;
        end
    end

    if bestRes < resOld
        b            = bestB;
        residualNorm = bestRes;
    else
        residualNorm = resOld;
        if options.Verbose
            fprintf('  Relaxation failed at j = %d, iter = %d\n', column, iter);
        end
        return;
    end
end
end
