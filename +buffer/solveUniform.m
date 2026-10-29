%SOLVEUNIFORM Newton solver for the equidistant grid of Figures 1 and 3.
%
%   [B, INFO] = BUFFER.SOLVEUNIFORM(B0, LAMBDA, NU, ETA, R, OPTIONS) requires
%   an equidistant R.  OPTIONS: SwitchRadius, Tolerance, MaxIterations,
%   Verbose.  INFO stores iterations, residual and the residual history.

function [b, info] = solveUniform(b0, lambda, nu, eta, r, options)

arguments
    b0 (1,1) double
    lambda (1,1) double {mustBePositive}
    nu (1,1) double
    eta (1,1) double
    r (:,1) double {mustBePositive}
    options.SwitchRadius double = []
    options.Tolerance (1,1) double {mustBePositive} = 1e-10
    options.MaxIterations (1,1) double {mustBeInteger, mustBePositive} = 25
    options.Verbose (1,1) logical = true
end

if numel(r) < 3
    error('buffer:invalidGrid', 'The grid needs at least three points.');
end

h = diff(r);
if max(abs(h - mean(h))) > 1e-9 * mean(h)
    error('buffer:nonUniformGrid', ...
        'BUFFER.SOLVEUNIFORM requires an equidistant grid.');
end

b = buffer.initialProfile(r, b0, lambda, nu, eta, options.SwitchRadius);

residualHistory = nan(options.MaxIterations, 1);
converged       = false;
iteration       = 0;
residualNorm    = inf;

for iter = 1:options.MaxIterations
    iteration = iter;
    [F, J] = buffer.residualJacobianUniform(b, r, b0, lambda, nu, eta);
    b = b - J \ F;

    residualNorm = norm(F, inf);
    residualHistory(iter) = residualNorm;

    if options.Verbose
        fprintf('iteration %3d: ||F||_inf = %.3e\n', iter, residualNorm);
    end

    if residualNorm < options.Tolerance
        converged = true;
        if options.Verbose
            fprintf('converged after %d iterations\n', iter);
        end
        break;
    end
end

info = struct( ...
    'iterations',      iteration, ...
    'residual',        residualNorm, ...
    'converged',       converged, ...
    'residualHistory', residualHistory(1:iteration));
end
