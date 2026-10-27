function [b, info] = solve_buffer_bvp(b0, lambda, nu, eta, r, options)
%SOLVE_BUFFER_BVP Solve the buffer boundary value problem by Newton iteration.
%
%   B = SOLVE_BUFFER_BVP(B0, LAMBDA, NU, ETA, R) solves
%
%       b'' + (2/r) b' - (1/lambda) (b - 1) (nu b + eta) - b / (lambda r) = 0
%
%   on the uniform grid R with b(r(1)) = B0 and the truncated far-field
%   condition b(r(end)) = 1 - q/r(end) + ... . The nonlinear system obtained
%   from second-order central differences is solved by Newton's method with
%   the exact tridiagonal Jacobian.
%
%   B = SOLVE_BUFFER_BVP(..., 'SwitchRadius', R_SWITCH) selects the radius at
%   which the initial profile switches from the small-r to the large-r
%   expansion (default 1).
%
%   B = SOLVE_BUFFER_BVP(..., 'Tolerance', TOL) sets the convergence
%   threshold on the infinity norm of the residual (default 1e-10).
%
%   B = SOLVE_BUFFER_BVP(..., 'MaxIterations', N) sets the iteration limit
%   (default 25).
%
%   B = SOLVE_BUFFER_BVP(..., 'Verbose', FALSE) suppresses the iteration log
%   (default TRUE).
%
%   [B, INFO] = SOLVE_BUFFER_BVP(...) also returns INFO with the fields
%   ITERATIONS, RESIDUAL, CONVERGED and RESIDUAL_HISTORY.
%
%   See also BUFFER_EXPANSION, BUFFER_RESIDUAL_JACOBIAN.

arguments
    b0     (1,1) double
    lambda (1,1) double
    nu     (1,1) double
    eta    (1,1) double
    r      (:,1) double
    options.SwitchRadius  (1,1) double  = 1
    options.Tolerance     (1,1) double  = 1e-10
    options.MaxIterations (1,1) double  = 25
    options.Verbose       (1,1) logical = true
end

b = initial_profile(r, b0, lambda, nu, eta, options.SwitchRadius);

residualHistory = nan(options.MaxIterations, 1);
converged = false;

for iteration = 1:options.MaxIterations
    [residual, jacobian] = buffer_residual_jacobian(b, r, lambda, nu, eta, b0);
    b = b - jacobian \ residual;

    residualNorm = norm(residual, inf);
    residualHistory(iteration) = residualNorm;

    if options.Verbose
        fprintf('iteration %3d: ||F||_inf = %.3e\n', iteration, residualNorm);
    end

    if residualNorm < options.Tolerance
        converged = true;
        if options.Verbose
            fprintf('converged after %d iterations\n', iteration);
        end
        break
    end
end

info = struct( ...
    'iterations',      iteration, ...
    'residual',        residualNorm, ...
    'converged',       converged, ...
    'residualHistory', residualHistory(1:iteration));
end

function b = initial_profile(r, b0, lambda, nu, eta, switchRadius)
%INITIAL_PROFILE Piecewise initial profile for the Newton iteration.
%   Uses the small-r expansion below the switch radius and the large-r
%   expansion above it.
[bInner, ~, bOuter, ~] = buffer_expansion(r, b0, lambda, nu, eta);

b = bOuter;
inner = r < switchRadius;
b(inner) = bInner(inner);
end
