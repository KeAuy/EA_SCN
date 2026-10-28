function [b, info] = solveUniform(b0, lambda, nu, eta, r, options)
%SOLVEUNIFORM Solve the buffer BVP on a uniform grid by Newton iteration.
%
%   [B, INFO] = BUFFER.SOLVEUNIFORM(B0, LAMBDA, NU, ETA, R, OPTIONS) solves
%
%       b'' + (2/r) b' - (1/lambda)(b - 1)(nu b + eta) - b/(lambda r) = 0
%
%   on the uniform grid R with b(r(1)) = B0 and the four-term far-field
%   condition at r(end).
%
%   OPTIONS:
%       SwitchRadius   fixed switch radius for the initial profile, or []
%                      to choose the point where the expansions are closest
%       Tolerance      residual infinity norm target (default 1e-10)
%       MaxIterations  iteration limit (default 25)
%       Verbose        print iteration progress (default true)
%
%   INFO has fields ITERATIONS, RESIDUAL, CONVERGED and RESIDUALHISTORY.
%
%   See also BUFFER.RESIDUALJACOBIANUNIFORM, BUFFER.INITIALPROFILE.

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
