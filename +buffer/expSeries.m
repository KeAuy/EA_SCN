function b = expSeries(r, lambda, nu, eta, S)
%EXPSERIES Single exponential series of the buffer solution.
%
%   B = BUFFER.EXPSERIES(R, LAMBDA, NU, ETA, S) evaluates
%
%       b = 1 + q (exp(-alpha r) - 1)/r,
%       alpha = (sqrt(S^2 + lambda/q) - S)/lambda,
%
%   with q = 1/(nu + eta).  S may be a scalar or one value per nu.  S = 0
%   gives the linearised solution, S = 1/2 the Exp-Ser approximation, and
%   the Exp-Var / Exp-Global variants use their own S(q).
%
%   See also BUFFER.APPROXIMATIONS.

r  = r(:);
nu = nu(:)';
q  = 1 ./ (eta + nu);

if isscalar(S)
    S = repmat(S, 1, numel(q));
else
    S = S(:).';
    if numel(S) ~= numel(q)
        error('buffer:sizeMismatch', ...
            'S must be scalar or have one entry per nu.');
    end
end

alpha = (sqrt(S.^2 + lambda ./ q) - S) ./ lambda;
b     = 1 + q .* (exp(-alpha .* r) - 1) ./ r;
end
