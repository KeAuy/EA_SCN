%EXPSERIES Single-exponential series 1 + q (exp(-alpha r) - 1)/r.
%
%   ALPHA = (sqrt(S^2 + lambda/q) - S)/lambda.  S is a scalar (Linear and
%   ExpSer panels) or one value per nu (ExpVar and ExpGlobal panels).

function b = expSeries(r, lambda, nu, eta, S)

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
