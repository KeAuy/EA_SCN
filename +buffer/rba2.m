%RBA2 Rational approximation RBA2.
%
%   The first term is the rational boundary-layer approximation, the second
%   carries the O(lambda) correction.

function b = rba2(r, lambda, nu, eta)

r  = r(:);
nu = nu(:)';
q  = 1 ./ (eta + nu);

disc = (1 + q ./ r).^2 - 4 * nu .* q.^2 ./ r;
bRpa = 1 - (1 ./ (2 * nu .* q)) .* (1 + q ./ r - sqrt(disc));

b = bRpa + 2 * lambda * eta ./ (((1 + r ./ q).^2 - 4 * nu .* r).^2);
end
