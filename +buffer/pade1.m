%PADE1 Pade [1/1] approximation of b(r)
%
%       b = 1 - q/(r + (sqrt(q (q + 8 lambda)) + q)/2).

function b = pade1(r, lambda, nu, eta)

r = r(:);
q = 1 ./ (eta + nu(:)');
b = 1 - q ./ (r + (sqrt(q .* (q + 8 * lambda)) + q) / 2);
end
