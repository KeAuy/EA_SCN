%DOUBLEEXPPROFILE Double-exponential profile built from a given alpha
%
%       b = 1 + q (exp(-alpha r) - 1)/r
%             - eta q^3 (exp(-alpha r) (1 + alpha r) - 1)/r^2.

function b = doubleExpProfile(r, q, alpha, eta)

r     = r(:);
q     = reshape(q, 1, []);
alpha = reshape(alpha, 1, []);

decay = exp(-alpha .* r);
b = 1 + q .* (decay - 1) ./ r ...
    - eta * q.^3 .* (decay .* (1 + alpha .* r) - 1) ./ r.^2;
end
