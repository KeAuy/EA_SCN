function b = doubleExpProfile(r, q, alpha, eta)
%DOUBLEEXPPROFILE Double-exponential profile for one alpha per column.
%
%   B = BUFFER.DOUBLEEXPPROFILE(R, Q, ALPHA, ETA) evaluates
%
%       b = 1 + q (exp(-alpha r) - 1)/r
%           - eta q^3 (exp(-alpha r)(1 + alpha r) - 1)/r^2.
%
%   Q and ALPHA are row vectors, R is a column vector.
%
%   See also BUFFER.DOUBLEEXPSERIES, BUFFER.DOUBLEEXPVARIANT.

r     = r(:);
q     = reshape(q, 1, []);
alpha = reshape(alpha, 1, []);

decay = exp(-alpha .* r);
b = 1 + q .* (decay - 1) ./ r ...
    - eta * q.^3 .* (decay .* (1 + alpha .* r) - 1) ./ r.^2;
end
