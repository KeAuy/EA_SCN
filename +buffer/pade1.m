function b = pade1(r, lambda, nu, eta)
%PADE1 [1/1] Pade approximation of the buffer solution.
%
%   B = BUFFER.PADE1(R, LAMBDA, NU, ETA) evaluates
%
%       b = 1 - q / (r + (sqrt(q (q + 8 lambda)) + q)/2),   q = 1/(nu+eta).
%
%   See also BUFFER.APPROXIMATIONS, BUFFER.PADE2.

r = r(:);
q = 1 ./ (eta + nu(:)');
b = 1 - q ./ (r + (sqrt(q .* (q + 8 * lambda)) + q) / 2);
end
