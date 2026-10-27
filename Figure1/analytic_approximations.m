function approx = analytic_approximations(r, lambda, nu, eta)
%ANALYTIC_APPROXIMATIONS Closed-form approximations of the buffer solution.
%
%   APPROX = ANALYTIC_APPROXIMATIONS(R, LAMBDA, NU, ETA) evaluates the
%   closed-form approximations of b(r) used in Figure 1 at the radii R:
%
%       APPROX.PADE              [1/1] Pade approximation
%       APPROX.PADE2             [2/2] Pade approximation
%       APPROX.EXP_SERIES        exponential series
%       APPROX.EXP_PADE          exponential-Pade approximation
%       APPROX.DOUBLE_EXP_SERIES double exponential series
%       APPROX.RBA2              rational approximation RBA2
%       APPROX.LINEARIZED        solution of the linearised equation
%
%   All fields have the same shape as R. Some approximations are complex for
%   the parameters of Figure 1A; plot REAL(...) of those.
%
%   See also FIG1A, FIG1B, FIG1C.

q = 1 / (eta + nu);

approx = struct();
approx.pade              = pade_approximation(r, q, lambda);
approx.pade2             = pade2_approximation(r, q, lambda, nu, eta);
approx.exp_series        = exp_series_approximation(r, q, lambda);
approx.exp_pade          = exp_pade_approximation(r, q, lambda, eta);
approx.double_exp_series = double_exp_series_approximation(r, q, lambda, eta);
approx.rba2              = rba2_approximation(r, q, lambda, nu, eta);
approx.linearized        = 1 + q * (exp(-r / sqrt(q * lambda)) - 1) ./ r;
end

function b = pade_approximation(r, q, lambda)
%PADE_APPROXIMATION [1/1] Pade approximation.
%   Matched to the near-field and far-field expansions of b(r).
b = 1 - q ./ (r + (sqrt(q * (q + 8 * lambda)) + q) / 2);
end

function b = pade2_approximation(r, q, lambda, nu, eta)
%PADE2_APPROXIMATION [2/2] Pade approximation.
%   b(r) = (r^2 + a1 r + a2) / (r^2 + b1 r + b2). The coefficients follow
%   from matching the far-field expansion; the intermediate quantities below
%   (r0, k0, ... , g0) are the steps of that algebraic derivation.
p  = q * nu;
r0 = 24 * lambda - 3 * p * q * (p - 2);
k0 = 8 * lambda^2 - 2 * lambda * q * (2 * p^2 - 5 * p - 2) - p * q^2 * (p - 2);
q0 = k0 + 24 * lambda^2 - 4 * lambda * p * q * (2 * p - 3);
v0 = 2 * lambda * (p + 6) + p * q * (p - 3) * (p - 2);
f0 = 6 * lambda^2 - lambda * q * (p + 3) * (p - 2) - eta * p * q^3 * (2 * p - 3);
j0 = 24 * lambda^2 - 2 * lambda * p * q * (5 * p - 8) - eta * p^2 * q^3 * (p - 2);
e0 = 12 * lambda^2 + eta * p * q^2 * (8 * lambda + eta * p * q^2);
h0 = p^2 * q^4 * (p - 2)^2 + 1600 * lambda^4 ...
     - 104 * lambda^3 * q * (p - 2) * (9 * p + 4) ...
     + 4 * lambda^2 * q^2 * (39 * p^4 - 87 * p^3 - 63 * p^2 + 184 * p + 4) ...
     - 4 * lambda * p * q^3 * (p - 2) * (6 * p^3 - 21 * p^2 + 20 * p + 2);
w0 = 12000 * lambda^5 ...
     - 4 * lambda^4 * q * (2687 * p^2 - 2966 * p - 553) ...
     + 4 * lambda^3 * q^2 * (729 * p^4 - 1626 * p^3 + 486 * p^2 + 421 * p + 70) ...
     - lambda^2 * q^3 * (243 * p^6 - 810 * p^5 + 609 * p^4 + 210 * p^3 - 21 * p^2 - 240 * p - 4) ...
     - 2 * lambda * eta * p * q^5 * (p - 2) * (20 * p^2 - 13 * p - 1) ...
     + p^2 * eta^2 * q^7;

g0 = q0^3 + 2 * lambda * r0 * (3 * q0 * f0 - 9 * lambda * e0 * r0 + v0 * sqrt(-3 * q * w0));

b2 = (q / r0) * (q0 + g0^(1/3) + h0 / g0^(1/3));
b1 = (1 / v0) * ((b2 * (b2 * r0 - 3 * q * k0)) / (6 * lambda * q^2) - j0);
a1 = b1 - q;
a2 = b2 - q * (b1 - eta * q^2);

b = (r.^2 + a1 * r + a2) ./ (r.^2 + b1 * r + b2);
end

function b = exp_series_approximation(r, q, lambda)
%EXP_SERIES_APPROXIMATION Exponential series.
%   b = 1 + q (exp(-alpha r) - 1) / r, with alpha fixed by S = 1/2.
S     = 1 / 2;
alpha = (sqrt(S^2 + lambda / q) - S) / lambda;

b = 1 + q * (exp(-alpha * r) - 1) ./ r;
end

function b = exp_pade_approximation(r, q, lambda, eta)
%EXP_PADE_APPROXIMATION Exponential-Pade approximation.
%   b = 1 + q (exp(-alpha r) - 1)/r + eta q^3 / (beta + r^2), where alpha and
%   beta follow from a cubic construction.
p0 = 6 + eta * q^2 * (1 - eta * q) / lambda;
q0 = (1 / p0) * (6 + eta * q^2 / (2 * lambda));
r0 = q0 * (1 - 2 * lambda / q) + 2 * eta^2 * q^2 / p0;
u0 = r0 / 3 - q0^2 / 2;
v0 = q0 * (r0 - q0^2) + 12 * lambda / (q * p0);
w0 = q0^2 * (r0 / 3 - q0^2 / 4) ...
     + 4 * lambda * (6 * (q0 + lambda / q) - eta^2 * q^2) / (3 * q * p0);
e0 = v0^2 * (v0^2 / 4 + u0 * (u0^2 - 3 * w0)) - w0 * (w0 - 3 * u0^2)^2;
g0 = v0^2 / 2 + u0 * (u0^2 - 3 * w0) + sqrt(e0);
h0 = (w0 + u0^2) / g0^(1/3) + g0^(1/3) - 2 * u0;

alpha = (1 / (2 * lambda)) * (sqrt(h0) - q0 + sqrt(2 * v0 / sqrt(h0) - h0 - 6 * u0));
beta  = eta * q^3 / (q * alpha * (lambda * alpha + 1) - 1);

b = 1 + q * (exp(-alpha * r) - 1) ./ r + q^3 * eta ./ (beta + r.^2);
end

function b = double_exp_series_approximation(r, q, lambda, eta)
%DOUBLE_EXP_SERIES_APPROXIMATION Double exponential series.
%   b = 1 + q (exp(-alpha r) - 1)/r
%       - eta q^3 (exp(-alpha r) (1 + alpha r) - 1)/r^2,
%   where alpha is the decaying root of the cubic below.
alpha = double_exp_series_alpha(q, lambda, eta);

% Drop a numerically vanishing imaginary part so that the result stays real
% for parameter sets with three real roots; genuinely complex alpha (as in
% Figure 1A) are kept.
if abs(imag(alpha)) < 1e-12 * abs(alpha)
    alpha = real(alpha);
end

decay = exp(-alpha * r);
b = 1 + q * (decay - 1) ./ r ...
    - eta * q^3 * (decay .* (1 + alpha * r) - 1) ./ r.^2;
end

function alpha = double_exp_series_alpha(q, lambda, eta)
%DOUBLE_EXP_SERIES_ALPHA Decaying root of the double-exponential cubic.
%   eta q^2 P alpha^3 - Q alpha^2 - R alpha + 1/q = 0 with
%   P = 2 lambda/3, Q = lambda - eta q^2/2 and R = 1.
p0 = 2 * lambda / 3;
q0 = lambda - eta * q^2 / 2;
r0 = 1;

w0 = 3 * eta * q^2 * p0;
e0 = q0^2 + w0 * r0;
f0 = (3 / 2) * (q0 * r0 - w0 / q) * w0 + q0^3;

% The three roots are obtained by rotating the principal cube root with the
% cube roots of unity; the first branch is the physically relevant one.
cubeRoots = [1/2 + sqrt(3)/2 * 1i, 1/2 - sqrt(3)/2 * 1i, -1];
g0 = cubeRoots(1) * (sqrt(f0^2 - e0^3) + f0)^(1/3);

alpha = (q0 - g0 - e0 / g0) / w0;
end

function b = rba2_approximation(r, q, lambda, nu, eta)
%RBA2_APPROXIMATION Rational approximation RBA2.
bRpa = 1 - (1 / (2 * q * nu)) * ...
       (1 + q ./ r - sqrt((1 + q ./ r).^2 - 4 * nu * q^2 ./ r));

b = bRpa + 2 * lambda * eta * ((1 + r / q).^2 - 4 * nu * r).^(-2);
end
