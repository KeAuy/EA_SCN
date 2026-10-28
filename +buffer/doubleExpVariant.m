function b = doubleExpVariant(r, lambda, nu, eta, type)
%DOUBLEEXPVARIANT Exp-Var and Exp-Global double exponential approximations.
%
%   B = BUFFER.DOUBLEEXPVARIANT(R, LAMBDA, NU, ETA, TYPE) evaluates the
%   double exponential series with the integral-averaged coefficients used
%   in Figures 3/4/5.  TYPE is 'var' or 'global'.
%
%   See also BUFFER.DOUBLEEXPALPHA, BUFFER.DOUBLEEXPPROFILE.

arguments
    r double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    type (1,:) char {mustBeMember(type, {'var', 'global'})}
end

r = r(:);
q = 1 ./ (eta + nu(:)');

switch type
    case 'var'
        pCoef = lambda * (8 * log(2) - 5) ...
                + ((4 / 3) * q.^2 * eta) .* (1 - eta * q) ...
                  .* (1 - 3 * log(4 / 3));
        qCoef = lambda + ((2 / 3) * q.^2 * eta) ...
                .* (1 - 6 * log(9 / 8) ...
                    + (2 * q * eta) .* (1 - 6 * log(4 / 3)));
        rCoef = (q * eta + 2) / 3;
    case 'global'
        pCoef = 2 * lambda * (1 - log(2)) ...
                + q.^2 * eta .* (1 - q * eta) .* (log(3) - 1);
        qCoef = lambda - 2 * q.^2 * eta ...
                .* (1 - log(81 / 32) + 2 * q * eta * log(9 / 8));
        rCoef = q * eta + 2 * (1 - q * eta) * log(3 / 2);
end

alpha = buffer.doubleExpAlpha(q, eta, pCoef, qCoef, rCoef).';
b     = buffer.doubleExpProfile(r, q, alpha, eta);
end
