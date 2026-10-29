function [pCoef, qCoef, rCoef, q] = doubleExpCoefficients(lambda, nu, eta, type)
%DOUBLEEXPCOEFFICIENTS Coefficients of the Exp-Var / Exp-Global cubic.
%
%   [P, Q, R, QV] = BUFFER.DOUBLEEXPCOEFFICIENTS(LAMBDA, NU, ETA, TYPE)
%   returns the cubic coefficients used by the double-exponential variant
%   approximations.  TYPE is 'var' or 'global'.  All outputs are column
%   vectors; QV = 1./(NU + ETA).
%
%   See also BUFFER.DOUBLEEXPVARIANT, BUFFER.DOUBLEEXPALPHA.

arguments
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    type (1,:) char {mustBeMember(type, {'var', 'global'})}
end

nu = nu(:);
q  = 1 ./ (eta + nu);

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
end
