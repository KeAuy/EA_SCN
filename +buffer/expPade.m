function b = expPade(r, lambda, nu, eta)
%EXPPADE Exponential-Pade approximation of the buffer solution.
%
%   B = BUFFER.EXPPADE(R, LAMBDA, NU, ETA) evaluates the exponential-Pade
%   approximation used in Figures 1, 2 and 5.  The result may be complex for
%   some parameter values; callers plot the real part or take ABS of errors.
%
%   See also BUFFER.APPROXIMATIONS.

r = r(:);
q = 1 ./ (eta + nu(:)');

pP = 6 + eta * q.^2 .* (1 - eta * q) / lambda;
pQ = (1 ./ pP) .* (6 + eta * q.^2 / (2 * lambda));
pR = pQ .* (1 - 2 * lambda ./ q) + 2 * eta^2 * q.^2 ./ pP;
pU = pR / 3 - pQ.^2 / 2;
pV = pQ .* (pR - pQ.^2) + 12 * lambda ./ (q .* pP);
pW = pQ.^2 .* (pR / 3 - pQ.^2 / 4) ...
     + 4 * lambda .* ((6 * (pQ + lambda ./ q) - eta^2 * q.^2) ...
                      ./ (3 * q .* pP));
pE = pV.^2 .* (pV.^2 / 4 + pU .* (pU.^2 - 3 * pW)) ...
     - pW .* (pW - 3 * pU.^2).^2;
pG = pV.^2 / 2 + pU .* (pU.^2 - 3 * pW) + sqrt(pE);
pH = (pW + pU.^2) ./ pG.^(1/3) + pG.^(1/3) - 2 * pU;

alpha = (sqrt(pH) - pQ ...
         + sqrt(2 * pV ./ sqrt(pH) - pH - 6 * pU)) / (2 * lambda);
beta  = eta * q.^3 ./ (q .* alpha .* (lambda * alpha + 1) - 1);

b = 1 + q .* (exp(-alpha .* r) - 1) ./ r + q.^3 * eta ./ (beta + r.^2);
end
