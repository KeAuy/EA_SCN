function [b, B1, B2] = pade2(r, lambda, nu, eta, options)
%PADE2 [2/2] Pade approximation of the buffer solution.
%
%   B = BUFFER.PADE2(R, LAMBDA, NU, ETA, OPTIONS) evaluates
%
%       b = (r^2 + a1 r + a2) / (r^2 + b1 r + b2).
%
%   OPTIONS:
%       Refine             refine b1, b2 with FSOLVE (default false)
%       RealCoefficients   drop the imaginary part of b1, b2 (default true)
%
%   Figure 1 uses the raw algebraic coefficients; Figure 2/4/5 use the
%   FSOLVE-refined coefficients because the paper's algebraic system is not
%   exactly satisfied by the closed-form construction.  Figure 3 uses the
%   raw coefficients with their real part.
%
%   See also BUFFER.APPROXIMATIONS, BUFFER.PADE1.

arguments
    r double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    nu double
    eta (1,1) double
    options.Refine (1,1) logical = false
    options.RefineMethod (1,:) char ...
        {mustBeMember(options.RefineMethod, {'fsolve', 'newton'})} = 'fsolve'
    options.RealCoefficients (1,1) logical = true
end

r  = r(:);
nu = nu(:)';
q  = 1 ./ (eta + nu);
p  = q .* nu;

cR = 24 * lambda - 3 * p .* q .* (p - 2);
cK = 8 * lambda^2 ...
     - 2 * lambda * q .* (2 * p.^2 - 5 * p - 2) ...
     - p .* q.^2 .* (p - 2);
cQ = cK + 24 * lambda^2 - 4 * lambda * p .* q .* (2 * p - 3);
cV = 2 * lambda * (p + 6) + p .* q .* (p - 3) .* (p - 2);
cF = 6 * lambda^2 ...
     - lambda * q .* (p + 3) .* (p - 2) ...
     - eta * p .* q.^3 .* (2 * p - 3);
cJ = 24 * lambda^2 ...
     - 2 * lambda * p .* q .* (5 * p - 8) ...
     - eta * p.^2 .* q.^3 .* (p - 2);
cE = 12 * lambda^2 + eta * p .* q.^2 .* (8 * lambda + eta * p .* q.^2);
cH = p.^2 .* q.^4 .* (p - 2).^2 + 1600 * lambda^4 ...
     - 104 * lambda^3 .* q .* (p - 2) .* (9 * p + 4) ...
     + 4 * lambda^2 .* q.^2 ...
       .* (39 * p.^4 - 87 * p.^3 - 63 * p.^2 + 184 * p + 4) ...
     - 4 * lambda * p .* q.^3 .* (p - 2) ...
       .* (6 * p.^3 - 21 * p.^2 + 20 * p + 2);
cW = 12000 * lambda^5 ...
     - 4 * lambda^4 .* q .* (2687 * p.^2 - 2966 * p - 553) ...
     + 4 * lambda^3 .* q.^2 ...
       .* (729 * p.^4 - 1626 * p.^3 + 486 * p.^2 + 421 * p + 70) ...
     - lambda^2 .* q.^3 ...
       .* (243 * p.^6 - 810 * p.^5 + 609 * p.^4 ...
          + 210 * p.^3 - 21 * p.^2 - 240 * p - 4) ...
     - 2 * lambda * eta * p .* q.^5 .* (p - 2) .* (20 * p.^2 - 13 * p - 1) ...
     + p.^2 .* eta^2 .* q.^7;

cG = cQ.^3 + 2 * lambda * cR ...
     .* (3 * cQ .* cF - 9 * lambda * cE .* cR + cV .* sqrt(-3 * q .* cW));

B2raw = (q ./ cR) .* (cQ + cG.^(1/3) + cH ./ cG.^(1/3));

if options.Refine
    B2guess = real(B2raw);
    B1guess = (1 ./ cV) ...
        .* ((B2guess .* (B2guess .* cR - 3 * q .* cK)) ...
            ./ (6 * lambda * q.^2) - cJ);
    B1guess = real(B1guess);

    if strcmpi(options.RefineMethod, 'fsolve')
        B1 = zeros(size(q));
        B2 = zeros(size(q));
        fsolveOptions = optimoptions('fsolve', 'Display', 'none', ...
            'TolFun', 1e-14, 'TolX', 1e-14, 'MaxIterations', 500);

        for i = 1:numel(q)
            qj = q(i);
            padeSystem = @(coef) [ ...
                2 * lambda * coef(1)^2 + coef(1) * coef(2) - coef(2)^2 / qj ...
                - 2 * lambda * qj^2 * coef(1) - (qj^2 + 2 * lambda) * coef(2);
                (1 + qj) * coef(1)^2 - 2 * coef(1) * coef(2) / qj ...
                + 2 * (3 * lambda - qj^3) * coef(1) + (1 - qj) * coef(2) ...
                - qj^2 * (6 * lambda + qj^2 * (1 - qj))];

            initialGuess = [B1guess(i); B2guess(i)];
            [coef, ~, exitFlag] = fsolve(padeSystem, initialGuess, fsolveOptions);

            if exitFlag <= 0 && i > 1
                initialGuess = [B1guess(i-1); B2guess(i-1)];
                coef = fsolve(padeSystem, initialGuess, fsolveOptions);
            end

            B1(i) = coef(1);
            B2(i) = coef(2);
        end
    else
        % Vectorized Newton refinement over all nu columns.  This is the
        % refinement used by the lambda-nu parameter-grid figures.
        B1 = B1guess;
        B2 = B2guess;

        for itP = 1:40
            fP1 = 2 * lambda * B1.^2 + B1 .* B2 - B2.^2 ./ q ...
                  - 2 * lambda * q.^2 .* B1 ...
                  - (q.^2 + 2 * lambda) .* B2;
            fP2 = (1 + q) .* B1.^2 - 2 * B1 .* B2 ./ q ...
                  + 2 * (3 * lambda - q.^3) .* B1 ...
                  + (1 - q) .* B2 ...
                  - q.^2 .* (6 * lambda + q.^2 .* (1 - q));
            resP = max(abs([fP1(:); fP2(:)]));
            if ~isfinite(resP) || resP < 1e-10
                break;
            end

            jP11 = 4 * lambda .* B1 + B2 - 2 * lambda * q.^2;
            jP12 = B1 - 2 * B2 ./ q - (q.^2 + 2 * lambda);
            jP21 = 2 * (1 + q) .* B1 - 2 * B2 ./ q ...
                   + 2 * (3 * lambda - q.^3);
            jP22 = -2 * B1 ./ q + (1 - q);
            detP = jP11 .* jP22 - jP12 .* jP21;

            B1 = B1 + (-fP1 .* jP22 + fP2 .* jP12) ./ detP;
            B2 = B2 + (-fP2 .* jP11 + fP1 .* jP21) ./ detP;
        end
    end
else
    B2 = B2raw;
    if options.RealCoefficients
        B2 = real(B2);
    end

    B1 = (1 ./ cV) ...
        .* ((B2 .* (B2 .* cR - 3 * q .* cK)) ./ (6 * lambda * q.^2) - cJ);
    if options.RealCoefficients
        B1 = real(B1);
    end
end

if options.Refine && options.RealCoefficients
    B1 = real(B1);
    B2 = real(B2);
end

A1 = B1 - q;
A2 = B2 - q .* (B1 - eta * q.^2);

b = (r.^2 + A1 .* r + A2) ./ (r.^2 + B1 .* r + B2);
end
