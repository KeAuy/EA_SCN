function residual = bc(ya, yb, lambda, eta, nu)
%BC Two-point boundary conditions of the buffer equation.
%
%   RESIDUAL = BUFFER.BC(YA, YB, LAMBDA) imposes
%
%       b'(r_min) = b(r_min)/(2 lambda),    b(r_max) = 1.
%
%   RESIDUAL = BUFFER.BC(YA, YB, LAMBDA, ETA, NU) replaces the outer
%   condition by the two-term far-field expansion
%
%       b(r_max) = 1 - q/r_max + eta q^3/r_max^2,   q = 1/(nu + eta).
%
%   See also BUFFER.ODE, BUFFER.GUESS, BUFFER.GENERATEB0.

if nargin < 5
    residual = [ya(2) - ya(1) / (2 * lambda);
                yb(1) - 1];
else
    q    = 1 / (eta + nu);
    rMax = 1e2;
    residual = [ya(2) - ya(1) / (2 * lambda);
                yb(1) - (1 - q / rMax + eta * q^3 / rMax^2)];
end
end
