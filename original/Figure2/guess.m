function yinit = guess(b0, lambda, nu, eta, r)
%GUESS  Initial guess for BVP5C: blend of two asymptotic expansions.
%
%   YINIT = GUESS(B0, LAMBDA, NU, ETA, R) returns the state [b; db/dr] at
%   radius R, built from
%
%   the near-origin (Frobenius) expansion, r -> 0,
%       b   = b0 + c1*r + c2*r^2,
%       c1  = b0/(2*lambda),
%       c2  = ((b0 - 1)*(nu*b0 + eta) + c1)/(6*lambda),
%
%   and the far-field expansion, r -> Inf, with q = 1/(nu + eta),
%       b   = 1 - q/r + eta*q^3/r^2 + eta*(1 - 2*eta*q)*q^4/r^3.
%
%   The near-origin branch is used where the two expansions agree to within
%   10%, and also for r < 1; elsewhere the far-field branch is used.
%
%   See also ODEFUN, BCFUN.

q    = 1/(eta + nu);
c1   = b0/(2*lambda);
beta = ((b0 - 1)*(nu*b0 + eta) + c1)/(6*lambda);

% Near-origin expansion and its derivative.
bNear  = b0 + c1*r + beta*r^2;
dbNear = c1 + 2*beta*r;

% Far-field expansion in x = 1/r, then convert db/dx to db/dr.
x     = 1/r;
bFar  = 1 - q*x + eta*q^3*x^2 + eta*(1 - 2*eta*q)*q^4*x^3;
dbFar = q*x^2 - 2*eta*q^3*x^3 - 3*eta*(1 - 2*eta*q)*q^4*x^4;
dbFar = -dbFar/r^2;                 % chain rule: d/dr = -(1/r^2)*d/dx

relativeGap = abs(bNear - bFar)/(abs(bNear) + abs(bFar) + 1e-14);

if relativeGap < 0.1 || r < 1
    yinit = [bNear; dbNear];
else
    yinit = [bFar; dbFar];
end

end
