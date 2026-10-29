function F = compute_F_Fig2(b0, b, r, hLeft, hRight, N, lambda, nu, eta)
%COMPUTE_F_FIG2  Residual of the discretised radial buffer-gas equation.
%
%   F = COMPUTE_F_FIG2(B0, B, R, HLEFT, HRIGHT, N, LAMBDA, NU, ETA) evaluates
%   the finite-difference residual of
%
%       lambda*(b'' + 2*b'/r) = (b - 1)*(nu*b + eta) + b/r
%
%   on the (non-uniform) grid R(1:N). HLEFT(k) = R(k) - R(k-1) and
%   HRIGHT(k) = R(k+1) - R(k) are the local grid spacings; the scheme is
%   second-order accurate in s = log(r).
%
%   Boundary conditions:
%       F(1) = b(1) - b0                                  (inner cut-off)
%       F(N) = b(N) - (1 - q/R(N) + eta*q^3/R(N)^2)       (far field, q = 1/(nu+eta))
%
%   NOTE: imposing b(N) = 1 at a finite rMax is only first-order accurate and
%   lifts the whole error curve (it removes the dip of Fig. 2C). If b(N) = 1
%   is preferred, choose rMax >= 1e6.
%
%   See also COMPUTE_J_FIG2.

q = 1/(eta + nu);

F    = zeros(N, 1);
F(1) = b(1) - b0;

for k = 2:N-1
    hMean = 0.5*(hLeft(k) + hRight(k));
    drSum = hLeft(k) + hRight(k);

    % Second derivative (central, non-uniform grid).
    d2b = (b(k+1) - b(k))/hRight(k) - (b(k) - b(k-1))/hLeft(k);

    F(k) = lambda*d2b/hMean ...
         + 2*lambda*(b(k+1) - b(k-1))/(drSum*r(k)) ...
         - (b(k) - 1)*(nu*b(k) + eta) ...
         - b(k)/r(k);
end

F(N) = b(N) - (1 - q/r(N) + eta*q^3/r(N)^2);

end
