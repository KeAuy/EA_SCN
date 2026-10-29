function J = compute_J_Fig2(b, r, hLeft, hRight, N, lambda, nu, eta)
%COMPUTE_J_FIG2  Jacobian of the residual returned by COMPUTE_F_FIG2.
%
%   J = COMPUTE_J_FIG2(B, R, HLEFT, HRIGHT, N, LAMBDA, NU, ETA) assembles the
%   sparse N-by-N Jacobian of the finite-difference residual. The matrix is
%   tridiagonal apart from the two diagonal entries of the Dirichlet rows.
%
%   See also COMPUTE_F_FIG2.

J = spalloc(N, N, 3*N);
J(1, 1) = 1;
J(N, N) = 1;

for k = 2:N-1
    hMean = 0.5*(hLeft(k) + hRight(k));
    drSum = hLeft(k) + hRight(k);
    rk    = r(k);

    % d/db(k-1), d/db(k), d/db(k+1) of lambda*(b'' + 2*b'/r).
    J(k, k-1) = lambda/hLeft(k)/hMean - 2*lambda/(drSum*rk);
    J(k, k)   = (-lambda/hRight(k) - lambda/hLeft(k))/hMean ...
                - (2*nu*b(k) - nu + eta) ...   % d/db of -(b-1)*(nu*b+eta)
                - 1/rk;
    J(k, k+1) = lambda/hRight(k)/hMean + 2*lambda/(drSum*rk);
end

end
