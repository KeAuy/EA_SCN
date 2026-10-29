function J = compute_J_Fig2(b, r, hLeft, hRight, N, lambda, nu, eta)

J = spalloc(N, N, 3*N);
J(1, 1) = 1;
J(N, N) = 1;

for k = 2:N-1
    hMean = 0.5*(hLeft(k) + hRight(k));
    drSum = hLeft(k) + hRight(k);
    rk    = r(k);

    J(k, k-1) = lambda/hLeft(k)/hMean - 2*lambda/(drSum*rk);
    J(k, k)   = (-lambda/hRight(k) - lambda/hLeft(k))/hMean ...
                - (2*nu*b(k) - nu + eta) ...
                - 1/rk;
    J(k, k+1) = lambda/hRight(k)/hMean + 2*lambda/(drSum*rk);
end

end
