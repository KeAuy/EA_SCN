function F = compute_F_Fig2(b0, b, r, hLeft, hRight, N, lambda, nu, eta)

q = 1/(eta + nu);

F    = zeros(N, 1);
F(1) = b(1) - b0;

for k = 2:N-1
    hMean = 0.5*(hLeft(k) + hRight(k));
    drSum = hLeft(k) + hRight(k);

    d2b = (b(k+1) - b(k))/hRight(k) - (b(k) - b(k-1))/hLeft(k);

    F(k) = lambda*d2b/hMean ...
         + 2*lambda*(b(k+1) - b(k-1))/(drSum*r(k)) ...
         - (b(k) - 1)*(nu*b(k) + eta) ...
         - b(k)/r(k);
end

F(N) = b(N) - (1 - q/r(N) + eta*q^3/r(N)^2);

end
