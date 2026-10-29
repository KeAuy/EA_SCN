function J = compute_J_Fig2(b, r, hL, hR, N, lambda, nu, eta)

J = spalloc(N, N, 3*N);
J(1,1) = 1;
J(N,N) = 1;

for k = 2:N-1
    hb = 0.5*(hL(k)+hR(k));
    J(k,k-1) = lambda/hL(k)/hb - 2*lambda/((hL(k)+hR(k))*r(k));
    J(k,k)   = (-lambda/hR(k) - lambda/hL(k))/hb ...
               - (2*nu*b(k) - nu + eta) - 1/r(k);
    J(k,k+1) = lambda/hR(k)/hb + 2*lambda/((hL(k)+hR(k))*r(k));
end
end
