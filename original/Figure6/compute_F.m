function F = compute_F_Fig6(b0, b, r, hL, hR, lambda, nu, eta)

[N, M] = size(b);
if isscalar(lambda), lambda = lambda*ones(1, M); end
if isscalar(b0),     b0     = b0*ones(1, M);     end

q = 1/(eta + nu);

F = zeros(N, M);
F(1,:) = b(1,:) - b0;

k   = (2:N-1).';
hb  = 0.5*(hL(k) + hR(k));

F(k,:) = lambda .* ((b(k+1,:) - b(k,:))./hR(k) - (b(k,:) - b(k-1,:))./hL(k))./hb ...
       + 2*lambda .* (b(k+1,:) - b(k-1,:))./((hL(k) + hR(k)).*r(k)) ...
       - (b(k,:) - 1).*(nu*b(k,:) + eta) - b(k,:)./r(k);

F(N,:) = b(N,:) - (1 - q/r(N) + eta*q^3/r(N)^2);
end
