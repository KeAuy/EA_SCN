function [aL, aD, aR] = compute_J_Fig6(b, r, hL, hR, lambda, nu, eta)

[N, M] = size(b);
if isscalar(lambda), lambda = lambda*ones(1, M); end

aL = zeros(N, M);
aD = ones(N, M);
aR = zeros(N, M);

k  = (2:N-1).';
hb = 0.5*(hL(k) + hR(k));

aL(k,:) = lambda./hL(k)./hb - 2*lambda./((hL(k) + hR(k)).*r(k));
aR(k,:) = lambda./hR(k)./hb + 2*lambda./((hL(k) + hR(k)).*r(k));
aD(k,:) = (-lambda./hR(k) - lambda./hL(k))./hb ...
        - (2*nu*b(k,:) - nu + eta) - 1./r(k);
end
