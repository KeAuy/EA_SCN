function yinit = guess(b0, lambda, nu, eta, r)

q    = 1/(eta + nu);
c1   = b0/(2*lambda);
beta = ((b0 - 1)*(nu*b0 + eta) + c1)/(6*lambda);

bNear  = b0 + c1*r + beta*r^2;
dbNear = c1 + 2*beta*r;

x     = 1/r;
bFar  = 1 - q*x + eta*q^3*x^2 + eta*(1 - 2*eta*q)*q^4*x^3;
dbFar = q*x^2 - 2*eta*q^3*x^3 - 3*eta*(1 - 2*eta*q)*q^4*x^4;
dbFar = -dbFar/r^2;

relativeGap = abs(bNear - bFar)/(abs(bNear) + abs(bFar) + 1e-14);

if relativeGap < 0.1 || r < 1
    yinit = [bNear; dbNear];
else
    yinit = [bFar; dbFar];
end

end
