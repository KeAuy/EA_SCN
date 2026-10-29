function res = bcfun(ya, yb,lambda)

res1 = ya(2) - ya(1)/(2*lambda);

res2 = yb(1) - 1;

res = [res1; res2];

end


