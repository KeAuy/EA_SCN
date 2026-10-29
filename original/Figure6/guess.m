function yinit = guess(b0, lambda, nu, eta, r)

    q1   = 1/(eta + nu);
    beta = ((b0-1)*(nu*b0+eta) + b0/(2*lambda)) / (6*lambda);

    b1  = b0 + (b0/(2*lambda))*r + beta*r^2;
    Db1 = b0/(2*lambda) + 2*beta*r;

    x  = 1/r;
    b2  = 1 - q1*x + eta*q1^3*x^2 ...
            + eta*(1-2*q1*eta)*q1^4*x^3;
    Db2 = q1*x^2 - 2*eta*q1^3*x^3 ...
            - 3*eta*(1-2*q1*eta)*q1^4*x^4;

    Db2 = -Db2 / r^2;

    if abs(b1-b2)/(abs(b1)+abs(b2)+1e-14) < 0.1

        yinit = [b1; Db1];
    else

        if r < 1
            yinit = [b1; Db1];
        else
            yinit = [b2; Db2];
        end
    end
end
