%ODE Right-hand side of the buffer equation in first-order form, y = [b; db].
%
%       d b/dr  = db
%       d db/dr = -(2/r) db + (b - 1)(nu b + eta)/lambda + b/(lambda r).

function dydr = ode(r, y, lambda, nu, eta)

b  = y(1);
db = y(2);

dydr = [db;
        -(2 / r) * db + (b - 1) * (nu * b + eta) / lambda + b / (lambda * r)];
end
