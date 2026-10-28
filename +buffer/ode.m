function dydr = ode(r, y, lambda, nu, eta)
%ODE Right-hand side of the buffer equation in first-order form.
%
%   The second-order equation is written as y = [b; db/dr]:
%
%       d b /dr = db
%       d db/dr = -(2/r) db + (b - 1)(nu b + eta)/lambda + b/(lambda r)
%
%   See also BUFFER.BC, BUFFER.GUESS, BUFFER.GENERATEB0.

b  = y(1);
db = y(2);

dydr = [db;
        -(2 / r) * db + (b - 1) * (nu * b + eta) / lambda + b / (lambda * r)];
end
