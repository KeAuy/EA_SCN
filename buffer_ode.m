function dydr = buffer_ode(r, y, lambda, nu, eta)
%BUFFER_ODE Right-hand side of the buffer equation in first-order form.
%
%   The second-order equation
%
%       b'' + (2/r) b' - (1/lambda) (b - 1) (nu b + eta) - b / (lambda r) = 0
%
%   is written as the first-order system y = [b; db/dr]:
%
%       d b /dr = db
%       d db/dr = -(2/r) db + (b - 1) (nu b + eta) / lambda + b / (lambda r)
%
%   This is the ODE function used by BVP5C, see SHOOTING_FIG1_FIG3.
%
%   See also BUFFER_BC, BUFFER_GUESS.

b  = y(1);
db = y(2);

dydr = [db;
        -(2 / r) * db + (b - 1) * (nu * b + eta) / lambda + b / (lambda * r)];
end
