function residual = bc(ya, yb, lambda)
%BC Two-point boundary conditions of the buffer equation.
%
%   Regularity at r -> 0:  b'(r_min) = b(r_min)/(2 lambda)
%   Far field:             b(r_max) = 1
%
%   See also BUFFER.ODE, BUFFER.GUESS, BUFFER.GENERATEB0.

residual = [ya(2) - ya(1) / (2 * lambda);
            yb(1) - 1];
end
