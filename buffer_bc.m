function residual = buffer_bc(ya, yb, lambda)
%BUFFER_BC Two-point boundary conditions of the buffer equation.
%
%   Residual(1): regularity at r -> 0,
%
%       b'(r_min) = b(r_min) / (2 lambda).
%
%   Residual(2): far-field condition,
%
%       b(r_max) = 1.
%
%   YA and YB are the solution vectors at r_min and r_max, as supplied by
%   BVP5C, see SHOOTING_FIG1_FIG3.
%
%   See also BUFFER_ODE, BUFFER_GUESS.

residual = [ya(2) - ya(1) / (2 * lambda);
            yb(1) - 1];
end
