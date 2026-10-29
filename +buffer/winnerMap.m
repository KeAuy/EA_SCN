%WINNERMAP Index and value of the smallest error along the method dimension.
%
%   [WINNER, BEST] = BUFFER.WINNERMAP(ERRORS) takes an NLAMBDA-by-NNU-by-NMETHOD
%   array and returns the 1-based method index with the smallest error and the
%   corresponding error value at every grid point.

function [winner, best] = winnerMap(errors)

[best, winner] = min(errors, [], 3);
end
