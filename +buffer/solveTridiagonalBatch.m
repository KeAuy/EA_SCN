%SOLVETRIDIAGONALBATCH Thomas algorithm for many independent tridiagonal systems.
%
%   LOWER, MAIN, UPPER and RHS are N-by-NCOLUMN; LOWER(1) and UPPER(N) are
%   unused and may be zero.  The elimination is vectorized over the columns.

function x = solveTridiagonalBatch(lower, main, upper, rhs)
arguments
    lower double
    main double
    upper double
    rhs double
end

[n, nColumn] = size(main);

if ~isequal(size(lower), [n, nColumn]) ...
        || ~isequal(size(upper), [n, nColumn]) ...
        || ~isequal(size(rhs), [n, nColumn])
    error('buffer:sizeMismatch', ...
        'LOWER, MAIN, UPPER and RHS must have the same size.');
end

diagonal = zeros(n, nColumn);
solution = zeros(n, nColumn);

diagonal(1, :) = upper(1, :) ./ main(1, :);
solution(1, :) = rhs(1, :) ./ main(1, :);

for i = 2:n
    denominator = main(i, :) - lower(i, :) .* diagonal(i-1, :);
    if i < n
        diagonal(i, :) = upper(i, :) ./ denominator;
    end
    solution(i, :) = (rhs(i, :) - lower(i, :) .* solution(i-1, :)) ...
        ./ denominator;
end

x = solution;
for i = n-1:-1:1
    x(i, :) = solution(i, :) - diagonal(i, :) .* x(i+1, :);
end
end
