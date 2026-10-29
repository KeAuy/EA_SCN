%MEANERRORGRID Radial mean of the error on a lambda-nu grid.
%
%   ERR = BUFFER.MEANERRORGRID(REFERENCE, APPROX, MODE, OPTIONS) averages along
%   the first dimension and returns NLAMBDA-by-NNU.
%
%   OPTIONS.ValidMask marks the entries allowed in the mean; entries outside
%   the mask take OPTIONS.InvalidValue (default Inf) before averaging, which
%   reproduces the "invalid method" convention of Figure 6B.

function err = meanErrorGrid(reference, approx, mode, options)

arguments
    reference double
    approx double
    mode (1,:) char {mustBeMember(mode, {'abs', 'log'})} = 'abs'
    options.ValidMask = []
    options.InvalidValue (1,1) double = Inf
end

if ~isequal(size(reference), size(approx))
    error('buffer:sizeMismatch', ...
        'REFERENCE and APPROX must have the same size.');
end

if strcmpi(mode, 'log')
    pointwise = abs(log(approx) - log(reference));
else
    pointwise = abs(approx - reference);
end

if ~isempty(options.ValidMask)
    if ~isequal(size(options.ValidMask), size(pointwise))
        error('buffer:sizeMismatch', ...
            'VALIDMASK must have the same size as the error array.');
    end
    pointwise(~options.ValidMask) = options.InvalidValue;
end

averaged = mean(pointwise, 1);
err = reshape(averaged, size(pointwise, 2), size(pointwise, 3));
end
