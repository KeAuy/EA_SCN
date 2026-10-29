%MAKEPROGRESS Text progress bar callback for a PARFOR loop over lambda.
%
%   CALLBACK = BUFFER.MAKEPROGRESS(TOTAL, LAMBDAVALUES, TSTART) returns a
%   function handle suitable for AFTEREACH on a PARALLEL.POOL.DATAQUEUE.  Each
%   call increments the completed count and redraws a one-line bar with the
%   current lambda, the elapsed time and the estimated time left.
%
%   The bar is throttled to roughly twelve redraws per second and prints a
%   final newline once all tasks have reported.

function callback = makeProgress(total, lambdaValues, tStart)

state.done = 0;
state.last = 0;
callback   = @update;

    function update(index)
        state.done = state.done + 1;
        elapsed    = toc(tStart);

        if elapsed - state.last < 0.08 && state.done < total
            return
        end
        state.last = elapsed;

        fraction = state.done / total;
        width    = 32;
        filled   = min(width, round(fraction * width));
        bar      = [repmat('=', 1, filled) repmat(' ', 1, width - filled)];
        left     = elapsed / max(fraction, eps) * (1 - fraction);

        fprintf('\r[%s] %5.1f%%  %4d/%d  lambda=%.2e  elapsed %s  left %s   ', ...
            bar, 100 * fraction, state.done, total, lambdaValues(index), ...
            hms(elapsed), hms(left));

        if state.done == total
            fprintf('\n');
        end
    end
end

function text = hms(seconds)
seconds = max(0, round(seconds));
if seconds < 60
    text = sprintf('%ds', seconds);
elseif seconds < 3600
    text = sprintf('%d:%02d', floor(seconds / 60), mod(seconds, 60));
else
    text = sprintf('%d:%02d:%02d', floor(seconds / 3600), ...
        mod(floor(seconds / 60), 60), mod(seconds, 60));
end
end
