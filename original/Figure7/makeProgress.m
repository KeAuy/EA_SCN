function cb = makeProgress(total, lambda_v, t0)
state.done  = 0;
state.tlast = 0;
cb = @update;

    function update(i)
        state.done = state.done + 1;
        el = toc(t0);
        if el - state.tlast < 0.08 && state.done < total
            return
        end
        state.tlast = el;

        frac  = state.done / total;
        nb    = 32;
        nf    = min(nb, round(frac*nb));
        bar   = [repmat('=',1,nf) repmat(' ',1,nb-nf)];
        eta_s = el / max(frac, eps) * (1 - frac);

        fprintf('\r[%s] %5.1f%%  %4d/%d  lam=%.2e  已用 %s  剩 %s   ', ...
            bar, 100*frac, state.done, total, lambda_v(i), hms(el), hms(eta_s));
        if state.done == total
            fprintf('\n');
        end
    end
end

function s = hms(t)
t = max(0, round(t));
if t < 60
    s = sprintf('%ds', t);
elseif t < 3600
    s = sprintf('%d:%02d', floor(t/60), mod(t,60));
else
    s = sprintf('%d:%02d:%02d', floor(t/3600), mod(floor(t/60),60), mod(t,60));
end
end