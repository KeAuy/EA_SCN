function report_Serial_vs_12workers(result)
    speedup       = result.serialMedianSec / result.parallelMedianSec;
    efficiency    = speedup / result.workers;
    endToEnd      = result.parallelMedianSec + result.poolStartSec;
    speedupEnd    = result.serialMedianSec / endToEnd;
    efficiencyEnd = speedupEnd / result.workers;

    fprintf('\n================ serial vs %d workers ================\n', ...
        result.workers);
    fprintf('serial median        : %.3f s\n', result.serialMedianSec);
    fprintf('parallel median      : %.3f s\n', result.parallelMedianSec);
    fprintf('pool startup         : %.3f s\n', result.poolStartSec);
    fprintf('speedup (solve only) : %.3f x\n', speedup);
    fprintf('efficiency           : %.1f %%\n', 100 * efficiency);
    fprintf('end-to-end time      : %.3f s\n', endToEnd);
    fprintf('speedup (end-to-end) : %.3f x\n', speedupEnd);
    fprintf('efficiency end-to-end: %.1f %%\n', 100 * efficiencyEnd);
end
