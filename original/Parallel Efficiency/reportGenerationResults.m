function reportGenerationResults(results)
    T1 = results(1).runtimeSec;

    fprintf('\n================ b0 generation 100 x 100 ================\n');
    fprintf('%8s %12s %10s %12s %12s %12s\n', ...
        'Workers', 'Runtime(s)', 'Speedup', 'Efficiency', 'Pool(s)', 'ActualWorker');

    for k = 1:numel(results)
        speedup = T1 / results(k).runtimeSec;
        efficiency = speedup / results(k).workers;

        fprintf('%8d %12.3f %10.3f %12.3f %12.2f %12d \n', ...
            results(k).workers, results(k).runtimeSec, speedup, ...
            efficiency, results(k).poolStartSec, ...
            results(k).actualWorkers);
    end
end
