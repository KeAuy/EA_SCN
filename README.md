# Paper Reproduction

Reproduction paper: **SCN.pdf**. 

The repository contains seven figure sets. We use the paper's methods to reproduce the figures in the paper, which provide numerical and analytic approximate methods for solving the upper semilinear elliptic partial differential equation:

$$
\lambda \Delta b = (b - 1)(\nu b + \eta) + \frac{b}{\mathbf{r}}
$$

All output lands in `figures/`.

## Requirements

- MATLAB (developed with R2022b)
- Optimization Toolbox 
- Parallel Computing Toolbox (optional; the grid of Figures 6/7 is built with `parfor`)

## Quick start

```matlab
run('tests/Test.m')      
run('Figure1/figure1A.m') 
```

## Parallel Computing

### Serial vs. 12-workers

| Metric                               | Result    |
| ------------------------------------ | --------- |
| Serial median runtime                | 113.489 s |
| Parallel median runtime (solve only) | 21.696 s  |
| Parallel pool startup                | 17.943 s  |
| Speedup (solve only)                 | 5.231×    |
| Parallel efficiency (solve only)     | 43.6%     |
| End-to-end time (startup + solve)    | 39.638 s  |
| Speedup (end-to-end)                 | 2.863×    |
| Efficiency (end-to-end)              | 23.9%     |


### Worker-Count Scalability

|Workers|Runtime (s)|Speedup|Efficiency|Pool Startup (s)|Actual Workers|
|---|---|---|---|---|---|
|1|110.731|1.000×|100.0%|0.00|1|
|2|61.032|1.814×|90.7%|13.29|2|
|4|36.697|3.017×|75.4%|9.21|4|
|8|24.500|4.520×|56.5%|12.66|8|
|12|21.432|5.167×|43.1%|17.04|12|

Source code and data can be found at this path:`original/Parallel Efficiency`

## Notes

-  The current code framework was generated with assistance from `DeepSeek-V4.1-flash 1017`. For reference, the original handwritten implementation is available in the `original`.
