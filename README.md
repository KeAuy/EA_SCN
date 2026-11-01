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

## Notes

-  The current code framework was generated with assistance from `DeepSeek-V4.1-flash 1017`. For reference, the original handwritten implementation is available in the `original`.
