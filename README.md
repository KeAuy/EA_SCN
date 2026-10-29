# Paper Reproduction

Reproduction of the figures of SCN.pdf. 

The repository contains seven independent figure sets; every panel is a standalone script that solves the buffer equation:

$$
\lambda \Delta b = (b - 1)(\nu b + \eta) + \frac{b}{\mathbf{r}}
$$

We use numerical method to test the efficiency of the approximation. All output lands in `figures/`.
## Requirements

- MATLAB (developed with R2022b)
- Optimization Toolbox 
- Parallel Computing Toolbox (optional; the `b0` grid of Figures 6/7 is built with `parfor` when it is available and falls back to a serial loop otherwise)

## Quick start

```matlab
run('tests/Test.m')        % lightweight interface and solver test
run('Figure1/figure1A.m')  % single panel, writes figures/figure1A.png
```

## Notes

-  The current code framework was generated with the assistance of `DeepSeek-V4.1-flash 1017`. For reference, the original handwritten implementation is available in the `original/` directory.
