# Reprodue paper

Reproducing paper:  SCN.pdf.

The repository contains 7 independent figure sets. Shared numerics, closed-form approximations, error measures and plotting helpers live in the
MATLAB package [`+buffer`](+buffer/), so the figure scripts only contain parameters, calls and plot configuration.

## Requirements

- MATLAB
- Optimization Toolbox (`fsolve`)

## Quick start

Open this folder in MATLAB, then run one of the figure scripts:

```matlab
figure1A
figure2A
figure3A
figure4A
figure5A
```

A lightweight interface test is available:

```matlab
run('tests/Test.m')
```



## API

All helpers are called with the `buffer.` namespace.

| Function | Purpose |
| --- | --- |
| `buffer.solveUniform` | Newton solver on a uniform grid (Figures 1 and 3) |
| `buffer.solveGeometric` | Newton solver with relaxation on a geometric grid (Figures 2, 4 and 5) |
| `buffer.residualJacobianUniform` | Residual and exact Jacobian for the uniform solver |
| `buffer.residualGeometric`, `buffer.jacobianGeometric` | Residual and Jacobian for the geometric solver |
| `buffer.geometricGrid` | Geometric grid and local spacings |
| `buffer.initialProfile` | Piecewise small-r/far-field starting profile |
| `buffer.expansion`, `buffer.farField` | Small-r and truncated far-field expansions |
| `buffer.approximations` | Evaluate a selected set of closed-form approximations |
| `buffer.pade1`, `buffer.pade2`, `buffer.rba2`, `buffer.expSeries`, `buffer.expPade`, `buffer.doubleExpSeries`, `buffer.doubleExpVariant` | Individual approximations |
| `buffer.chargeDensity`, `buffer.meanError` | Figure 5 charge density and radial error mean |
| `buffer.plotProfiles`, `buffer.plotErrors`, `buffer.palette` | Shared plotting helpers |
| `buffer.generateB0` | Regenerate `b0(nu)` data with BVP5C |
| `buffer.ode`, `buffer.bc`, `buffer.guess` | BVP5C interface used by `generateB0` |



## Code Structure and AI Assistance

The current code framework was generated with the assistance of `DeepSeek-V4.1-flash 1017`. For reference, the original handwritten implementation is available in the `original/` directory.

