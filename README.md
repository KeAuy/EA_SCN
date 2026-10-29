# Paper Reproduction

Reproduction of the figures of SCN.pdf. The repository contains seven
independent figure sets; every panel is a standalone script that solves the
buffer equation itself and writes its PNG into `figures/`.

## Requirements

- MATLAB (developed with R2022b)
- Optimization Toolbox (`fsolve` is used when refining the Padé2 coefficients)
- Parallel Computing Toolbox (optional; the `b0` grid of Figures 6/7 is built
  with `parfor` when it is available and falls back to a serial loop otherwise)

## Quick start

```matlab
run('tests/Test.m')        % lightweight interface and solver test
run('Figure1/figure1A.m')  % single panel, writes figures/figure1A.png
```

## Figures

| Folder | Panels | Content |
| --- | --- | --- |
| `Figure1` | A, B, C | b(r) from the uniform-grid solver vs. the approximations |
| `Figure2` | A, B, C | r-averaged b error vs. nu |
| `Figure3` | A, B, C | Exp / DblExp approximations vs. the numerical solution |
| `Figure4` | A, B, C | b-error comparison with the minimum marked |
| `Figure5` | A, B, C | charge-density c = nu (b - 1) + 1/r error |
| `Figure6` | A1/A2, B1/B2 | best-method map and error map on a 100x100 grid |
| `Figure7` | A, B | rule-based method regions and the error maps they imply |

All output lands in `figures/`.

## Data

Nothing is stored on disk: no `.mat`, no cache, no generated data of any kind.
The inner boundary value `b0 = b(r_min)` is computed in memory at run time with
BVP5C and handed straight to the solver:

- `buffer.generateB0` for one lambda and a vector of nu (Figures 2, 4, 5);
- `buffer.generateB0Grid` for a lambda-by-nu grid (Figures 6, 7).

So every figure is reproducible from the code alone, and the repository stays
small enough for GitHub. The trade-off is run time: the 100-by-100 `b0` grid is
rebuilt on every run, which takes roughly half a minute on a parallel pool and
a few minutes serially.

## Layout

```
+buffer/      solver, approximations, plotting helpers
Figure1-7/    one script per panel
original/     pre-refactor implementations, kept for reference
tests/Test.m  smoke test over the public interface
figures/      generated PNGs
```

## Notes

- Figures 6 and 7 use 100x100 parameter grids. The `b0` rows are independent per
  lambda, so they are built in a `parfor` loop (one row per worker, progress bar
  from a DataQueue); the Newton solve itself then batches all lambda values of
  one nu and is quick.
- The legacy `original/` scripts are reference material: they are kept
  comment-free and are not part of the automated run. `original/Figure6` and
  `original/Figure7` hold the pre-refactor implementation; there
  `shooting_fig6.m` rebuilds `b0_Fig6.mat` with a PARFOR loop (one lambda row
  per worker, progress bar from a DataQueue) and the panel scripts read that
  file. `*.mat` is git-ignored, so the generated file never enters the
  repository.
