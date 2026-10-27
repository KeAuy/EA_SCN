# Simple buffer case — Figure 1

Numerical and closed-form solutions of the nonlinear boundary value problem

```
b'' + (2/r) b' - (1/λ) (b - 1)(ν b + η) - b/(λ r) = 0,    0 < r < ∞
```

with a regular singular point at `r = 0` and `b(r) → 1` as `r → ∞`. The three scripts reproduce the panels of Figure 1 of the paper and compare the numerical solution with the asymptotic approximations.

## Requirements

MATLAB.

## Usage

Open the folder in MATLAB (or add it to the path) and run

```matlab
fig1a   % λ = 0.1, ν = 0.1, η = 1
fig1b   % λ = 0.1, ν = 10,  η = 1
fig1c   % λ = 1,   ν = 10,  η = 1
```

Each script prints the Newton iterations and opens one figure. The vector graphics `Figure1_A.eps`, `Figure1_B.eps` and `Figure1_C.eps` are the stored outputs of the three panels.

## Files

| File | Description |
| --- | --- |
| `fig1a.m`, `fig1b.m`, `fig1c.m` | Figure scripts: parameters, solve, plot |
| `solve_buffer_bvp.m` | Finite-difference/Newton solver on a uniform grid |
| `buffer_residual_jacobian.m` | Residual and exact tridiagonal Jacobian |
| `buffer_expansion.m` | Small-r (Frobenius) and large-r expansions of `b(r)` |
| `analytic_approximations.m` | Padé, Padé2, Exp-Ser, EXP-Padé, DblExp-Ser, RBA2, linearised solution |
| `buffer_ode.m`, `buffer_bc.m`, `buffer_guess.m` | BVP5C interface: ODE, boundary conditions, initial guess |
| `shooting_fig1_fig3.m` | Reference solution obtained with BVP5C |

## Parameters

| Panel | λ | ν | η | b(r_min) | switch radius | Newton tolerance | max iterations |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1A | 0.1 | 0.1 | 1 | 0.131 | 3.0 | 1e-8 | 125 |
| 1B | 0.1 | 10 | 1 | 0.37 | 1.0 | 1e-11 | 25 |
| 1C | 1 | 10 | 1 | 0.735 | 1.0 | 1e-10 | 25 |

The grid is uniform on `[1e-3, 1e2]` with 50 000 points in all three cases. The near-field value `b(r_min)` is the one used in the paper; the Newton solver imposes it as a Dirichlet condition at the inner boundary.

## Numerical method

1. Second-order central differences on a uniform grid turn the ODE into a nonlinear algebraic system `F(b) = 0`.
2. The system is solved by Newton's method with the exact tridiagonal Jacobian (`buffer_residual_jacobian.m`), using a sparse direct solve.
3. The initial profile switches from the small-r expansion to the large-rexpansion at `switchRadius` (`solve_buffer_bvp.m`).
4. The outer boundary condition is the truncated large-r expansion evaluated at `r_max`. Convergence is reached in 5–6 Newton iterations for all three panels.

## Closed-form approximations

`analytic_approximations.m` returns a struct with these fields:

| Field | Approximation | Used in |
| --- | --- | --- |
| `pade` | [1/1] Padé | 1A, 1B, 1C |
| `pade2` | [2/2] Padé | 1A, 1B, 1C |
| `exp_series` | single exponential series, `S = 1/2` | 1A, 1B, 1C |
| `exp_pade` | exponential-Padé | 1A |
| `double_exp_series` | double exponential series (cubic for `α`) | 1A, 1B, 1C |
| `rba2` | rational approximation RBA2 | 1A |
| `linearized` | solution of the linearised equation | 1C |

Some approximations are complex for the parameters of panel 1A; the scripts plot `real(...)` of those.

## Notes

- `solve_buffer_bvp.m` imposes `b(r_min) = b0` at the inner boundary, while the BVP5C interface (`buffer_bc.m`) imposes the regularity condition `b'(r_min) = b(r_min)/(2λ)`. Both discretise the same singular solution, and the figure scripts use the Dirichlet version.
- The large-r expansions in `buffer_expansion.m` and in the outer boundary condition of `buffer_residual_jacobian.m` are truncated series; their coefficients are kept as in the original implementation.
