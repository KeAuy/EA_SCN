# Numerical and analytical approximations for solving an ODE
- rewritting paper: SCN.pdf

MATLAB code for the numerical solution of the nonlinear boundary-value problem

    b'' + (2/r) b' - (1/λ)(b - 1)(ν b + η) - b/(λ r) = 0,    0 < r < ∞,
    b(r) → 1  as  r → ∞,

and for the closed-form approximations studied in the paper.

## Figures

| # | Folder | What it shows | Scripts |
|---|--------|---------------|---------|
| 1 | `Figure1/` | b(r) profiles and analytical approximations, 3 parameter sets | `fig1a.m`, `fig1b.m`, `fig1c.m` |
| 2 | `Figure2/` | log10 of the r-averaged error vs ν, λ = 0.02 / 2 / 20 | `Simple_*_Fig2_A/B/C.m` |
| 3 | `Figure3/` | b(r) vs approximations for the Figure-1 cases | `simple_buffer_case_Fig3_A/B/C.m` |
| 4 | `Figure4/` | error comparison incl. ExpVar / ExpGlobal / DblExpVar / DblExpGlobal | `simple_buffer_case_Fig4_A/B/C.m` |
| 5 | `Figure5/` | same comparison for the concentration c(r) | `simple_buffer_case_Fig5_A/B/C.m` |
| 6 | `Figure6/` | best-approximation maps over (λ, ν) | `simple_buffer_case_Fig6_*.m` |
