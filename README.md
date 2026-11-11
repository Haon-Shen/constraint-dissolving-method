# A Constraint Dissolving Method for Nonlinearly Constrained Optimization over Nonconvex Sets: MATLAB code

This repository contains the MATLAB implementation of the constraint dissolving method **CD** and the scripts that reproduce all tables and figures of the numerical section of the paper

> *A Constraint Dissolving Method for Nonlinearly Constrained Optimization over Nonconvex Sets* (submitted to Mathematical Programming Computation).

The methods solve problems of the form

    min f(x)   s.t.   G(x) ∈ C,   x ∈ D,

where C is closed and convex and D is closed and possibly nonconvex; the methods access D only through a projection onto it. CD is compared with

* **ALM**: the safeguarded augmented Lagrangian method of Jia, Kanzow, Mehlitz, Wachsmuth, *Math. Program.* 199 (2023) 1365–1415, Algorithm 4.1, and
* **PD**: the inexact penalty decomposition method of Kanzow, Lapucci, *Comput. Optim. Appl.* 85 (2023) 937–971, Algorithms 1 and 2.

## Requirements

MATLAB, without any toolbox. The results of the paper were computed with MATLAB R2025b on a MacBook Air (Apple M4, 16 GB). The code uses `maxk`, `tiledlayout` and `exportgraphics`, so MATLAB R2020a or later is needed.

## Quick start

In MATLAB, from the root folder of the repository,

```matlab
run_all
```

reproduces all experiments except PD on the MaxCut problems and, unless the stock return data have been added (see [Data](#data)), the experiment of Section 5.2.2 (about 11 minutes in total; see `run_all.m` for the running times of the individual experiments and for how to include PD on MaxCut). Every experiment can also be run on its own from its folder, e.g.

```matlab
cd MPCC
run_mpcc        % CD and ALM
run_mpcc_pd     % PD (reads the starting points written by run_mpcc)
```

The drivers add `../common` to the MATLAB path themselves.

Each experiment writes its output to the subfolder `results/` of its folder. The outputs reported in the paper, together with the underlying per-instance data (`.mat`), are in the subfolders `reference/`. The tables and figures can be redrawn from these data without running the solvers, e.g. `maxcut_table_paper('reference')` in `MaxCut/`. The iteration numbers, objective values and solutions are deterministic; the running times depend on the machine.

## Tables and figures of the paper

| Paper | Experiment | Driver (folder) | Output in `results/` |
|---|---|---|---|
| Table 1 | Section 5.1, complementarity-constrained Rosenbrock problem | `run_mpcc`, `run_mpcc_pd` (`MPCC/`) | `mpcc_rosenbrock.tex` (rows of CD and ALM), `mpcc_rosenbrock_pd.tex` (row of PD), `*.txt` (summaries) |
| Table 2, Figure 1 | Section 5.2.1, interval-budget portfolio model, orl and pard instances | `run_portfolio` (`Portfolio/`), which calls `portfolio_table_paper` and `plot_portfolio_paper` | `portfolio_table.tex`, `portfolio_table.txt`, `fig_portfolio.pdf`; per-instance results in `orl/` and `pard/` |
| Tables 3 and 4, Figures 2 and 3 | Section 5.2.2, mean–variance model on market data | `run_mv` (`PortfolioMV/`), which calls `mv_write_output` and `plot_mv_paper` | `mv_tables.tex`, `mv_summary.txt`, `fig_mv_out.pdf` (Figure 2), `fig_mv_regimes.pdf` (Figure 3) |
| Table 5, Figures 4 and 5 | Section 5.3, MaxCut problems, rudy collection | `run_maxcut` (`MaxCut/`), followed by `maxcut_table_paper`, `plot_maxcut_paper`, `plot_maxcut_profiles` | `maxcut_table.tex`, `maxcut_table.txt`, `fig_maxcut.pdf` (Figure 4), `fig_maxcut_profiles.pdf` (Figure 5) |

The `.tex` files contain the table rows as they appear in the paper; the `.txt` files contain the same numbers and the further counts quoted in the text.

## Parameters

* All parameters shared by the experiments are set in `common/default_params.m`.
* ALM uses the parameters of Jia et al. (2023, Section 6): growth factor 10, decrease factor 0.8, and the initial penalty parameter given by the formula in `common/initial_penalty.m`. PD uses the growth factor 1.1 and the same initial penalty parameter, as well as the same inner tolerance as CD and ALM.
* The parameters of CD (τ, η, α, the scaling c of Q_D and, except for the Rosenbrock problem, a fixed initial penalty parameter β₀) are set for each problem class at the beginning of the corresponding driver.
* All methods share the starting points, the upper bound 1e8 on the penalty parameter, the budget of 1e7 oracle evaluations and the inner tolerance 1e-4/√(k+1).

## Repository layout

```
run_all.m                    reproduces all experiments
common/                      solvers
  cdm.m                      CD, the constraint dissolving method of the paper
  cdm_penalty.m              constraint dissolving penalty function and its gradient
  safeguarded_alm.m          ALM [Jia et al. 2023, Algorithm 4.1]
  pd.m                       PD [Kanzow, Lapucci 2023, Algorithms 1 and 2]
  spectral_gradient.m        subproblem solver of CD and ALM [Jia et al. 2023, Algorithm 3.1]
  default_params.m, initial_penalty.m, proj_box.m, proj_sparse_nonneg.m
MPCC/                        Section 5.1
  run_mpcc.m, run_mpcc_pd.m  drivers
  mpcc_rosenbrock_problem.m  problem data (f, G, C, D, Q_D, Q_C)
  mpcc_starting_points.m, mpcc_run_batch.m, proj_complementarity.m
Portfolio/                   Section 5.2.1
  run_portfolio.m            driver
  portfolio_problem.m        problem data
  portfolio_table_paper.m, plot_portfolio_paper.m
  portfolio_load_instances.m, read_portfolio_instance.m, portfolio_warm_start.m,
  portfolio_run_batch.m, portfolio_summary.m, proj_budget.m
  data/                      instances of Frangioni and Gentile (see below)
PortfolioMV/                 Section 5.2.2
  run_mv.m                   driver
  mv_problem.m               problem data
  mv_rolling.m               rolling-window experiment for one data set
  mv_metrics.m, mv_write_output.m, plot_mv_paper.m
  mv_load_datasets.m, mv_dataset_info.m, mv_warm_start.m
  data/                      place for the stock return data (not distributed, see below)
MaxCut/                      Section 5.3
  run_maxcut.m               driver
  maxcut_problem.m           problem data
  maxcut_load.m, maxcut_read.m, maxcut_sdp_start.m, maxcut_cut.m, proj_rank1_psd.m
  maxcut_table_paper.m, plot_maxcut_paper.m, plot_maxcut_profiles.m
  data/                      rudy graphs of the Biq Mac library (see below)
```

## Using CD for other problems

A problem is described by a struct `prob`; its fields are listed in the headers of `common/cdm_penalty.m` (f, G, their derivatives, and the projective mappings Q_D and Q_C) and `common/cdm.m` (the projections onto D and C). The files `*_problem.m` are examples. Then

```matlab
addpath common
par = default_params();
[x, s, out] = cdm(prob, x0, par);
```

## Data

* `Portfolio/data/`: the orl and pard instances (n = 200, 300, 400) of the collection of A. Frangioni and C. Gentile, *Oper. Res. Lett.* 35 (2007) 181–185, in their original format (see `Portfolio/data/README.txt`). Only the covariance matrices and the expected returns are used.
* `PortfolioMV/data/`: the stock return data of Section 5.2.2 are **not distributed with this repository**; they are available from the corresponding author on reasonable request. The experiment uses daily returns of five data sets: `init_stock_data.mat` contains NASDAQ100 and FTSE100 as used by Wu, Sun, Ge, Allen-Zhao, Zeng, *Eur. J. Oper. Res.* 319 (2024) 820–833; `SnPday_Bear.mat`, `SnPday_COVID_19.mat` and `SnPday_Bull.mat` contain the S&P500 data sets of a declining market, the COVID-19 period and a rising market. After the four files have been placed in `PortfolioMV/data/`, `run_mv` reproduces Tables 3 and 4 and Figures 2 and 3; without them, `run_all` skips this experiment. The results reported in the paper are in `PortfolioMV/reference/` in any case. The time periods are listed in `PortfolioMV/mv_dataset_info.m`. Stocks without a price at the beginning of a period are removed by `mv_load_datasets.m`.
* `MaxCut/data/`: the 130 graphs of the rudy collection of the Biq Mac library (http://biqmac.aau.at/biqmaclib.html) and their optimal cut values (`rudy_opt.txt`).

The data in the folders `Portfolio/data/` and `MaxCut/data/` come from the sources listed above and remain subject to the terms of those sources.
