# Portfolio Replication via Randomized Low-Rank Matrix Approximations

This repository contains the MATLAB codebase and presentation materials for a project focused on Numerical Linear Algebra applied to Quantitative Finance. The core objective is to construct an index-replicating tracking portfolio (S&P 100) using randomized pivoting algorithms for low-rank matrix decomposition (CUR and Interpolative Decomposition - ID).

## Executive Summary
Replicating a financial index or target portfolio with a sparse subset of constituents (e.g., 15 to 30 assets) is a central problem in portfolio management aimed at minimizing transaction costs while matching benchmark performance. 

This project benchmarks randomized numerical linear algebra (RandNLA) algorithms—specifically randomized LUPP, CPQR, DEIM, Leverage Score sampling, and Rangefinder—applied to asset cross-correlation matrices. By extracting column/row skeletons via randomized sketching and power iteration, the pipeline identifies dominant latent factors driving index dynamics and solves a constrained least-squares problem to compute optimal replication weights.

## Key Features & Empirical Findings
* **Randomized Pivoting Algorithms:** Implementation of randomized LUPP (LU with partial pivoting) and CPQR (Column Pivoting QR) with optional power iteration steps to accelerate column/row selection from large-scale covariance/correlation structures.
* **Comparative Algorithm Benchmarking:** Evaluates DEIM (Discrete Empirical Interpolation Method), Leverage Score sampling, and deterministic baselines against randomized pivoting.
* **Index Tracking Performance:** Construct sparse tracking portfolios for the S&P 100 index across multiple rolling horizons (2, 3, and 4 months).
* **Empirical Results:** Demonstrates that randomized LUPP methods consistently outperform Leverage Score sampling in out-of-sample tracking error ($L_1$-norm) while maintaining low computational overhead.

## Quick Start (MATLAB)

### 1. Set Workspace Path
Open MATLAB in the project root directory and add the source code and packages to your path:
```matlab
addpath(genpath('src'));
```

### 2. Core Matrix Approximation & Portfolio Construction
Compute CUR decomposition and calculate index replication weights for a target asset universe:
```matlab
% Select asset subset (e.g., 15 or 30 stocks) via randomized LUPP
[C, R, U, J_s] = CUR_LUPP(A, k);

% Calculate portfolio replication weights
weights = portfolios.calcolaCoefficientiReplica(A, J_s);
```

### 3. Execution & Benchmarking
Run evaluation scripts from the `+plotting` package namespace:
```matlab
% Benchmark sketching algorithms on S&P 100 correlation matrices
plotting.benchmarkSketching();

% Evaluate out-of-sample index tracking performance
plotting.graficaRepliche();
```

## Repository Structure
```text
.
├── data/                                    # Dataset workspace
│   ├── environment.mat                     # Pre-loaded workspace environments
│   └── sp100_daily_2025-01-01_to_2025-09-10.mat # Historical daily price series
├── docs/                                    # Presentation materials
│   └── presentazione_MPALN.pdf              # Beamer presentation slides
└── src/                                     # Core MATLAB codebase
    ├── CUR_CPQR.m                           # CUR decomposition via Column Pivoting QR
    ├── CUR_DEIM.m                           # CUR decomposition via DEIM algorithm
    ├── CUR_LeverageScore.m                  # Leverage score sampling CUR
    ├── CUR_LUPP.m                          # CUR decomposition via LUPP
    ├── CUR_U_eval.m                         # Intersection matrix U evaluator
    ├── rangefinder.m                        # Randomized rangefinder algorithm
    ├── srcur.m                              # Subsampled randomized CUR routine
    ├── +plotting/                           # Benchmarking and portfolio visualization
    │   ├── benchmarkSketching.m
    │   ├── graficaRepliche.m
    │   ├── plot_sp100_range.m
    │   ├── test_CUR.m
    │   ├── test_CUR_large.m
    │   └── test_CUR_large1.m
    ├── +utils/                              # Helper functions, RSVD, and sketch selectors
    │   ├── cur_algos.m
    │   ├── embed.m
    │   ├── RSVD.m
    │   ├── scegliAlgoritmo.m
    │   ├── scegliSketch.m
    │   ├── SNN.m
    │   └── SNN1.m
    ├── chebfun/                             # Open-source Chebfun library dependency
    └── portfolios/                          # Portfolio replication routines
        ├── calcolaCoefficientiReplica.m
        └── calcolaCorrelazioneSto.m
```

## Requirements
* **MATLAB** (R2021a or newer recommended)
  * *Statistics and Machine Learning Toolbox*
  * *Optimization Toolbox*