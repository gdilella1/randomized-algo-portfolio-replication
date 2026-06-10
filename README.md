# Portfolio Replication via Randomized Low-Rank Matrix Approximations

This repository contains the MATLAB codebase and the Beamer presentation for a project focused on Numerical Linear Algebra applied to Quantitative Finance. The core objective is to construct an index-replicating portfolio using randomized probabilistic algorithms for low-rank matrix decomposition.

## Executive Summary
Replicating a financial index or a target portfolio with a limited subset of assets is a classic problem in portfolio management. This project leverages **Randomized Numerical Linear Algebra (RandNLA)** to accelerate and regularize the selection and weighting of assets. By computing a probabilistic low-rank approximation (such as Randomized SVD or Interpolative Decomposition) of the asset returns covariance matrix, we filter out market noise and identify the dominant latent factors driving the index performance.

## Key Features
* **Probabilistic Low-Rank Approximation:** MATLAB implementation of randomized algorithms to compute fast, low-rank approximations of large-scale financial matrices.
* **Portfolio Replication:** Optimization pipeline that utilizes the low-rank structures to construct tracking portfolios with minimized tracking error and controlled turnover.
* **Beamer Presentation:** The complete theoretical framework, error bounds analysis, and empirical results are detailed in the included LaTeX Beamer slides.

## Repository Structure
* `src/`: Core MATLAB functions and the main execution script (`main.m`).
  * `src/algorithms/`: Randomized matrix decomposition routines.
  * `src/portfolio/`: Replication and optimization scripts.
  * `src/plotting/`: Scripts for tracking error and asset weight visualization.
* `docs/`: The Beamer presentation PDF and LaTeX source files.

## Theoretical Background
Traditional methods like standard Principal Component Analysis (PCA) via deterministic SVD can be computationally expensive on large-scale high-frequency financial data and prone to overfitting under market regime shifts. Randomized algorithms solve this by drawing random samples from the matrix column space, offering:
1. Significant computational speedups.
2. Robustness against high-frequency noise.
3. Strict, predictable probabilistic error bounds.

## How to Run
1. Open MATLAB and navigate to the root directory.
2. Add the source folders to your path: `addpath(genpath('src/'));`
3. Execute `main.m` to run the simulation, construct the replicating portfolio, and output the tracking performance plots.
