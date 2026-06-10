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
