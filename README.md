
# ParaOptus

ParaOptus is a collection of source files dedicated to numerical experiments for the ParaOpt algorithm applied to unstable problems. The project focuses on optimization problems arising in dynamical systems, particularly in regimes where instability or sensitivity significantly affects convergence and performance. The repository is intended to reproduce and support the numerical results presented in the associated scientific work (which will appear in ESIAM). The first version can be found at (HAL)[https://inria.hal.science/hal-04966081v1].

## Scientific Objectives

The primary goal of this repository is to reproduce the numerical experiments presented in specific sections of the reference study.

The following functions generate the figures corresponding to each section:

* fig_section_4_1() — reproduces the numerical experiments of Section 4.1

* fig_section_4_2() — reproduces the numerical experiments of Section 4.2

* fig_section_4_3() — reproduces the numerical experiments of Section 4.3

Each function runs the necessary simulations and post-processing steps to generate the figures associated with its respective section.

## Repository Structure

fig_section_4_*. files: main entry points for reproducing figures.

Auxiliary files: contain supporting routines (model definitions, optimization procedures, numerical solvers, plotting utilities) required by the fig_section_* functions.

All additional files are dependencies used internally by the figure-generation scripts.

## Requirements

To run the experiments, ensure to have MATLAB installed. We have used MATALB R2025b for our simulation.

## Reproducibility

For full reproducibility of the reported results, users should:

Use consistent parameter settings as defined in the figure functions.

Ensure compatible software versions.

Run computations with sufficient numerical precision.

## Usage

To reproduce a specific set of experiments, execute the corresponding function: fig_section_4_1(), fig_section_4_2(), fig_section_4_3(). Figures will be generated according to the configuration specified in each script.

## Purpose

This repository is designed for verification of published numerical results.



