# OdeMHPlanner.jl

Welcome to the documentation for **OdeMHPlanner.jl**, a Julia package for uncertainty-aware learning and planning in dynamical systems with unknown parameters and infrequent output measurements.

OdeMHPlanner implements the method described in:
> **Learning Dynamics from Infrequent Output Measurements for Uncertainty-Aware Optimal Control**
> *Robert Lefringhausen, Theodor Springer, Sandra Hirche*
> arXiv:2512.08013 (2025)
> <https://arxiv.org/abs/2512.08013>

## Overview

The package targets control problems in which the system dynamics are partially unknown and the state is only indirectly observed through infrequent and noisy output measurements. Rather than identifying a single nominal model, OdeMHPlanner explicitly represents uncertainty over both the model parameters and the latent state trajectory, and propagates this uncertainty into the control design.

The approach follows a Bayesian workflow:

1. **Learning:** A Metropolis–Hastings (MH) sampler, equipped with a numerical ODE solver, draws samples from the posterior distribution over unknown parameters and latent state trajectories, conditioned on infrequent input–output measurements.
2. **Planning:** The posterior samples are used to formulate a scenario-based optimal control problem (OCP), yielding control inputs that explicitly account for the inferred uncertainty.

This enables principled uncertainty quantification and safer decision-making compared to point-estimate-based approaches.

## Installation

This package is not registered in the General registry. Clone the repository and instantiate the environment locally:

```bash
git clone https://github.com/TUM-ITR/ode-mh-planner.git
cd ode-mh-planner
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Getting Started

The [Experiments](@ref experiments) section provides a structured entry point into the package. The framework is demonstrated on a glucose regulation task for Type 1 diabetes patients, and the section is organized as follows:

- **[Experiments overview](@ref experiments):** Simulation setup, model description, prior distributions, cost functional, and baseline definitions.
- **[Inference and Sampler Tuning](@ref sampling):** How to infer unknown dynamics and latent state trajectories from infrequent input–output data using the MH sampler, including tuning and diagnostics for reliable posterior exploration.
- **[Optimal Control](@ref optimal-control):** How to formulate and solve a scenario-based OCP using the inferred posterior samples.
- **[Monte Carlo Study](@ref monte-carlo):** Statistical evaluation of the method across 100 independent runs.

Each subpage combines a scientific explanation with practical code examples and guidance.