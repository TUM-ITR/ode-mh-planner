# API Reference

This page documents the main functions and types provided by `OdeMHPlanner`.

## Sampling

Functions and types related to Bayesian learning of system dynamics and latent state trajectories using the Metropolis–Hastings (MH) sampler.

```@docs
MH_sample
ODE_MH
staged_ODE_MH
```

## Analysis and Diagnostics

Utilities for analyzing and diagnosing the MCMC chains produced by the MH sampler.

```@docs
compute_autocorrelation
compute_ess
compute_gelman_rubin
```

## Optimal Control

Functions for formulating and solving the scenario-based optimal control problem using posterior samples obtained from MH.

```@docs
solve_MH_OCP
```
