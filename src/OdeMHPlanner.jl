"""
    OdeMHPlanner

Framework for uncertainty-aware learning and scenario-based optimal control
of nonlinear dynamical systems with infrequent output measurements.

The package implements a Metropolis–Hastings (MH) sampler with ODE-integrated 
latent trajectories and uses posterior samples to formulate
scenario-based optimal control problems that account for uncertainty in the 
learned dynamics.

Main components:
- Bayesian inference of unknown dynamics and latent states from sparse measurements
- Numerical integration of ODE dynamics within the MH sampler
- Scenario-based optimal control using posterior samples

See the documentation for examples on sampler tuning and optimal control.
"""
module OdeMHPlanner
using DifferentialEquations
using Distributions
using Ipopt
using JuMP
using LinearAlgebra
using Printf
using Random
using Statistics
using StatsBase

export MH_sample, ODE_MH, staged_ODE_MH, compute_ess, compute_gelman_rubin, compute_autocorrelation, solve_MH_OCP

"""
MH sample

Fields:
- `theta`: parameters
- `x_t`: state at current time step ``t=0``
- `x_init`: initial state of the training trajectory, i.e., state at ``t=-T``; this is used to adapt the proposal distribution in the staged sampler
"""
mutable struct MH_sample
    theta::Array{Float64} # parameters
    x_t::Array{Float64} # state at current time step t=0
    x_init::Array{Float64} # initial state of the training trajectory, i.e., state at t=-T; this is used to adapt the proposal distribution in the staged sampler
end

include("sampler.jl")
include("diagnostics.jl")
include("RK4.jl")
include("optimal_control.jl")
end