# [Inference and Sampler Tuning](@id sampling)

This page walks through the Bayesian inference step of the framework: running the MH sampler to obtain posterior samples over unknown parameters and latent state trajectories, and tuning the sampler for reliable posterior exploration. The complete script is available at `experiments/sampler_tuning.jl`.

For the model definition, measurement setup, and prior distributions used in this example, see the [Experiments overview](@ref experiments).

## Data Generation

The script generates synthetic data from a fixed ground-truth parameter realization and initial state:

```julia
theta_true = [
    log(0.015),     # log(p2)
    log(2e-6),      # log(p3)
    log(0.21)       # log(n)
]

x_init_true = [78.0, 5e-4, 8.0]   # [G(-T), X(-T), I(-T)]
```

The ground-truth trajectory is simulated over the full 18-hour window using the Tsitouras 5/4 adaptive solver. A total of ``M = 300`` measurement instants are sampled uniformly at random: ``M_{\text{train}} = 200`` from the training window (6 am–6 pm) and ``M_{\text{test}} = 100`` from the test window (6 pm–12 am). The glucose output at each measurement time is corrupted by additive Gaussian noise with standard deviation ``\sigma = 8`` mg/dL.

```julia
M = 300
M_train = 200
M_test = 100

t_train = sort(-T_train .+ rand(M_train) .* T_train)
t_test = sort(rand(M_test) .* T_test)

# Simulate ground-truth trajectory and generate noisy measurements.
prob = ODEProblem(ode_rhs, x_init_true, t_span, theta_true)
sol = solve(prob, Tsit5(); saveat=t_m)
x = Array(sol)
y = zeros(n_y, M)
for m in 1:M
    y[:, m] = g_theta(theta_true, x[:, m], u_t_bolus(t_m[m]), t_m[m]) + sample_w_theta(theta_true, 1)
end
```

The training measurements are used for inference; the test measurements are used afterwards to evaluate how well the inferred posterior models predict future system behavior.

## ODE Solver Configuration

Each iteration of the MH sampler requires integrating the ODE forward to evaluate the likelihood. The ODE integration is handled via [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/), which allows using any of its available solvers.

```julia
rk4_step_size = 0.5       # step size in minutes
ODE_solver = RK4()
ODE_solver_opts = (dt=rk4_step_size, adaptive=false)
```

Here we use the fixed-step RK4 solver with a step size of 0.5 min to ensure consistency with the discretization used in the optimal control problem downstream. For inference-only tasks or longer training windows, a variable-step solver with appropriate tolerances can reduce runtime significantly.

## Prior Specification

The prior encodes domain knowledge about physiologically plausible parameter ranges and initial conditions. Since we work with parameters on a log scale to ensure positivity, the Gaussian priors on ``\theta = (\log p_2, \log p_3, \log n)`` correspond to lognormal priors on the original parameters:

```julia
const theta_mean = [
    -4.26,      # log(p2)
    -13.27,     # log(p3)
    -1.66       # log(n)
]

const theta_var = [
    0.18^2,
    0.28^2,
    0.23^2
]
```

The initial state is modeled as Gaussian, centered near a basal fasting equilibrium:

```julia
const x_init_mean = [80.0, 0.0, 7.0]         # [G(-T), X(-T), I(-T)]
const x_init_var  = [10.0^2, 0.001^2, 2.0^2]
```

See the [Experiments overview](@ref experiments) for the rationale behind these choices.

## Running the Staged MH Sampler

The sampler uses a staged strategy to build up an effective proposal distribution before drawing the final posterior samples. The key tuning parameters are:

```julia
K       = Int(1e5)  # number of MH samples in final stage
k_d     = 0         # thinning factor (0 = no thinning, for diagnostic purposes)
K_b     = 200       # burn-in length per stage
M_chunk = 5         # measurements added per stage
K_stage = 500       # samples drawn per stage
alpha   = 0.85      # proposal covariance scaling
```

At each stage, the sampler incorporates ``M_{\text{chunk}} = 5`` additional measurements into the likelihood and draws ``K_{\text{stage}} = 500`` samples after a burn-in of ``K_b = 200``. The proposal covariance is initialized from the prior and then updated at each stage using the empirical covariance of the chain, scaled by ``\alpha``. This continues until all ``M_{\text{train}} = 200`` measurements are included (40 stages), after which the final sampling run draws ``K = 10^5`` samples without thinning for diagnostic analysis.

```julia
MH_samples, acceptance_ratio, runtime_sampling = staged_ODE_MH(
    u_t_bolus, t_train, y_train, (-T_train, 0.0),
    K, K_b, k_d,
    f_theta!, g_theta!, log_pdf_w_theta,
    log_pdf_theta, theta_0,
    log_pdf_x_init, x_init_0,
    proposal_z_cov_0, M_chunk, K_stage, alpha;
    regularizer=regularizer,
    ODE_solver=ODE_solver, ODE_solver_opts=ODE_solver_opts,
    print_progress=true
)
```

The proposal scaling ``\alpha`` should yield an acceptance rate of approximately 25%. If the observed acceptance rate deviates substantially, adjust ``\alpha`` accordingly: increase it if acceptance is too high (proposals are too conservative), decrease it if acceptance is too low (proposals overshoot).

## Determining the Thinning Factor

The ``K = 10^5`` samples from the diagnostic run above are used to determine an appropriate thinning interval. We compute the normalized autocorrelation functions (ACFs) of the model parameters and initial state components up to a maximum lag of 50:

```julia
max_lag = 50
autocorrelation = compute_autocorrelation(MH_samples; max_lag=max_lag)
```

![Autocorrelation functions of the MH samples](../assets/autocorrelation.svg)

The figure shows the normalized ACFs of the model parameters ``(p_2, p_3, n)`` (red) and the initial state components ``\boldsymbol{x}(-T_{\text{train}})`` (green). The autocorrelation of all components decays substantially within a lag of 25. Based on this, we apply a thinning factor of ``k_d = 25`` in all subsequent experiments, so that the retained samples exhibit negligible empirical autocorrelation.

The effective sample size (ESS) provides a complementary summary of sampling efficiency:

```julia
ess = compute_ess(MH_samples; max_lag=100)
@printf("Minimum ESS: %.1f (= %.2f / s)\n", minimum(ess), minimum(ess) / runtime_sampling)
```

The goal of tuning is to maximize the ESS per unit of computation time. If the ESS is low relative to the total number of iterations, this typically indicates poor mixing — consider adjusting the proposal scaling, increasing the number of samples per stage, or reducing the number of measurements added per stage.

## Posterior Validation

As a sanity check, each posterior sample is simulated forward over the test window (6 pm–12 am) and compared against the 100 held-out test measurements:

```julia
t_pred = t_test
x_pred = Array{Float64}(undef, n_x, length(t_pred), K)
for k in 1:K
    ode_rhs_pred(dx, x, p, t) = f_theta!(dx, MH_samples[k].theta, x, u_t_bolus(t), t)
    prob_pred = ODEProblem(ode_rhs_pred, MH_samples[k].x_t, t_span, MH_samples[k].theta)
    sol_pred = solve(prob_pred, ODE_solver; ODE_solver_opts..., saveat=t_pred)
    x_pred[:, :, k] .= Array(sol_pred)
end
```

![Posterior predictions vs. held-out test data](../assets/prediction.svg)

The posterior mean (solid line) should track the held-out glucose measurements reasonably well, and the prediction band (shaded region) should capture the test data without being excessively wide. If the posterior mean systematically deviates from the test data or the prediction band collapses, this suggests issues with the prior specification or insufficient mixing during sampling.

## Additional Diagnostics

The script also produces trace plots and posterior histograms for each sampled component. The trace plots should appear stationary after burn-in with no long-term trends. The posterior histograms are overlaid with the prior distributions and the true parameter values — if the data is informative, the posterior should be tighter than the prior and centered near the true value. Note, however, that these marginal plots do not show correlations between parameters; the joint posterior may still yield accurate predictions even if individual marginals appear broad.

For a more rigorous convergence assessment, the script includes an optional section that runs multiple independent chains and computes the Gelman–Rubin statistic ``\hat{R}``. Values close to 1 (typically ``\hat{R} < 1.05``) indicate good convergence across chains.