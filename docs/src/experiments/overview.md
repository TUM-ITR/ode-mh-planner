# [Overview](@id experiments)

This page describes the simulation setup used to evaluate the proposed Bayesian inference and control framework. For the simulation we consider glucose regulation in Type 1 diabetes. This setting provides a natural test case for the proposed method: the dynamics are nonlinear and patient-specific, only the glucose level is directly measurable, and measurements typically occur at irregular sampling times. At the same time, glucose regulation is a safety-critical task: inadequate control actions can induce hypo- or hyperglycemia, both of which carry significant medical risks. These characteristics make this setting well-suited for evaluating the proposed framework, which infers parameters and latent states from infrequent measurements, quantifies the associated uncertainty, and uses this information to compute control inputs that remain robust under the inferred uncertainty. An in-depth description of the setup is given below; results are presented on the subpages linked at the bottom.

All experiment scripts are located in the `experiments/` directory of the repository.

## Glucose–Insulin Dynamics

We employ a version of the Bergman minimal model ([Bergman et al., 1981](https://doi.org/10.1172/jci110398)), extended with meal disturbances and exogenous insulin infusion as in [Ali et al., 2011](https://doi.org/10.1002/oca.920). The model captures the interaction between blood glucose ``G(t)`` [mg/dL], plasma insulin ``I(t)`` [mU/L], and a remote insulin action state ``X(t)`` [1/min] that represents the delayed effect of insulin on glucose uptake. The continuous-time dynamics are given by

```math
\begin{aligned}
\dot{G}(t) &= -p_1\left(G(t) - G_b\right) -X(t) G(t) + D(t),\\
\dot{X}(t) &= -p_2 X(t) + p_3\left(I(t) - I_b\right),\\
\dot{I}(t) &= -n\left(I(t) - I_b\right) + u(t),
\end{aligned}
```

where ``t`` denotes time (in minutes), the parameters ``p_1``, ``p_2``, and ``p_3`` characterize insulin action dynamics, ``n`` is the insulin clearance rate, and ``G_b`` and ``I_b`` denote basal levels. The input ``u(t)`` represents exogenous insulin infusion delivered by an insulin pump, while ``D(t)`` denotes glucose appearance resulting from meals. We consider a Type~1 diabetes scenario, where the glucose effectiveness parameter ``p_1`` is physiologically negligible and therefore omitted, resulting in ``\dot{G}(t) = -X(t) G(t) + D(t)`` ([Ali et al., 2011](https://doi.org/10.1002/oca.920)).

## Meal Disturbances

Meal-induced glucose appearance is modeled using an exponential profile

```math
D(t) = S_{\text{meal}}\, B\, e^{-B\,(t - t_{\text{meal}})}, \quad t \geq t_{\text{meal}},
```

where ``S_{\text{meal}}`` is the meal size, ``t_{\text{meal}}`` is the meal time, and ``B = 0.05`` min``{}^{-1}`` controls the decay rate. Three meals are administered during the 18-hour simulation window (6 am–12 am):

| Time   | Meal size ``S_{\text{meal}}`` [mg/dL] |
|:-------|--------------------------------------:|
| 8 am   | 60                                    |
| 1 pm   | 90                                    |
| 7 pm   | 80                                    |

In practice, meal timing and size are rarely known precisely, but relaxing this assumption is beyond the scope of the present study.

## Simulation Protocol

The system is simulated over an 18-hour period (6 am–12 am) split into two phases. The first 12 hours (6 am–6 pm) serve as the training window and include the meals at 8 am and 1 pm. During this phase, a simple insulin dosing strategy is applied: each meal triggers an insulin amount proportional to the meal size, delivered over a one-hour window. This keeps the simulated trajectory physiologically plausible without requiring a controller.

The remaining 6 hours (6 pm–12 am) constitute the control horizon and contain the 7 pm meal. This phase is used to evaluate the control strategy computed from the posterior samples obtained during the training window.

The ground-truth trajectory is simulated using the Tsitouras 5/4 method (the default adaptive solver in [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)).

## Measurement Model

In practice, glucose measurements can occur at infrequent and irregular time points due to sensor limitations and intermittent reporting. To mimic this setting, we sample ``M = 200`` measurement instants uniformly at random from the training window. The corresponding glucose observations are corrupted by additive zero-mean Gaussian noise:

```math
y_m = G(t_m) + v_m, \quad v_m \sim \mathcal{N}(0, \sigma^2),
```

with standard deviation ``\sigma = 8`` mg/dL. Only glucose is directly measured; the states ``X(t)`` and ``I(t)`` remain latent.

## Prior Distributions

A key advantage of the proposed method is that the prior is formulated directly in state-space form. Since physiological knowledge about glucose–insulin regulation is typically available in this representation, models such as the Bergman minimal model can be incorporated naturally, and only the model parameters must be inferred from data. We assume the basal insulin level is known and fix it at ``I_b = 7`` mU/L ([Ali et al., 2011](https://doi.org/10.1002/oca.920)) and estimate the parameters ``p_2``, ``p_3``, and ``n`` together with the latent state trajectory.

To ensure positivity and reflect prior uncertainty around physiologically plausible values, we place independent lognormal priors on the model parameters. The prior means and variances are chosen such that approximately 95% of the prior mass lies within the physiological ranges reported in [Ali et al., 2011](https://doi.org/10.1002/oca.920). The initial state at 6 am is modeled as Gaussian and centered at physiologically plausible fasting levels. The exact prior distributions are:

| Parameter  | Prior distribution                                  |
|:-----------|:----------------------------------------------------|
| ``p_2``    | ``\text{LogNormal}(-4.26,\; 0.18^2)``              |
| ``p_3``    | ``\text{LogNormal}(-13.27,\; 0.28^2)``             |
| ``n``      | ``\text{LogNormal}(-1.66,\; 0.23^2)``              |
| ``G(-T)``  | ``\mathcal{N}(80.0,\; 8.0^2)``                     |
| ``X(-T)``  | ``\mathcal{N}(0.0,\; 0.001^2)``                    |
| ``I(-T)``  | ``\mathcal{N}(7.0,\; 2.0^2)``                      |

## Inference Setup

The MH sampler (Algorithm 1 in the paper) uses a fourth-order Runge–Kutta integrator (RK4) with a step size of 0.5 min for propagating the ODE dynamics within each likelihood evaluation. To improve mixing, we employ an adaptive proposal strategy (Remark 1 in the paper): the sampler starts with 5 measurements and increases the number of data points by 5 per stage. At each stage, the proposal covariance is adapted based on the empirical covariance of the chain, allowing it to track the evolving posterior geometry.

The thinning interval is determined empirically by running the sampler for ``10^5`` iterations on a representative instance and inspecting the autocorrelation functions (ACFs) of the model parameters and initial state; see [Inference and Sampler Tuning](@ref sampling) for details. Based on this analysis, a thinning factor of ``k_d = 25`` is applied in all subsequent experiments.

Using this thinning factor, the sampler is then run to produce ``K = 100`` posterior samples for use in the scenario-based optimal control problem.

## Optimal Control Problem

The control objective is to maintain glucose close to a target of ``G_{\text{ref}} = 80`` mg/dL while compensating for meal disturbances and respecting safety constraints. The cost functional over the prediction horizon ``H = 6`` hours is

```math
J = W_{G_f} \left(G(H) - G_{\text{ref}}\right)^2 + \int_0^H \left( W_G \left(G(t) - G_{\text{ref}}\right)^2 + W_U\, u(t)^2 \right) dt,
```

with weights ``W_G = 1``, ``W_{G_f} = 10``, and ``W_U = 10^{-3}``. The objective is thus dominated by glucose regulation, while the regularization term discourages unnecessarily aggressive insulin delivery.

To reduce the risk of hypo- and hyperglycemia and thus ensure patient safety, we enforce the clinically recommended glucose bounds

```math
70 \;\text{mg/dL} \leq G(t) \leq 180 \;\text{mg/dL},
```

and the insulin pump limits ``0 \leq u(t) \leq 20`` mU/min.

The OCP is discretized using the same RK4 scheme (step size 0.5 min) and solved with [JuMP](https://jump.dev/), [IPOPT](https://coin-or.github.io/Ipopt/), and the [HSL](https://www.hsl.rl.ac.uk/) linear solvers.

## Nominal Baseline

For comparison, we evaluate a nominal-model baseline. In this approach, the parameters are fixed at the geometric mean of their respective prior distributions, and the latent state is estimated using an extended Kalman filter (EKF) equipped with the same RK4 integration scheme to accommodate infrequent measurements. The optimal input trajectory is then computed using this nominal model and the EKF state estimate as the initial condition. This baseline represents a standard approach that does not account for parametric or state uncertainty beyond point estimation.

## Experiment Subpages

Results and analysis are organized across the following subpages:

- **[Inference and Sampler Tuning](@ref sampling):** Sampler configuration, autocorrelation analysis, thinning factor selection, and posterior validation.
- **[Optimal Control](@ref optimal-control):** Single-run demonstration of the full inference-to-control pipeline, with trajectory comparisons against the nominal baseline.
- **[Monte Carlo Study](@ref monte-carlo):** Statistical evaluation over 100 independent runs, including cost and constraint violation statistics.
