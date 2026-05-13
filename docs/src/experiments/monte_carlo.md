# [Monte Carlo Study](@id monte-carlo)

This page presents the statistical evaluation of the proposed MH-scenario method across 100 independent simulation runs. The Monte Carlo study assesses robustness with respect to variability in patient dynamics and initial conditions by repeatedly sampling new patients from the prior and running the full inference-to-control pipeline on each.

The experiment scripts are located in `experiments/montecarlo/`. Detailed instructions for running the study on an HPC cluster using Docker, Apptainer, and Slurm are provided in `experiments/montecarlo/README.md`.

## Experiment Design

Each of the 100 runs proceeds as follows:

1. A true patient is sampled by drawing the physiological parameters ``(p_2, p_3, n)`` and the initial state from the prior distributions (see [Overview](@ref experiments)).
2. The ground-truth trajectory is simulated over the training window (6 am–6 pm), and ``M = 200`` noisy glucose measurements are generated at random times.
3. The staged MH sampler produces ``K = 100`` posterior samples with thinning factor ``k_d = 25``.
4. The scenario-based OCP is solved using the posterior samples, and the resulting control input is applied to the true system over the 6-hour control horizon (6 pm–12 am).
5. The realized glucose trajectory is evaluated for cost and constraint satisfaction.

The same procedure is repeated for the Nominal + EKF baseline, which uses the prior geometric mean as parameters and an extended Kalman filter for state estimation.

## Results

| Method | Cost (×10⁴) | Violations |
|:-------|------------:|-----------:|
| MH-scenario | 8.98 ± 3.50 | 0 / 100 |
| Nominal + EKF | 16.55 ± 13.27 | 44 / 100 |

The Nominal + EKF baseline violates the prescribed glucose safety bounds (70–180 mg/dL) in 44 of the 100 runs, highlighting its inability to maintain safe glucose levels under patient variability. In contrast, the MH-scenario approach satisfies all constraints in every run, demonstrating strong robustness to uncertainty in both the model parameters and the latent state.

The proposed method also achieves substantially lower cost (mean and standard deviation). By jointly inferring the state and model parameters and optimizing over multiple posterior scenarios, it obtains more accurate predictions than a single nominal model, resulting in more effective control actions and improved overall performance.

## Computation Times

The table below reports the computation times for the inference and planning steps, averaged over the 100 Monte Carlo runs. Each run was executed on a single CPU core of a dual AMD EPYC 7542 (32-core) system.

| Step | Time |
|:-----|-----:|
| Inference (staged MH) | 271 ± 28 s |
| Optimal control (scenario OCP) | 753 ± 144 s |
| Total | 1024 ± 149 s |

The total computation time of approximately 17 minutes per run is not yet sufficient for real-time application, even for slowly varying dynamics such as glucose–insulin regulation. However, the current implementation is not optimized for runtime, and the results demonstrate that the method produces reliable and safe control inputs across a wide range of patient models. Reducing computation times to a practically viable level — through solver warm-starting, parallelization, or more efficient sampling strategies — is an important direction for future work.