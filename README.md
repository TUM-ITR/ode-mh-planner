<p align="center">
<img width="800" height="100" src="assets/logo_with_text.png">
</p>

[![Dev](https://img.shields.io/badge/docs-dev-blue?logo=Julia&logoColor=white)](https://TUM-ITR.github.io/ode-mh-planner/dev/)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.18183541-blue)](https://doi.org/10.5281/zenodo.18183541)

**OdeMHPlanner** provides a framework for **uncertainty-aware learning and planning** in dynamical systems with unknown dynamics and infrequent output measurements. 

It implements the method described in:

> **Learning Dynamics from Infrequent Output Measurements for Uncertainty-Aware Optimal Control**
> *Robert Lefringhausen, Theodor Springer, Sandra Hirche*
> arXiv:2512.08013 (2025)
> <https://arxiv.org/abs/2512.08013>

## Overview

The package targets control problems in which the system dynamics are unknown and the state is only partially observed through infrequent output measurements. Instead of relying on a single identified model, OdeMHPlanner explicitly represents uncertainty over both the dynamics and the latent state trajectory.

The approach combines:

1. Bayesian learning via a Metropolis–Hastings (MH) sampler equipped with a numerical ODE solver to sample from the posterior distribution over dynamics and latent states.
2. Scenario-based optimal control, where posterior samples are propagated to compute control inputs that are robust to model uncertainty.

This enables principled uncertainty quantification and safer decision-making compared to point-estimate-based methods.

For detailed instructions and examples, please refer to the **[Documentation](https://TUM-ITR.github.io/ode-mh-planner/dev/)**.

## Installation

This package is not registered. Clone and instantiate it locally:

```bash
git clone https://github.com/TUM-ITR/ode-mh-planner.git
cd ode-mh-planner
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Solvers

The results reported in the paper were obtained using proprietary **HSL linear solvers** (in particular `ma57`) for improved numerical performance. HSL solvers are available from <https://www.hsl.rl.ac.uk/> and can be obtained under a free academic license; installation instructions are provided on their website.

If HSL solvers are not available, the code automatically falls back to standard open-source solvers (e.g. MUMPS). While functional, this configuration may be slower and has not been exhaustively tested.

## Erratum

After the paper had been accepted and published, we found a bookkeeping error in the pseudocode of the Metropolis–Hastings sampler (Algorithm 1 in the paper). The error was also present in the code up to v0.2.1: the chain index was only advanced on acceptance, so rejected proposals were not recorded as repeated states. As a result, the stored samples are not guaranteed to follow the posterior, as states with a high acceptance probability are over-represented; in our experiments this bias turned out to be small. The description of the sampler in the text of the paper is correct. The corrected pseudocode is given below: every iteration produces a sample (line 7), and burn-in and thinning refer to MH iterations.

```text
Algorithm 1 (corrected): Metropolis–Hastings sampler

Require: observations D, model f_θ(·), g_θ(·), p_W(· | θ),
         priors p(θ), p(x(-T)), proposal density q(z' | z),
         initial value z[1], number of samples K,
         burn-in K_b, thinning parameter k_d
Ensure:  K samples from p(θ, x(0) | D)

 1: set N ← K_b + 1 + (K - 1)(k_d + 1)
 2: compute likelihood p(y_{1:M} | z[1])
 3: for k = 1, ..., N - 1 do
 4:     propose z' ~ q(z' | z[k])
 5:     simulate latent trajectory under z' and evaluate p(y_{1:M} | z')
 6:     compute acceptance probability α
 7:     set z[k+1] ← z' with probability α,
        otherwise z[k+1] ← z[k]
 8: end for
 9: discard burn-in and perform thinning
10: compute x[k](0) = Φ(0; z[k], u(·)) for each sample
```

The fix is included from v0.3.0. Since burn-in and thinning now refer to MH iterations
rather than accepted proposals, the sampler settings were adjusted accordingly (burn-in
`K_b = 800`, samples per stage `K_stage = 2000`, thinning `k_d = 75`), and all experiments
were rerun with the corrected sampler. Updated results are available in the
[documentation](https://tum-itr.github.io/ode-mh-planner/dev/); the qualitative conclusions
of the paper are unchanged. The code used for the published results is available under tag
[v0.2.1](https://github.com/TUM-ITR/ode-mh-planner/tree/v0.2.1).

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@article{lefringhausen2025learning,
  title={Learning Dynamics from Infrequent Output Measurements for Uncertainty-Aware Optimal Control},
  author={Lefringhausen, Robert and Springer, Theodor and Hirche, Sandra},
  journal={arXiv preprint arXiv:2512.08013},
  year={2025}
}
```
