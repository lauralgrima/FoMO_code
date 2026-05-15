AQUA: Adaptive Rate Q-learning of Utilities for Apportioning Actions

AQUA is a reinforcement learning framework developed to explain rapid multi-option foraging behavior in freely moving mice navigating large spatial environments. The model is a simple implementation of an RL agent that combines principles from reinforcement learning, optimal foraging theory, and adaptive learning-rate optimization to capture how animals rapidly learn and exploit multiple reward sources.

This repository contains MATLAB code for simulating and fitting the AQUA model on concatenated behavioral sessions from the FoMO (Foraging across Multiple Options) task described in:

Grima, L. L., Guo, Y., Narayan, L., Hermundstad, A. M., & Dudman, J. T. (2026). A global dopaminergic learning rate enables adaptive foraging across many options. Neuron. Open Access.
 https://www.sciencedirect.com/science/article/pii/S0896627326002795

⸻

Overview

Animals foraging in natural environments must continually decide whether to:

* remain at a currently sampled option,
* transition to another option sensitive to the distance to the option,
* update beliefs about option quality (probability of reward) from experience, 
* balance exploration against exploitation.

AQUA addresses this problem using four key components:

1. Separate stay vs. leave decisions
    * The model independently estimates the utility of remaining at the current option versus transitioning to that option from elsewhere.
2. Dynamic global learning rate
    * Learning rates decay adaptively across experience, enabling rapid early learning and stable long-term behavior.
3. Online reward probability estimation
    * Option values are continuously updated from experienced rewards.
4. Distance-weighted utility computation
    * Choice probabilities are scaled by movement costs between spatially distributed options.
5. Option choice is proprotional to expected utility with some exploration
    * We refer to this action selection policy as epsilon-proportional (in analogy to epsilon-greedy).

The model reproduces:

* matching behavior,
* rapid learning,
* stable multi-option allocation,
* session-to-session adaptation,
* and experimentally observed dopamine-related learning dynamics.

⸻

Repository Contents

.
├── HX_model_session_forAlphaConcat.m
├── data/
├── analysis/
├── figures/
├── utilities/
└── README.md

Core files

HX_model_session_forAlphaConcat.m

Example implementation of the AQUA model applied to concatenated behavioral sessions.

This function:

* simulates choice allocation across six spatially distributed options,
* updates reward beliefs dynamically,
* incorporates travel-cost weighting,
* applies adaptive learning-rate schedules,
* and compares simulated behavior to empirical data.

Inputs include:

Variable	Description
alpha	Dynamic learning-rate vector across visits
visit_matrix	Empirical visit timestamps and port identities
cost_per_port	Pairwise movement-cost matrix
rew_sched	Ground-truth reward schedule
income	Empirical reward/income trace
session_ids	Session boundaries for concatenated data

Outputs include:

Variable	Description
trans_r2	Transition matrix similarity metric
income_r2	Income prediction accuracy
p_reward	Estimated reward probability for each option
income_model	Simulated reward income trace

⸻

Installation

Requirements

* MATLAB R2023b or newer (recommended)

No external toolboxes beyond standard MATLAB functionality are currently required.
The TONIC repository (custom code and scripts from DudmanLab) can be useful. 
If you would like a copy please contact.

⸻

Quick Start

Example usage of the model

[trans_r2, income_r2, visits_for_LL, rewards_for_LL, p_reward, income_model] = ...
    HX_model_session_forAlphaConcat(
        alpha,
        visit_matrix,
        cost_per_port,
        rew_sched,
        income,
        session_ids);

A script illustrating use of the minimal workflow is HX_final_script.m

1. Load behavioral session data.
2. Construct reward schedules and visit matrices (using HX_load* functions).
3. Define movement-cost matrix.
4. Define adaptive learning-rate schedule (akin to a learning rate scheduler or can be derived from optimization).
5. Run AQUA simulation.
6. Compare simulated behavior against empirical behavior using the output variables.

⸻

Model Description

Reward update rule

AQUA updates reward expectations according to:

P(R|O_n)_t = \alpha(\tau) \cdot R_t(O_n) + (1 - \alpha(\tau)) \cdot P(R|O_n)_{t-1}

where:

* P(R|O_n) is the estimated reward probability for option n,
* R_t(O_n) is the experienced reward outcome,
* α(τ) is a dynamic global learning rate over experience.

Utility-weighted action selection

Choice probabilities are normalized by spatial transition costs:

U(O_n) = \frac{P(R|O_n)}{D(O_c, O_n)}

where:

* U(O_n) is option utility,
* D(O_c, O_n) is the movement cost from the current option to the target option.

Exploration policy

The model includes an epsilon-greedy exploration policy that relaxes:

\epsilon(t) = \epsilon_0 + e^{-t/\tau}

allowing high exploration early in learning and more stable exploitation later.

⸻

Behavioral Paradigm

The FoMO task consists of:

* six spatially distributed lick spouts,
* deterministic interval reward schedules,
* freely moving mice,
* unconstrained transitions between options,
* and no explicit trial structure.

Mice rapidly learn to allocate choices proportionally to reward probability while incorporating movement costs between options.

The AQUA model was specifically developed to account for this form of many-option spatial foraging behavior.

⸻

Reproducing Results

The original study examined:

* matching behavior,
* learning dynamics,
* transition matrices,
* reward collection efficiency,
* adaptive learning-rate estimation,
* and dopamine-related modulation of learning.

Typical analyses include:

* transition matrix similarity (r²),
* reward-income prediction,
* matching sensitivity,
* KL-divergence over learning,
* and session-to-session adaptation.

⸻

Data Format

Behavioral data are expected to be organized as:

* visit-by-time matrices,
* reward schedule matrices,
* movement-cost matrices,
* and session-index vectors.

Additional preprocessing scripts may be added in future updates.

⸻

Citation

If you would like to cite our original description of AQUA in published work:

@article{GRIMA2026,
title = {A global dopaminergic learning rate enables adaptive foraging across many options},
journal = {Neuron},
year = {2026},
issn = {0896-6273},
doi = {https://doi.org/10.1016/j.neuron.2026.04.010},
url = {https://www.sciencedirect.com/science/article/pii/S0896627326002795},
author = {Laura L. Grima and Yipei Guo and Lakshmi Narayan and Ann M. Hermundstad and Joshua T. Dudman},
}

⸻

Related Concepts

AQUA combines ideas from:

* temporal-difference reinforcement learning,
* Q-learning,
* matching theory,
* marginal value theory,
* adaptive optimization,
* and spatial foraging theory.

⸻

Future Directions

Planned additions may include:

* Python implementation,
* model fitting pipelines,
* behavioral preprocessing utilities,
* detailed descriptions of photometry analysis integration,
* probabilistic schedule simulations,
* and expanded benchmarking against alternative RL models.

⸻

License

CC-BY 4.0

⸻

Contact

For questions regarding the model or experimental paradigm:

* Laura L. Grima — grimal@janelia.hhmi.org
* Joshua T. Dudman — dudmanj@janelia.hhmi.org

⸻

Acknowledgments

This work was developed at the Janelia Research Campus, Howard Hughes Medical Institute
where Joshua T. Dudman is a Senior Group Leader.