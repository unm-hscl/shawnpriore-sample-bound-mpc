# Sample Bound MPC
Code for the paper, "Chance Constrained Stochastic Optimal Control Based on Sample Statistics With Almost Surely Probabilistic Guarantees," Under Review, IEEE TAC.

arXiv: [https://arxiv.org/abs/2303.16981](https://arxiv.org/abs/2303.16981)

## Requirements
* CVX [http://cvxr.com/cvx/](http://cvxr.com/cvx/)

## Examples
### Disturbance from Gravitational Effects
We use a 6d CWH system with additive noise representing the effect of The J2 gravitational harmionic, and solar and lunar gravity to model the dynamics of a satellite rendezvous operation. This example compares the efficiency of the proposed approach with other sampling based methods.

### Gaussian Disturbance
We use a 6d CWH system with additive Gaussian noise to model the dynamics of a three satellite rendezvous operation. This example compares the efficacy of the proposed approach with with the Cantelli's inequality (the analytic counterpart to sample based concentration inequality derived in this work).
