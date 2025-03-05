# SBI-Comparison

This repository contains the implementation of different **Simulation-Based Inference (SBI)** methods for comparative experiments. The methods included are:

- **BayesFlow**
- **Flow Matching** 
- **Sequential Neural Likelihood (SNL)** 

The repository contains two experimental setups:
- `gaussian`: Gaussian likelihood experiments.
- `disease`: Disease modeling experiments.


## Installation

Ensure you have Python 3.11.10 installed. To set up the environment, run:

```bash
git clone https://github.com/xudongwu-0/SBI-comparison.git
cd SBI-comparison
```

## Simulation Environment

The experiments were conducted using:

BayesFlow: BayesFlow Version 1.1.6, implemented in the BayesFlow framework for amortized Bayesian workflows.
Flow Matching: BayesFlow Version 2.0.0, applied for flow matching experiments.
Sequential Neural Likelihood (SNL): SBI Version 0.23.2, implemented in the sbi (Simulation-Based Inference) library.
Python: Version 3.11.10.
