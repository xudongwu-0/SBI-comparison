# SBI-Comparison

This repository provides implementations of various **Simulation-Based Inference (SBI)** methods for comparative experiments. The included methods are:

- **BayesFlow** (Amortized Bayesian workflows)
- **Flow Matching** (Flow-based inference)
- **Sequential Neural Likelihood (SNL)** (Implemented in the `sbi` library)

Each method contains two experimental setups:
- `gaussian`: Gaussian likelihood experiments.
- `disease`: Disease modeling experiments.


## Installation

Ensure you have Python 3.11.10 installed. To set up the environment, run:

```bash
git clone https://github.com/xudongwu-0/SBI-comparison.git
cd SBI-comparison
```

## Fully Interactive Jupyter Notebooks

All experiments in this repository are designed to be **fully interactive** and are provided as **Jupyter notebooks** (`.ipynb`). Simply open any notebook, run the code, and modify parameters as needed to explore different settings in real-time.

### Running the Notebooks

To launch the Jupyter environment, navigate to the repository and run:

```bash
pip install jupyter
jupyter notebook

```
Then, open any .ipynb file from the corresponding method folder (BayesFlow/, AFM/, or SNL/). These notebooks allow you to:

- Adjust prior settings, likelihood choices, and inference hyperparameters.
- Run experiments interactively without modifying core scripts.
- Visualize posterior distributions and calibration results in real-time.

This structure ensures that every experiment is reproducible and customizable with minimal effort.


## Package Requirements

The experiments were conducted using:

- Python: Version 3.11.10.
- BayesFlow: BayesFlow Version 1.1.6, implemented in the BayesFlow framework for amortized Bayesian workflows.
- Flow Matching: BayesFlow Version 2.0.0, applied for flow matching experiments.
- Sequential Neural Likelihood (SNL): SBI Version 0.23.2, implemented in the sbi (Simulation-Based Inference) library.
