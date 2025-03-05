# SBI-Comparison

This repository contains the implementation of different **Simulation-Based Inference (SBI)** methods for comparative experiments. The methods included are:

- **BayesFlow** (Amortized Bayesian workflows)
- **Flow Matching** (Applied for flow-based inference)
- **Sequential Neural Likelihood (SNL)** (Implemented in the `sbi` library)

The repository contains two experimental setups:
- `gaussian`: Gaussian likelihood experiments.
- `disease`: Disease modeling experiments.

## 📌 Repository Structure
```
SBI-comparison/ │── AFM/ # Amortized Flow Matching (BayesFlow 2.0.0) │ ├── gaussian/ # Gaussian experiment │ ├── disease/ # Disease modeling experiment │── BayesFlow/ # BayesFlow-based inference (BayesFlow 1.1.6) │ ├── gaussian/
│ ├── disease/ │── SNL/ # Sequential Neural Likelihood (SBI 0.23.2) │ ├── gaussian/ │ ├── disease/ │── README.md # Project documentation │── requirements.txt # Dependencies │── .gitignore # Ignore unnecessary files
```

## 🛠 Installation

Ensure you have Python 3.11.10 installed. To set up the environment, run:

```bash
git clone https://github.com/xudongwu-0/SBI-comparison.git
cd SBI-comparison
pip install -r requirements.txt
```
