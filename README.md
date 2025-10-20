# Analysis of Higher Education Dropouts Dynamics through Multilevel Functional Decomposition of Recurrent Events in Counting Processes

This repository provides the code and resources necessary to reproduce the analyses and experiments presented in the paper:

> **Analysis of Higher Education Dropouts Dynamics through Multilevel Functional Decomposition of Recurrent Events in Counting Processes**  
> [arXiv preprint](https://arxiv.org/abs/2411.13370)
> by A. Ragni, C. Masci, A. M. Paganoni

The repository is divided into two main components, corresponding to the **simulation study** and the **real-world case study** described in the paper.

## Code Organization


### **1. Simulation Study (`/Simulation`)**

Contains scripts to replicate the simulation experiments described in **Section 3** of the paper.  
These experiments assess the performance and robustness of the proposed decomposition model under controlled conditions.

**Files:**
- `simulation.R` — Main script to run the complete simulation pipeline for different runs over different parameter settings
- `analysis.R` — Summarizes and visualizes simulation outcomes across simulations
- `PlotsSimulation.R` - Main functions for plotting the simulation output and reproducing the Figures in Section 3
- `simulate_once.R` — Runs a single instance of the simulation (one run, for a particular parameter set up) 
- `utils.R` — Utility functions for `simulate_once.R`

### **2. Case Study (`/Case study`)**

Includes scripts for the empirical analysis presented in **Section 4** of the paper.  
The case study focuses on modeling and predicting dropout dynamics in higher education using the proposed methodology.

**Files:**
- `1_DataPreprocessing.R` — Data cleaning, preparation, and transformation into event-counting format.  
- `2_FromCompensatorsToPCA.R` — Converts estimated compensator functions into a form suitable for functional PCA.  
- `3_Prediction.R` — Implements predictive modeling based on extracted functional components.  
- `4_Bootstrap.R` — Conducts bootstrap analyses for inference and uncertainty quantification.

