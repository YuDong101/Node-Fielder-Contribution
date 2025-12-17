# Nodal Fielder Contribution

This repository contains the MATLAB source code for analyzing the relationship between the node-wise Fiedler vector (in directed, weighted networks) and the local synchronization speed of dynamical systems.

**Reference:** Nodal Spectral Sensitivity and the Distribution of Dynamical Timescales in Directed Networks from *Yu et al., Phys. Rev. E Submission*.

---

## 1. File Structure

Please ensure the following five MATLAB files are located in the same directory. The workflow relies on the interaction between the main script and the auxiliary calculation modules.

* **`Feilder_network.m`**: The main script. It defines the network topology (Adjacency Matrix $A$), runs the coupled simulations, computes the Fiedler statistics, and generates comparison plots.
* **`fiedler_directed_weighted.m`**: A core function that computes the Fiedler value, left/right eigenvectors, and per-node contributions for **directed, weighted, non-symmetric** graphs.
* **`simulate_linear_and_lorenz.m`**: A simulation engine that simultaneously solves a Linear Diffusion system and a Coupled Lorenz network on the same graph topology.
* **`estimate_decay_tau.m`**: A tool to estimate the decay time constant ($\tau$) for the linear system using exponential fitting or $1/e$ decay methods.
* **`sync_time_lor.m`**: A tool to calculate the synchronization time for the Lorenz system based on "first-entry + dwell time" criteria relative to a threshold.

> **Note:** The main script is named `Feilder_network.m`. Do not rename these files or move them into separate subfolders, as the script calls them directly.

## 2. Usage Instructions

### 2.1 Prerequisites
* MATLAB (Recommended version: R2018b or later).
* Statistics and Machine Learning Toolbox (Recommended for robust fitting in `estimate_decay_tau`, though fallbacks exist).

### 2.2 Configuration (Optional)
In `Feilder_network.m`, you can modify the system parameters at the beginning of the script to explore different topologies or dynamics:
* **A (Adjacency Matrix)**: Uncomment different example matrices (Directed cycle, Random cycle, Hub motif) or define your own $N \times N$ matrix.
* **Coupling Strengths**:
    * `glin`: Coupling strength for the linear diffusion process (default `0.1`).
    * `glor`: Coupling strength for the Lorenz network (default `10`).
* **Simulation Parameters**:
    * `Tend`: Total simulation time.
    * `nRun`: Number of ensemble runs to average the error (for smoothing).

### 2.3 Running the Simulation
1. Open `Feilder_network.m` in the MATLAB editor.
2. Run the script by clicking the **Run** button or typing `Feilder_network` in the Command Window.

**What the program does:**
* Computes the generalized Fiedler vector and node contributions for the defined graph.
* Simulates the time evolution of both Linear and Lorenz systems.
* Extracts the synchronization time ($\tau$) for every node.
* Calculates the Pearson correlation between the Fiedler node contribution and the actual synchronization time.

### 2.4 Output
Upon completion, the program will generate multiple figures:
1. **Figure 101**: Side-by-side comparison of synchronization error decay ($E(t)$) for Linear vs. Lorenz systems, and bar charts of node-specific $\tau$.
2. **Figure 2/32**: Scatter plots and bar charts correlating the theoretical Fiedler prediction with the empirical $\tau_{lin}$ and $\tau_{Lorenz}$.
3. **Figure 4001/4002**: Heatmaps of the final system states.
4. **Console Output**: Prints the Pearson correlation coefficients ($r$) and p-values for both systems.

## 3. Notes

* **Theoretical Basis**: The `fiedler_directed_weighted.m` function calculates $L = D_{out} - A$ and selects the eigenvalue with the smallest non-zero real part ($\text{Re}(\lambda_2)$). It is specifically designed to handle **non-symmetric** matrices.
* **Numerical Stability**: The simulation uses `ode45`. If you modify the network size $N$ to be very large, computation time will increase significantly.
* **Data Fitting**: `estimate_decay_tau.m` attempts to use `nlinfit`. If the Statistics Toolbox is missing, it will automatically fallback to `lsqcurvefit` or `fminsearch`.
