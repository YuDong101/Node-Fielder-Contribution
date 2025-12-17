\documentclass[a4paper,11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{geometry}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{titlesec}
\usepackage{enumitem}

% Setup geometry
\geometry{left=2.5cm, right=2.5cm, top=2.5cm, bottom=2.5cm}

% Custom section formatting to mimic Markdown headers
\titleformat{\section}{\Large\bfseries\sffamily}{}{0em}{}[\hrule]
\titleformat{\subsection}{\large\bfseries\sffamily}{}{0em}{}

% Define colors
\definecolor{codebg}{RGB}{245,245,245}
\definecolor{linkcolor}{RGB}{0,0,238}

\title{\textbf{Nodal Spectral Sensitivity and the Distribution of Dynamical Timescales in Directed Networks}}
\author{Dong Yu$^*$, Yu-kang Shu, Zi-yi Chen, Hang Chen, Yi-rui Tang, Jiang Long, and Qi-ming Pei$^*$ \\ \textit{School of Physics and Optoelectronic Engineering, Yangtze University}}
\date{}

\begin{document}

\maketitle

[cite_start]\noindent This repository contains the MATLAB source code for the research paper: \textbf{Nodal Spectral Sensitivity and the Distribution of Dynamical Timescales in Directed Networks}[cite: 1, 3].

\vspace{1em}
\hrule
\vspace{1em}

\section*{1. File Structure}

Please ensure the following three MATLAB files are located in the same directory. [cite_start]The program relies on the interaction between the main script and these auxiliary functions to analyze structural sensitivity and dynamical relaxation[cite: 35, 173].

\begin{itemize}[leftmargin=*]
    \item \textbf{\texttt{main.m}}: The main script. [cite_start]It initializes network topologies (Watts-Strogatz or Tunable Scale-Free), computes spectral properties, iterates through simulation conditions, and controls the main execution loop[cite: 36, 615].
    
    \item \textbf{\texttt{calc\_spectral\_sensitivity.m}}: A function module used to compute the \textbf{Nodal Fiedler Contribution ($I_{n_i}$)}. [cite_start]It diagonalizes the Laplacian matrix to obtain the left and right eigenvectors associated with the algebraic connectivity ($\lambda_2$)[cite: 35, 52].
    
    [cite_start]\item \textbf{\texttt{dynamics\_solver.m}}: A function module used to numerically integrate the differential equations for \textbf{Linear Consensus} or \textbf{Coupled Lorenz Oscillators} and extract the characteristic relaxation timescales ($\tau_i$)[cite: 87, 190].
\end{itemize}

\begin{quote}
\textbf{Note:} Do not rename these files or move them into separate subfolders, as \texttt{main.m} calls the functions directly from the current working directory.
\end{quote}

\section*{2. Usage Instructions}

\subsection*{2.1 Prerequisites}
\begin{itemize}
    \item MATLAB (Recommended version: R2018b or later).
    \item Ensure the current working directory in MATLAB is set to the folder containing these files.
\end{itemize}

\subsection*{2.2 Configuration (Optional)}
In \texttt{main.m}, system parameters are defined within a structure (e.g., \texttt{params}). [cite_start]You may modify these values to explore different topological and dynamical regimes[cite: 44, 260]:

\begin{itemize}
    \item \textbf{Topological parameters:}
    \begin{itemize}
        [cite_start]\item \textbf{$N$}: Number of nodes (e.g., $N=100$)[cite: 44].
        [cite_start]\item \textbf{$p$}: Rewiring probability for Watts-Strogatz networks ($0 \le p \le 1$)[cite: 432].
        [cite_start]\item \textbf{$\alpha$ (alpha)}: Generative parameter interpolating between Growing Random Networks ($\alpha=0$) and Scale-Free networks ($\alpha=1$)[cite: 616].
    \end{itemize}
    
    \item \textbf{Dynamical parameters:}
    \begin{itemize}
        [cite_start]\item \textbf{$g$}: Coupling strength ($g_{lin}$ for linear consensus, $g_{lor}$ for Lorenz systems)[cite: 163, 433].
        [cite_start]\item \textbf{Time span}: Simulation duration (typically $t \in [0, 100]$)[cite: 175].
    \end{itemize}
\end{itemize}

\subsection*{2.3 Running the Simulation}
\begin{enumerate}
    \item Open \texttt{main.m} in the MATLAB editor.
    \item Run the script by clicking the \textbf{Run} button or typing \texttt{main} in the Command Window.
\end{enumerate}

\noindent \textbf{What the program does:}
\begin{itemize}
    \item Iterates over the defined network configurations (varying $p$ or $\alpha$).
    [cite_start]\item Calls \texttt{calc\_spectral\_sensitivity.m} to calculate the nodal Fiedler contribution $I_{n_i}$ via eigenvalue perturbation theory[cite: 35, 79].
    [cite_start]\item Calls \texttt{dynamics\_solver.m} to simulate the time evolution of the system states[cite: 173].
    [cite_start]\item Computes the normalized error trajectories $\sigma_i(t)$ and fits the relaxation timescales[cite: 187].
\end{itemize}

\subsection*{2.4 Output}
Upon completion, the program will generate visualization data and display figures corresponding to the manuscript's results:
\begin{enumerate}
    [cite_start]\item \textbf{Error Trajectories}: Plots the temporal evolution of synchronization error $\sigma_i(t)$, color-coded by the rank of $I_{n_i}$ to illustrate the ``Anchor-Bottleneck'' mechanism[cite: 37, 272].
    [cite_start]\item \textbf{Timescale Correlation Scatter}: Visualizes the relationship between the structural metric $I_{n_i}$ and the dynamical timescales ($\tau_i$), showing the linear fit and correlation coefficients[cite: 437].
\end{enumerate}

\section*{3. Notes}

\begin{itemize}
    \item \textbf{Path Settings:} Ensure MATLAB’s ``Current Folder'' matches the repository location; otherwise, the script will fail to locate the function files.
    [cite_start]\item \textbf{Numerical Precision:} The differential equations are integrated using a standard adaptive Runge-Kutta scheme with relative tolerance $10^{-6}$ and absolute tolerance $10^{-8}$[cite: 174]. Modifying these tolerances may affect the accuracy of the tail distributions.
    [cite_start]\item \textbf{Performance:} Sweeping a large ensemble (e.g., $M=100$ realizations with $N_{run}=1000$) can be computationally intensive[cite: 175]. For quick tests, consider reducing the number of realizations in \texttt{main.m}.
\end{itemize}

\end{document}
