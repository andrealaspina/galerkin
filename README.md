# ***galerkin***: A Versatile Finite Element Framework for MATLAB

<p align="center">
  <img src="https://github.com/user-attachments/files/18724203/logo.pdf" width="300">
</p>

## Overview 🌍

***galerkin*** is a powerful and extensible [**MATLAB**](https://www.mathworks.com/products/matlab.html) framework designed for developing and testing advanced **finite element** formulations. It provides an integrated environment for simulating a variety of **2D/3D** **linear/nonlinear** **single/multi-physics** **single/multi-scale** problems in the **time/frequency** domain.

During the early stages of my PhD, I struggled with the complexity of large-scale research codes and faced significant challenges in implementing advanced finite element formulations for fluid-structure interaction. This motivated me to create an **accessible** and **flexible** framework that facilitates rapid **development**, **testing**, and **deployment** of novel numerical methods. This code is particularly beneficial for both early-career and experienced researchers striving to advance the field of **computational science and engineering**.

## Features ✨

***galerkin*** comes with a range of features designed to streamline the finite element simulation process:

- Automated workflow covering the main **pre-processing**, **processing**, and **post-processing** tasks.

- Support for multiple finite element **discretizations**, including:

	- Continuous Galerkin (**CG**)

	- Hybridizable Discontinuous Galerkin (**HDG**)

	- Face-Centred Finite Volume (**FCFV**)

- Applicability to a broad range of **physical problems**.

- Seamless coupling of an arbitrary number of sub-problems for **multi-physics applications**.

- Strong emphasis on **high-order methods**.

- **Parallel computing** capabilities (limited) via MATLAB's `parfor`.

- Easy integration with **external tools** such as [**GMSH**](https://gmsh.info), [**ParaView**](https://www.paraview.org), [**DistMesh**](https://doi.org/10.1137/S0036144503429121), [**Advanpix**](https://www.advanpix.com), and [**Comsol**](https://www.comsol.com).

## Implemented Formulations 📚

***galerkin*** includes implementations for an extensive range of physics-based formulations. For mathematical details, refer to the `*.mlx` files in the `formulations/` folder or the references at the bottom of this file.

Here is a non-exhaustive list of the currently implemented **formulations**:

### Thermal Problems 🌡️

- `Thermal_CG`: Solves the **heat equation** with the CG method [^1].

- `Thermal_HDG`: Solves the **heat equation** with the HDG method [^1].

### Structural Mechanics 🔩

- `Elasticity_CG`: Solves the **equations of elastodynamics** with the CG method for linear and nonlinear (St. Venant-Kirchhoff and Neo-Hookean) elasticity models [^1].

- `ElasticityLinear_HDG`: Solves the **equations of linear elastodynamics** with the HDG method [^1].

- `ElasticityModal_CG`: Solves the eigenvalue problem of **linear elasticity** for modal analysis using the CG method.

### Fluid Dynamics 🌊

- `CompressibleFlow_HDG`: Solves the **compressible Euler equations** with the HDG method.

- `WeaklyCompressibleFlowDM_HDG`: Solves the **Navier-Stokes equations** for weakly compressible flows with the HDG method using a density-momentum formulation [^3].

- `WeaklyCompressibleFlowVP_HDG`: Solves the **Navier-Stokes equations** for weakly compressible flows with the HDG method using a velocity-pressure formulation [^3].

- `IncompressibleFlow_FCFV`: Solves the **Navier-Stokes equations** for incompressible flows with the FCFV method [^8].

### Electromagnetics 🧲

- `Electromagnetic_HDG_fast`: Solves the **Maxwell equations** in the time domain with the HDG method [^5].

- `ElectromagneticPML_HDG_fast`: Solves the **Maxwell equations** in the time domain with the HDG method including perfectly matched layers [^5].

- `ElectromagneticFD_HDG_fast`: Solves the **Maxwell equations** in the frequency domain with the HDG method [^5].

- `ElectromagneticPML_FD_HDG_fast`: Solves the **Maxwell equations** in the frequency domain with the HDG method including perfectly matched layers [^5].

- `Magnetic_HDG`: Solves the **magnetic induction equation** with the HDG method [^4].

- `MagneticCURLCURL_HDG`: Solves the **magnetic induction equation** with the HDG method using a curl-curl formulation [^4].

### Plasma Physics 🔥

- `MagnetohydrodynamicsCURL_HDG`: Solves the **equations of magnetohydrodynamics** with the HDG method [^4].

- `MagnetohydrodynamicsCURLCURL_HDG`: Solves the **equations of magnetohydrodynamics** with the HDG method using a curl-curl formulation [^4].

- `Plasma1FluidElectromagneticAdvanced_HDG`: Solves the **Euler–Maxwell plasma equations** with the HDG method including an electrostatic initialization and a projection-based divergence correction method to enforce the Gauss laws [^7].

### Porous Media 🧽

- `Darcy_CG`: Solves the **Darcy law** with the CG method for the macroscopic problem, informed by effective properties from a unit cell problem solving the **Navier-Stokes equations** with Comsol's LiveLink for MATLAB [^6].

- `Darcy2Phase_CG`: Solves the **two-phase Darcy law** with the CG method for the macroscopic problem, informed by effective properties from a unit cell problem solving the **Cahn-Hilliard-Navier-Stokes equations** with Comsol's LiveLink for MATLAB.

- `Darcy2PhaseRichards_CG`: Solves the **Richards’ equation** with the CG method for the macroscopic problem, informed by effective properties from a unit cell problem solving the **Cahn-Hilliard-Navier-Stokes equations** with Comsol's LiveLink for MATLAB.

- `Darcy2PhaseRichards_HDG`: Solves the **Richards’ equation** with the HDG method for the macroscopic problem, informed by effective properties from a unit cell problem solving the **Cahn-Hilliard-Navier-Stokes equations** with Comsol's LiveLink for MATLAB.

## Multi-physics Coupling 🔄

Two strategies are implemented to tackle **multi-physics** problems:

- **Volume-based coupling**: Fully integrated physics within a single formulation (e.g., `MagnetohydrodynamicsCURL_HDG`, `MagnetohydrodynamicsCURLCURL_HDG`, `Plasma1FluidElectromagneticAdvanced_HDG`).

- **Surface-based coupling**: Interaction between sub-problems exchanging information at interfaces. For example, fluid-structure interaction problems couple [^2]:

  - A **fluid** formulation (`WeaklyCompressibleFlowDM_HDG` or `WeaklyCompressibleFlowVP_HDG`)

  - A **structure** formulation (`Elasticity_CG`)

  - A **moving mesh** formulation (`Elasticity_CG`)

A key feature of ***galerkin*** is the ability to **monolithically couple** $N$ different formulations, by automatically assembling their left-hand side matrix (LHS) and right-hand side vector (RHS) at each Newton iteration as:
```math
\begin{bmatrix}
\mathbf{K}_{11} & \dots  & \mathbf{K}_{1N} \\
\vdots          & \ddots & \vdots          \\
\mathbf{K}_{N1} & \dots  & \mathbf{K}_{NN}
\end{bmatrix}
\begin{bmatrix}
\Delta \mathbf{u}_{1} \\
\vdots                \\
\Delta \mathbf{u}_{N}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{f}_{1} \\
\vdots         \\
\mathbf{f}_{N}
\end{bmatrix}
```

## Project Structure 📁

The project is organized into the following directories:

- `formulations/`: Implemented finite element formulations

- `functions/`: Auxiliary functions

- `geometry/`: Geometries, meshes, and generation scripts

- `input/`: Simulation input files

- `output/`: Simulation results

- `tests/`: Test cases for validation and verification

## Simulation Types 📌

The following simulation types are supported:

- `SingleSimulation`: Runs a **single simulation**.

- `ConvergenceSpace`: Conducts a **spatial convergence** study.

- `ConvergenceTime`: Conducts a **temporal convergence** study.

- `ParametricStudy`: Conducts a **parametric study** by varying a single parameter.

- `ScalingStrong`: Assesses **strong scaling** performance.

- `ScalingWeak`: Assesses **weak scaling** performance.

## Data Types 📊

The data are essentially organized in structure arrays.

The structs defined in the **input file** are:

- `Simulation`: Defines the simulation type and the physical problem, including optional settings for restart, parallel execution, and multiprecision computing.

- `Parameters(iD)`: Specifies the main discretization-dependent parameters, such as the formulation, the analytical solution, the polynomial degree, the physical coefficients, and formulation-related parameters. If $N$ discretizations are coupled, the parameters must be set for each ($\text{iD} = 1,\dots,N$).

- `Geometry(iD)`: Stores the geometries as `DiscreteGeometry` objects, leveraging MATLAB's [Partial Differential Equation Toolbox](https://www.mathworks.com/help/pde/index.html).

- `Mesh(iD)`: Contains the mesh data, including nodes' coordinates and the elements' connectivity.

- `System`: Holds the global LHS matrix and RHS vector, along with user-defined settings such as tolerance and maximum iterations.

- `Time`: Stores time-related data, including initial, current, and final times, time step size, and the selected BDF order.

- `Solver`: Contains the chosen solver and preconditioner along with solver-specific parameters.

- `Boundaries(iD)`: Specifies the boundary splitting (edges in 2D, faces in 3D) for imposing boundary conditions.

- `Options`: Enables various features, such as plotting the geometry (`Options.PlotGeometry='yes'`), mesh and boundary conditions (`Options.PlotMesh='yes'`), or solution (`Options.PlotSolution={'VariableName'}`). Additional options include error computation (`Options.ComputeError={'VariableName'}`) and result saving/exporting.

The structs created in the **main file** are:

- `BCs(iD)`: Stores the nodes (relative to the underlying linear mesh) for boundary condition visualization.

- `Block`: Contains block-dependent data, including LHS, RHS, and matrix assembly indices.

- `Elements`: Organizes data in an element-wise manner for optimized `parfor` execution.

- `Faces(iD1,iD2)`: Stores face information for each sub-problem and interface coupling between different discretizations.

- `Memory`: Tracks memory usage during different simulation phases.

- `RefElement(iD1,iD2)`: Stores reference element data for each discretization and their coupling, including node coordinates, Gauss points, weights, and shape functions.

- `Results(iD)`: Stores results at chosen time steps.

- `Sizes(iD)`: Stores discretization-dependent sizes, such as global/local component counts, number of elements, faces, and nodes.

- `Timer`: Records time spent on pre-processing, processing (evaluation, solution, local problems), and post-processing tasks.

## Formulation Class 🧩

To facilitate the implementation of new finite element methods, ***galerkin*** allows you to declare a new class (likely the only part of the code you need to touch) that inherits from the base `Formulation` class.

Its main **properties** are:

- `NumGlobalComp`: Number of components for the global problem

- `NumLocalComp`: Number of components for the local problems (for HDG and FCFV)

- `NumPostComp`: Number of components of the post-processed variable(s) (for HDG)

- `DiscretizationType`: Discretization type (CG, HDG, or FCFV)

- `TimeDerOrder`: Time derivative order (e.g., 1 for thermal problems and 2 for elastic problems)

- `Domain`: Problem domain (time or frequency)

Its main **methods** are:

-  `initializeUnknowns()`: Initializes the problem unknowns.

-  `computeInitialConditions()`: Computes the initial conditions.

- `buildBlock()`: Constructs the LHS (diagonal block) and RHS of the sub-problem as
```math
\begin{bmatrix}
      & \vdots          &       \\
\dots & \mathbf{K}_{ii} & \dots \\
      & \vdots          &
\end{bmatrix}
\text{\quad and \quad}
\begin{bmatrix}
\vdots         \\
\mathbf{f}_{i} \\
\vdots
\end{bmatrix}
```

- `doCoupling()`: Constructs the LHS coupling block (off-diagonal block) as
```math
\begin{bmatrix}
      & \vdots          &       \\
\dots & \mathbf{K}_{ij} & \dots \\
      & \vdots          &
\end{bmatrix}
```

- `doPostProcess()`: Constructs the LHS and RHS of the post-processing problem.

- `storeResults()`: Stores results for analysis and visualization.

## Numerical Methods 🔢

Here, a glimpse of the main numerical methods behind ***galerkin*** is given.

### Spatial Discretization 📐

- **Element type**: Triangles (2D) and Tetrahedra (3D)

- **Polynomial degree**: From order 1 to 8

- **Nodes distribution**: Uniform (equally spaced) and Fekete (for better numerical conditioning)

### Time Integration ⏳

- **Backward Differentiation Formulas** (BDF): From order 1 to 6. BDF2 is initialized with the backward Euler method (BDF1) at the first time step. Schemes with $\text{BDFOrder}$ > 2 are initialized (the analytical solution needs to be available!) with the solution at the times $t = −n\Delta t$ with $n = [1,2,\dots,\text{BDFOrder}−1]$

- **Predictor**: Arbitrary order (0 for constant prediction, 1 for linear extrapolation, etc.)

### Solvers & Preconditioners 🧮

***galerkin*** allows the selection of various MATLAB built-in **solvers** and **preconditioners** to efficiently solve the resulting linear system:

- **Solvers**:

	- [`backslash`](https://www.mathworks.com/help/matlab/ref/double.mldivide.html): MATLAB’s popular **\ operator** for general linear systems of equations

	- [`pcg`](https://www.mathworks.com/help/matlab/ref/pcg.html?s_tid=doc_ta): **Preconditioned conjugate gradient** method for symmetric positive definite matrices

	- [`minres`](https://www.mathworks.com/help/matlab/ref/minres.html?s_tid=doc_ta): **Minimum residual** method for symmetric and non-positive definite matrices

	- [`gmres`](https://www.mathworks.com/help/matlab/ref/gmres.html?s_tid=doc_ta): **Generalized minimum residual** method for non-symmetric and non-positive definite matrices

- **Preconditioners**:

	- [`ichol`](https://www.mathworks.com/help/matlab/ref/ichol.html?s_tid=doc_ta): **Incomplete Cholesky** factorization

	- [`ilu`](https://www.mathworks.com/help/matlab/ref/ilu.html?s_tid=doc_ta): **Incomplete LU** factorization

### Error Norms 📈

When developing novel finite element formulations, it is crucial to perform a comprehensive error assessment to validate **model accuracy**. The framework supports multiple **error norms**:

- `Number` (i.e., scalar) norm computed in relative terms as
```math
\|E\| = \left|\dfrac{u^h-u^\text{ref}}{u^\text{ref}}\right|
```

- `L2` norm computed as
```math
\|E\|_{L^2(\Omega)} = \sqrt{\sum_{e=1}^{n_\text{el}}\int_{\Omega_e} \|\boldsymbol{u}^h-\boldsymbol{u}^\text{ref}\|^2}
```

- `Hdiv` norm computed as
```math
\|E\|_{H^\text{div}(\Omega)} = \sqrt{\sum_{e=1}^{n_\text{el}}\int_{\Omega_e} \left(\|\boldsymbol{u}^h-\boldsymbol{u}^\text{ref}\|^2+\|\boldsymbol{\nabla}\cdot\boldsymbol{u}^h-\boldsymbol{\nabla}\cdot\boldsymbol{u}^\text{ref}\|^2\right)}
```

- `Hcurl` norm computed as
```math
\|E\|_{H^\text{curl}(\Omega)} = \sqrt{\sum_{e=1}^{n_\text{el}}\int_{\Omega_e} \left(\|\boldsymbol{u}^h-\boldsymbol{u}^\text{ref}\|^2+\|\boldsymbol{\nabla}\times\boldsymbol{u}^h-\boldsymbol{\nabla}\times\boldsymbol{u}^\text{ref}\|^2\right)}
```

## Parallel Computing ⚡

***galerkin*** supports (limited) parallel computing capabilities using MATLAB's `parfor` for element computations. It achieves up to ***~10x speedup*** on relatively large-scale problems.

For cluster execution, modify the `test.sh` script in the main folder.

## External Tools 🛠️

Although ***galerkin*** can manage the main pre-processing, processing, and post-processing tasks, it can also leverage the functionalities of some external tools:

- **Free** tools:

	- [**GMSH**](https://gmsh.info): An open-source software, ideal for **high-order mesh generation**

	- [**ParaView**](https://www.paraview.org): The world’s leading open-source **post-processing visualization engine**
 
	- [**DistMesh**](https://doi.org/10.1137/S0036144503429121): A simple yet powerful unstructured **mesh generator** for MATLAB

- **Licensed** software:

	- [**Advanpix**](https://www.advanpix.com): A **multiprecision computing toolbox** for MATLAB

	- [**Comsol Multiphysics**](https://www.comsol.com): A powerful and versatile **simulation software**

## Getting Started 🚀

Clone the **repository** by running:
```bash
git clone https://github.com/andrealaspina/galerkin.git
```

From the main folder, execute the following command in MATLAB's **Command Window**:
```matlab
run test.m
```
This command executes all scripts in the `tests/` folder in ~3 minutes.

If all tests pass successfully, start by exploring `main.m`, which is structured into three sections: **Pre-processing**, **Processing**, and **Post-processing**. Although lengthy (~1500 lines), this file helps track the sequence of tasks from start to end of the simulation.

Next, check out some **test files**, which are named descriptively for easy navigation.

For a simple one-way coupled multi-physics problem, check out the **thermo-mechanical** problem in `input/jet_engine_blade`, reproducing the [`Thermal Stress Analysis of Jet Engine Turbine Blade`](https://www.mathworks.com/help/pde/ug/thermal-stress-analysis-of-jet-engine-turbine-blade.html) example in MATLAB's Partial Differential Equation Toolbox:

<p align="center">
  <img src="https://github.com/user-attachments/files/19511914/jet_engine_blade.pdf" width="500">
</p>

For a more advanced problem, refer to the **fluid-structure interaction** benchmark in `input/fsi_benchmark*`. By refining the mesh and decreasing the time step, you should be able to reproduce the results from my PhD thesis[^3]:

<p align="center">
  <img src="https://github.com/user-attachments/files/18719601/fsi.pdf" width="500">
</p>

For a Machiavellian multi-physics example, check out the **diocotron instability** problem in `input/diocotron_instability`. Similarly, you can reproduce the results from my paper[^7]:

<p align="center">
  <img src="https://github.com/user-attachments/assets/6b5da1c5-2d57-48a8-9f96-79411d220962" width="500">
</p>

## Contact 📧

***galerkin*** is solo-developed and maintained by **Dr.-Ing. Andrea La Spina**.

For any questions, feedback, or inquiries, contact me at andrealaspina.galerkin@gmail.com.

## My Mantra 💬

_If it doesn’t converge to machine precision, it’s simply wrong_.

## References 📖

[^1]: **A. La Spina**, M. Giacomini, A. Huerta, [_Hybrid coupling of CG and HDG discretizations based on Nitsche's method_](https://doi.org/10.1007/s00466-019-01770-8), Comput. Mech. (2020).

[^2]: **A. La Spina**, M. Kronbichler, M. Giacomini, W. A. Wall, A. Huerta, [_A weakly compressible hybridizable discontinuous Galerkin formulation for fluid-structure interaction problems_](https://doi.org/10.1016/j.cma.2020.113392), CMAME (2020).

[^3]: **A. La Spina**, [_Coupling of continuous and hybridizable discontinuous Galerkin methods for weakly compressible fluid-structure interaction_](https://mediatum.ub.tum.de/1586346), PhD thesis (2021).

[^4]: **A. La Spina**, J. Fish, [_A superconvergent hybridizable discontinuous Galerkin method for weakly compressible magnetohydrodynamics_](https://doi.org/10.1016/j.cma.2021.114278), CMAME (2022).

[^5]: **A. La Spina**, J. Fish, [_Time- and frequency-domain hybridizable discontinuous Galerkin solvers for the calculation of the Cherenkov radiation_](https://doi.org/10.1016/j.cma.2022.115170), CMAME (2022).

[^6]: J. Cui, **A. La Spina**, J. Fish, [_Data-physics driven multiscale approach for high-pressure resin transfer molding (HP-RTM)_](https://doi.org/10.1016/j.cma.2023.116405), CMAME (2023).

[^7]: **A. La Spina**, J. Fish, [_A hybridizable discontinuous Galerkin formulation for the Euler-Maxwell plasma model_](https://doi.org/10.1016/j.jcp.2023.112535), JCP (2024).

[^8]: M. Giacomini, D. Cortellessa, L. M. Vieira, R. Sevilla, A. Huerta, [_A hybrid pressure formulation of the face-centred finite volume method for viscous laminar incompressible flows_](https://doi.org/10.1002/nme.70037), IJNME (2025).