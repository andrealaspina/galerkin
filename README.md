

## TODO

- Add a picture
- For more info on the formulations, refer to the associated .mlx files or my publications in the references
- Geometry and boundaries leveraging Matlab PDE toolbox
- Linear mesh increased in order or directly high-order mesh
- Useful tests to check with a short explanation for each of them
- Description of the main structs: `Simulation`, `Parameters`, `Geometry`, `Mesh`, (`MeshFile`), `System`, `Time`, `Solver`, `Boundaries`, `Options` + `BCs`, `Block`, `Elements`, `Faces`, `RefElement`, `Sizes`, `Timer`
- Description of the most important files: `test.m`, `main.m`
- Elements computations (core of each formulation) in `computeBlockElement()`

# Galerkin: A Flexible Finite Element Framework

## Overview 🌍

Galerkin is a powerful and extensible [**MATLAB**](https://www.mathworks.com/products/matlab.html) framework designed for developing and testing advanced **finite element** formulations. It provides an integrated environment for simulating a variety of **2D/3D** **linear/nonlinear** **single/multi-physics** **single/multi-scale** problems in the **time/frequency** domain.

At the beginning of my PhD, I struggled with the complexity of a large-scale research code, and I faced major obstacles in implementing complex finite element formulations for fluid-structure interaction. This motivated me to create an **accessible** and **flexible** framework to facilitate rapid **development**, **testing**, and **deployment** of novel numerical methods. This should be particularly beneficial for both early-career and experienced **researchers** striving to advance the field of **computational science and engineering**.

## Features ✨

Galerkin comes with a range of features designed to streamline the finite element simulation process:

- Automated workflow covering the main **pre-processing**, **processing**, and **post-processing** tasks.

- Support for multiple finite element **discretizations** (feel free to implement more!), including:

	- **Continuous Galerkin** (**CG**)

	- **Hybridizable Discontinuous Galerkin** (**HDG**)

	- Coupling between **CG** and **HDG**

- Broad applicability to various **physical models**.

- **Coupling** of an arbitrary number of sub-problems.

- Strong focus on **high-order** methods.

- (Limited) **parallel computing** capabilities leveraging MATLAB’s `parfor`.

- Seamless integration with **external tools** like [**GMSH**](https://gmsh.info), [**ParaView**](https://www.paraview.org), [**Comsol Multiphysics**](https://www.comsol.com), and [**Advanpix**](https://www.advanpix.com).

## Project Structure 📁

The project is organized into the following directories:

- `advanpix`: Contains a comprehensive library of routines for arbitrary precision computations using [**Advanpix**](https://www.advanpix.com) (you would need to purchase the product to use its functionalities).

- `formulation`: Collects the implemented finite element formulations.

- `functions`: Collects a variety of auxiliary functions.

- `geometry`: Stores all geometries and meshes and the scripts to generate them.

- `input`: Stores the simulation input files.

- `output`: Stores the simulation output files.

- `symbolic`: Collects symbolic computations for specific applications.

- `tests`: Collects all test cases for validation and verification.

## Implemented Formulations 📚

Galerkin includes implementations for a variety of physics-based formulations. Here is a (not complete!) list:

### Thermal Problems 🌡️

- `Thermal_CG`: Solves the **heat equation** with the CG method [^1].

- `Thermal_HDG`: Solves the **heat equation** with the HDG method [^1].

### Structural Mechanics 🔩

- `Elasticity_CG`: Solves the **equations of elastodynamics** with the CG method and covers both linear and nonlinear elasticity models (St. Venant-Kirchhoff and Neo-Hookean) [^1].

- `ElasticityLinear_HDG`: Solves the **equations of linear elastodynamics** with the HDG method [^1].

### Fluid Dynamics 🌊

- `CompressibleFlow_HDG`: Solves the **compressible Euler equations** with the HDG method.

- `WeaklyCompressibleFlowDM_HDG`: Solves the **Navier-Stokes equations** for weakly compressible flows with the HDG method using a density-momentum-based formulation [^3].

- `WeaklyCompressibleFlowVP_HDG`: Solves the **Navier-Stokes equations** for weakly compressible flows with the HDG method using a velocity-pressure-based formulation [^3].

### Electromagnetics 🧲

- `Electromagnetic_HDG_fast`: Solves the **Maxwell equations** in the time domain with the HDG method [^5].

- `ElectromagneticPML_HDG_fast`: Solves the **Maxwell equations** in the time domain with the HDG method including perfectly matched layers [^5].

- `ElectromagneticFD_HDG_fast`: Solves the **Maxwell equations** in the frequency domain with the HDG method [^5].

- `ElectromagneticPML_FD_HDG_fast`: Solves the **Maxwell equations** in the frequency domain with the HDG method including perfectly matched layers [^5].

- `Magnetic_HDG`: Solves the **magnetic induction equation** with the HDG method [^4].

- `MagneticCURLCURL_HDG`: Solves the **magnetic induction equation** with the HDG method using an alternative curl-curl formulation [^4].

### Plasma Physics 🔥

- `MagnetohydrodynamicsCURL_HDG`: Solves the **equations of magnetohydrodynamics** with the HDG method [^4].

- `MagnetohydrodynamicsCURLCURL_HDG`: Solves the **equations of magnetohydrodynamics** with the HDG method using an alternative curl-curl formulation [^4].

- `Plasma1FluidElectromagneticAdvanced_HDG`: Solves the **Euler–Maxwell plasma equations** with the HDG method including an electrostatic initialization and a projection-based divergence correction method to enforce the Gauss laws [^7].

### Porous Media 🧽

- `Darcy_CG`: Solves the **Darcy law** with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the **Navier-Stokes equations** with Comsol's LiveLink for Matlab [^6].

- `Darcy2Phase_CG`: Solves the **two-phase Darcy law** with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the **Cahn-Hilliard-Navier-Stokes equations** with Comsol's LiveLink for Matlab.

- `Darcy2PhaseRichards_CG`: Solves the **Richards’ equation** with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the **Cahn-Hilliard-Navier-Stokes equations** with Comsol's LiveLink for Matlab.

- `Darcy2PhaseRichards_HDG`: Solves the **Richards’ equation** with the HDG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the **Cahn-Hilliard-Navier-Stokes equations** with Comsol's LiveLink for Matlab.

## Multiphysics Coupling 🔄

Two strategies are implemented to facilitate multiphysics simulations:

- **Volume-based coupling**: Fully integrated physics within a single formulation (e.g., `MagnetohydrodynamicsCURL_HDG`, `MagnetohydrodynamicsCURLCURL_HDG`, `Plasma1FluidElectromagneticAdvanced_HDG`, etc.).

- **Surface-based coupling**: Interaction between sub-problems exchanging information on the interface. For instance, fluid-structure interaction problems are solved by coupling [^2]:

  - A formulation for the **fluid** sub-problem (`WeaklyCompressibleFlowDM_HDG` or `WeaklyCompressibleFlowVP_HDG`).

  - A formulation for the **structure** sub-problem (`Elasticity_CG`).

  - A formulation for the **moving mesh** algorithm (`Elasticity_CG`).

A key feature of Galerkin is the simplicity of **monolithically coupling** N different formulations, whose associated LHS and RHS are automatically combined at each Newton iteration as
```math
\begin{bmatrix}
\mathbf{K}_{11} & \dots  & \mathbf{K}_{1N} \\
\vdots          & \ddots & \vdots          \\
\mathbf{K}_{N1} & \dots  & \mathbf{K}_{NN}
\end{bmatrix}^{i}
\begin{bmatrix}
\Delta \mathbf{u}_{1} \\
\vdots                \\
\Delta \mathbf{u}_{N}
\end{bmatrix}^{i+1}
=
\begin{bmatrix}
\mathbf{f}_{1} \\
\vdots         \\
\mathbf{f}_{N}
\end{bmatrix}^{i}
```
## Simulation Types 📌

A variety of simulation types are readily supported:

- `SingleSimulation`: Runs a **single simulation**.

- `ConvergenceSpace`: Conducts a **spatial convergence** study.

- `ConvergenceTime`: Conducts a **temporal convergence** study.

- `ParametricStudy`: Conducts a **parametric study** by varying a parameter.

- `ScalingStrong`: Assesses **strong scaling** performance.

- `ScalingWeak`: Assesses **weak scaling** performance.

## Numerical Methods 🧮

### Spatial Discretization

- Element types: Triangles in 2D & Tetrahedra in 3D.

- Polynomial Degree: Up to 8th order.

- Node distributions: Uniform & Fekete.

### Time Integration

- **Backward Differentiation Formulas** (BDF): Supports up to 6th order.

- Predictor: Allows for high-order predictor schemes.

### Solvers & Preconditioners

Solvers:

- [`backslash`](https://www.mathworks.com/help/matlab/ref/double.mldivide.html): MATLAB’s **\** operator to solve general linear systems of equations.

- [`pcg`](https://www.mathworks.com/help/matlab/ref/pcg.html?s_tid=doc_ta): **Preconditioned conjugate gradient** method for symmetric positive definite matrices.

- [`minres`](https://www.mathworks.com/help/matlab/ref/minres.html?s_tid=doc_ta): **Minimum residual** method for symmetric and non-positive definite matrices.

- [`gmres`](https://www.mathworks.com/help/matlab/ref/gmres.html?s_tid=doc_ta): **Generalized minimum residual** method for non-symmetric and non-positive definite matrices.

Preconditioners:

- [`ichol`](https://www.mathworks.com/help/matlab/ref/ichol.html?s_tid=doc_ta): **Incomplete Cholesky** factorization.

- [`ilu`](https://www.mathworks.com/help/matlab/ref/ilu.html?s_tid=doc_ta): **Incomplete LU** factorization.

## Parallel Computing ⚡

Parallelized element computations utilizing `parfor`.

Up to **~10x performance improvement** on relatively large-scale problems.

Cluster execution support via `test.sh` for batch processing.

## Error Norms 📈

The framework supports multiple error norms:

- `Number` (i.e. scalar norm) computed in relative terms as
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

## External Tools 🛠️

Although Galerkin is able to manage the main*pre-processing, processing, and post-processing tasks, it can also leverage the functionalities of some external software:

- **Free** tools:

	- [**GMSH**](https://gmsh.info): An open-source software, ideal for **high-order mesh generation**.

	- [**ParaView**](https://www.paraview.org): The world’s leading open source **post-processing visualization engine**.

- **Paid** tools:

	- [**Advanpix**](https://www.advanpix.com): A **multiprecision computing toolbox** for MATLAB (supported for thermal problems only, but easily extendable).

	- [**Comsol Multiphysics**](https://www.comsol.com): A powerful and versatile **simulation software** (used for the solution of the unit cell problems for multi-scale computations).

## Getting Started 🚀

TODO

## Contributions & Feedback 🙌🏼

I encourage contributions and feedback! If this project proves helpful, consider starring ⭐ the repository, reporting issues, or submitting pull requests.

Happy coding! 👨🏻‍💻

## References

[^1]: **A. La Spina**, M. Giacomini, A. Huerta, [_Hybrid coupling of CG and HDG discretizations based on Nitsche's method_](https://doi.org/10.1007/s00466-019-01770-8), Comput. Mech. (2020).
[^2]: **A. La Spina**, M. Kronbichler, M. Giacomini, W. A. Wall, A. Huerta, [_A weakly compressible hybridizable discontinuous Galerkin formulation for fluid-structure interaction problems_](https://doi.org/10.1016/j.cma.2020.113392), CMAME (2020).
[^3]: **A. La Spina**, [_Coupling of continuous and hybridizable discontinuous Galerkin methods for weakly compressible fluid-structure interaction_](https://mediatum.ub.tum.de/1586346), Ph.D. thesis (2021).
[^4]: **A. La Spina**, J. Fish, [_A superconvergent hybridizable discontinuous Galerkin method for weakly compressible magnetohydrodynamics_](https://doi.org/10.1016/j.cma.2021.114278), CMAME (2022).
[^5]: **A. La Spina**, J. Fish, [_Time- and frequency-domain hybridizable discontinuous Galerkin solvers for the calculation of the Cherenkov radiation_](https://doi.org/10.1016/j.cma.2022.115170), CMAME (2022).
[^6]: J. Cui, **A. La Spina**, J. Fish, [_Data-physics driven multiscale approach for high-pressure resin transfer molding (HP-RTM)_](https://doi.org/10.1016/j.cma.2023.116405), CMAME (2023).
[^7]: **A. La Spina**, J. Fish, [_A hybridizable discontinuous Galerkin formulation for the Euler--Maxwell plasma model_](https://doi.org/10.1016/j.jcp.2023.112535), JCP (2024).






	
