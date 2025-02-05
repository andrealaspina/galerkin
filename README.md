## TODO

- Project structure
- Refer to mxl files and [**Google Scholar**](https://scholar.google.com/citations?user=5n7qnBgAAAAJ&hl=en)
- Restart
- Geometry and boundaries with PDE toolbox
- Linear mesh increased in order or directly high-order mesh
- Useful tests to check
- Structs: `Simulation`, `Parameters`, `Geometry`, `Mesh`, (`MeshFile`), `System`, `Time`, `Solver`, `Boundaries`, `Options` + `BCs`, `Block`, `Elements`, `Faces`, `RefElement`, `Sizes`, `Timer`
- Files: `test.m`, `main.m`
- Folders: (`advanpix`), `formulations`, `functions`, `geometry`, `input`, `output`, `symbolic`, `tests`

# Galerkin: A Flexible Finite Element Framework

## Overview 🌍

Galerkin is a powerful, adaptable, and extensible [**MATLAB**](https://www.mathworks.com/products/matlab.html) framework designed for the development and testing of advanced **finite element** formulations. It provides an integrated environment for simulating **2D/3D** **linear/nonlinear** **single/multi-physics** **single/multi-scale** problems in the **time/frequency** domain.

During my PhD, I grappled with the overwhelming complexity of a large-scale C++-based research code, facing major obstacles in implementing intricate hybrid finite element formulations for fluid-structure interaction. These frustrations fueled my determination to develop a more **accessible** approach to finite element research. My objective is therefore to provide a streamlined framework that enables rapid **development**, **testing**, and **publication** of new methodologies. This should be particularly beneficial for both early-career and experienced **researchers** striving to advance the field of **computational science and engineering**.

## Features ✨

Comprehensive workflow automation covering main **pre-processing**, **processing**, and **post-processing** tasks.
Support for multiple finite element discretizations, including:

- **Continuous Galerkin** (**CG**)

- **Hybridizable Discontinuous Galerkin** (**HDG**)

- Coupling between **CG** and **HDG**

Broad applicability to various physics models.

Coupling of an arbitrary number of sub-problems.

Flexible discretization, solver selection, and high-order methods.

Parallel computing capabilities leveraging MATLAB’s `parfor` construct.

Seamless integration with external tools like GMSH, ParaView, Comsol, and Advanpix.

## Implemented Formulations 📚

The framework currently incorporates multiple formulations tailored for a variety of physics-based applications, making it a versatile tool for researchers and engineers. Here are listed the main (but not all!) formulations currently implemented.

### Thermal Problems 🌡️

- `Thermal_CG`: solves the heat equation with the CG method (see [https://doi.org/10.1007/s00466-019-01770-8](https://doi.org/10.1007/s00466-019-01770-8)).

- `Thermal_HDG`: solves the heat equation with the HDG method (see [https://doi.org/10.1007/s00466-019-01770-8](https://doi.org/10.1007/s00466-019-01770-8)).

### Structural Mechanics 🔩

- `Elasticity_CG`: solves the equations of elastodynamics with the CG method and covers both linear and nonlinear elasticity models (St. Venant-Kirchhoff and Neo-Hookean) (see [https://doi.org/10.1007/s00466-019-01770-8](https://doi.org/10.1007/s00466-019-01770-8)).

- `ElasticityLinear_HDG`: solves the equations of linear elastodynamics with the HDG method (see [https://doi.org/10.1007/s00466-019-01770-8](https://doi.org/10.1007/s00466-019-01770-8)).

### Fluid Dynamics 🌊

- `CompressibleFlow_HDG`: solves the compressible Euler equations with the HDG method.

- `WeaklyCompressibleFlowDM_HDG`: solves the  Navier-Stokes equations for weakly compressible flows with the HDG method using a density-momentum-based formulation (see [https://mediatum.ub.tum.de/1586346](https://mediatum.ub.tum.de/1586346)).

- `WeaklyCompressibleFlowVP_HDG`: solves the  Navier-Stokes equations for weakly compressible flows with the HDG method using a velocity-pressure-based formulation (see [https://mediatum.ub.tum.de/1586346](https://mediatum.ub.tum.de/1586346)).

### Electromagnetics 🧲

- `Electromagnetic_HDG_fast`: solves the Maxwell equations in the time domain with the HDG method (see [https://doi.org/10.1016/j.cma.2022.115170](https://doi.org/10.1016/j.cma.2022.115170)).

- `ElectromagneticPML_HDG_fast`: solves the Maxwell equations in the time domain with the HDG method including perfectly matched layers (see [https://doi.org/10.1016/j.cma.2022.115170](https://doi.org/10.1016/j.cma.2022.115170)).

- `ElectromagneticFD_HDG_fast`: solves the Maxwell equations in the frequency domain with the HDG method (see [https://doi.org/10.1016/j.cma.2022.115170](https://doi.org/10.1016/j.cma.2022.115170)).

- `ElectromagneticPML_FD_HDG_fast`: solves the Maxwell equations in the frequency domain with the HDG method including perfectly matched layers (see [https://doi.org/10.1016/j.cma.2022.115170](https://doi.org/10.1016/j.cma.2022.115170)).

- `Magnetic_HDG`: solves the magnetic induction equation with the HDG method (see [https://doi.org/10.1016/j.cma.2021.114278](https://doi.org/10.1016/j.cma.2021.114278)).

- `MagneticCURLCURL_HDG`: solves the magnetic induction equation with the HDG method and an alternative curl-curl formulation (see [https://doi.org/10.1016/j.cma.2021.114278](https://doi.org/10.1016/j.cma.2021.114278)).

### Plasma Physics 🔥

- `MagnetohydrodynamicsCURL_HDG`: solves the equations of magnetohydrodynamics (MHD) with the HDG method (see [https://doi.org/10.1016/j.cma.2021.114278](https://doi.org/10.1016/j.cma.2021.114278)).

- `MagnetohydrodynamicsCURLCURL_HDG`: solves the equations of magnetohydrodynamics (MHD) with the HDG method and an alternative curl-curl formulation (see [https://doi.org/10.1016/j.cma.2021.114278](https://doi.org/10.1016/j.cma.2021.114278)).

- `Plasma1FluidElectromagneticAdvanced_HDG`: solves the Euler–Maxwell plasma equations with the HDG method including an electrostatic initialization and a projection-based divergence correction method to enforce the Gauss laws (see [https://doi.org/10.1016/j.jcp.2023.112535](https://doi.org/10.1016/j.jcp.2023.112535)).

### Porous Media 🧽

- `Darcy_CG`: solves the  Darcy law with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Navier-Stokes equations with Comsol's LiveLink for Matlab (see [https://doi.org/10.1016/j.cma.2023.116405](https://doi.org/10.1016/j.cma.2023.116405)).

- `Darcy2Phase_CG`: solves the two-phase Darcy law with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Cahn-Hilliard-Navier-Stokes equations with Comsol's LiveLink for Matlab.

- `Darcy2PhaseRichards_CG`: solves the Richards’ equation with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Cahn-Hilliard-Navier-Stokes equations with Comsol's LiveLink for Matlab.

- `Darcy2PhaseRichards_HDG`: solves the Richards’ equation with the HDG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Cahn-Hilliard-Navier-Stokes equations with Comsol's LiveLink for Matlab.

## Multiphysics Coupling 🔄

Two strategies are implemented to facilitate multiphysics simulations:

- Volume-based coupling: Fully integrated physics within a single formulation (e.g., `MagnetohydrodynamicsCURL_HDG`, `MagnetohydrodynamicsCURLCURL_HDG`, `Plasma1FluidElectromagneticAdvanced_HDG`).

- Surface-based coupling: Interaction between subsystems. For instance, fluid-structure interaction (FSI) problems (see [https://doi.org/10.1016/j.cma.2020.113392](https://doi.org/10.1016/j.cma.2020.113392)) are solved by coupling:

  - A formulation for the fluid sub-problem (`WeaklyCompressibleFlowDM_HDG` or `WeaklyCompressibleFlowVP_HDG`).

  - A formulation for the structural sub-problem (`Elasticity_CG`).

  - A formulation for the moving mesh algorithm (`Elasticity_CG`).

A key feature of **Galerkin** is the simplicity of monolithically coupling N different formulations, whose associated LHS and RHS are automatically combined as
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

A variety of simulation types are supported:

- `SingleSimulation`: Runs a single simulation.

- `ConvergenceSpace`: Conducts a spatial convergence study.

- `ConvergenceTime`: Conducts a temporal convergence study.

- `ParametricStudy`: Conducts a parametric study by varying a parameter.

- `ScalingStrong`: Assesses strong scaling performance.

- `ScalingWeak`: Assesses weak scaling performance.

## Numerical Methods 🧮

### Spatial Discretization

- Element types: Triangles in 2D & Tetrahedra in 3D.

- Polynomial orders: Up to 8th order.

- Node distributions: Uniform & Fekete.

### Time Integration

- **Backward Differentiation Formulas** (BDF): Supports up to 6th order.

- Predictor: Allows for high-order predictor schemes.

### Solvers & Preconditioners

Solvers:

- [`backslash`](https://www.mathworks.com/help/matlab/ref/double.mldivide.html): MATLAB’s \ operator to solve general linear systems of equations.

- [`pcg`](https://www.mathworks.com/help/matlab/ref/pcg.html?s_tid=doc_ta): Preconditioned conjugate gradient method for symmetric positive definite matrices.

- [`minres`](https://www.mathworks.com/help/matlab/ref/minres.html?s_tid=doc_ta): Minimum residual method for symmetric and non-positive definite matrices.

- [`gmres`](https://www.mathworks.com/help/matlab/ref/gmres.html?s_tid=doc_ta): Generalized minimum residual method for non-symmetric and non-positive definite matrices.

Preconditioners:

- [`ichol`](https://www.mathworks.com/help/matlab/ref/ichol.html?s_tid=doc_ta): Incomplete Cholesky factorization.

- [`ilu`](https://www.mathworks.com/help/matlab/ref/ilu.html?s_tid=doc_ta): Incomplete LU factorization.

## Parallel Computing ⚡

Parallelized element computations utilizing `parfor`.

Up to **~10x performance improvement** on relatively large-scale problems.

Cluster execution support via `test.sh` for batch processing.

## Error Norms 📈

The framework supports multiple error norms:

- `Number`, i.e. scalar norm, computed in relative terms as
```math
\|E\| = \left|\dfrac{u^h-u^\text{ref}}{u^\text{ref}}\right|
```

- `L2` norm computed as
```math
\|E\|_{L^2(\Omega)} = \sqrt{\sum_{e=1}^{n_\text{el}}\int_{\Omega^e} \|\boldsymbol{u}^h-\boldsymbol{u}^\text{ref}\|^2}
```

- `Hdiv` norm computed as
```math
\|E\|_{H^\text{div}(\Omega)} = \sqrt{\sum_{e=1}^{n_\text{el}}\int_{\Omega^e} \left(\|\boldsymbol{u}^h-\boldsymbol{u}^\text{ref}\|^2+\|\boldsymbol{\nabla}\cdot\boldsymbol{u}^h-\boldsymbol{\nabla}\cdot\boldsymbol{u}^\text{ref}\|^2\right)}
```

- `Hcurl` norm computed as
```math
\|E\|_{H^\text{curl}(\Omega)} = \sqrt{\sum_{e=1}^{n_\text{el}}\int_{\Omega^e} \left(\|\boldsymbol{u}^h-\boldsymbol{u}^\text{ref}\|^2+\|\boldsymbol{\nabla}\times\boldsymbol{u}^h-\boldsymbol{\nabla}\times\boldsymbol{u}^\text{ref}\|^2\right)}
```

## External Tools 🛠️

Free tools:

- [**GMSH**](https://gmsh.info): High-order mesh generation.

- [**ParaView**](https://www.paraview.org): Advanced visualization.

Paid tools:

- [**Comsol Multiphysics**](https://www.comsol.com): Used for multiscale modeling (e.g., unit cell problems).

- [**Advanpix**](https://www.advanpix.com): Arbitrary precision arithmetic (applied to thermal problems but extendable).

## Contributions & Feedback 🙌

I encourage contributions and feedback! If this project proves helpful, consider starring ⭐ the repository, reporting issues, or submitting pull requests.

Happy coding! 🚀
