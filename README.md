## TODO

- Project structure
- Refer to mxl files and Google Scholar
- Restart
- Useful tests to check

# Galerkin: A Flexible Finite Element Framework

## Overview 🌍

Galerkin is a powerful, adaptable, and extensible **MATLAB** framework designed for the development and testing of advanced **finite element** formulations. It provides an integrated environment for simulating **2D/3D** **linear/nonlinear** **single/multi-physics** **single/multi-scale** problems in the **time/frequency** domain.

During my PhD, I grappled with the overwhelming complexity of a large-scale C++-based research code, facing major obstacles in implementing intricate hybrid finite element formulations for fluid-structure interaction. These frustrations fueled my determination to develop a more **accessible** approach to finite element research. My objective is therefore to provide a streamlined framework that enables rapid **development**, **testing**, and **publication** of new methodologies. This should be particularly beneficial for both early-career and experienced **researchers** striving to advance the field of **computational science and engineering**.

## Features ✨

Comprehensive workflow automation covering main **pre-processing**, **processing**, and **post-processing** tasks.
Support for multiple finite element discretizations, including:

- **Continuous Galerkin** (**CG**)

- **Hybridizable Discontinuous Galerkin** (**HDG**)

- Coupling between CG and HDG

Broad applicability to various physics models.

Coupling of an arbitrary number of sub-problems.

Flexible discretization, solver selection, and high-order methods.

Parallel computing capabilities leveraging MATLAB’s `parfor` construct.

Seamless integration with external tools like GMSH, ParaView, Comsol, and Advanpix.

## Implemented Formulations 📚

The framework currently incorporates multiple formulations tailored for a variety of physics-based applications, making it a versatile tool for researchers and engineers. Here are listed the main (but not all!) formulations currently implemented.

### Thermal Problems 🌡️

- `Thermal_CG.m`: solves the heat equation with the CG method.

- `Thermal_HDG.m`: solves the heat equation with the HDG method.

### Structural Mechanics 🔩

- `Elasticity_CG.m`: solves the equations of elastodynamics with the CG method and covers both linear and nonlinear elasticity models (St. Venant-Kirchhoff and Neo-Hookean).

- `ElasticityLinear_HDG.m`: solves the equations of linear elastodynamics with the HDG method.

### Fluid Dynamics 🌊

- `CompressibleFlow_HDG.m`: solves the compressible Euler equations with the HDG method.

- `WeaklyCompressibleFlowDM_HDG.m`: solves the  Navier-Stokes equations for weakly compressible flows with the HDG method using a density-momentum-based formulation.

- `WeaklyCompressibleFlowVP_HDG.m`: solves the  Navier-Stokes equations for weakly compressible flows with the HDG method using a velocity-pressure-based formulation.

### Electromagnetics 🧲

- `Electromagnetic_HDG_fast.m`: solves the Maxwell equations in the time domain with the HDG method.

- `ElectromagneticPML_HDG_fast.m`: solves the Maxwell equations in the time domain with the HDG method including perfectly matched layers.

- `ElectromagneticFD_HDG_fast.m`: solves the Maxwell equations in the frequency domain with the HDG method.

- `ElectromagneticPML_FD_HDG_fast.m`: solves the Maxwell equations in the frequency domain with the HDG method including perfectly matched layers.

- `Magnetic_HDG.m`: solves the magnetic induction equation with the HDG method.

- `MagneticCURLCURL_HDG.m`: solves the magnetic induction equation with the HDG method and an alternative curl-curl formulation.

### Plasma Physics 🔥

- `MagnetohydrodynamicsCURL_HDG.m`: solves the equations of magnetohydrodynamics (MHD) with the HDG method.

- `MagnetohydrodynamicsCURLCURL_HDG.m`: solves the equations of magnetohydrodynamics (MHD) with the HDG method and an alternative curl-curl formulation.

- `Plasma1FluidElectromagneticAdvanced_HDG.m`: solves the Euler–Maxwell plasma equations with the HDG method including an electrostatic initialization and a projection-based divergence correction method to enforce the Gauss laws.

### Porous Media 🧽

- `Darcy_CG.m`: solves the  Darcy law with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Navier-Stokes equations with Comsol's LiveLink for Matlab.

- `Darcy2Phase_CG.m`: solves the two-phase Darcy law with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Cahn-Hilliard-Navier-Stokes equations with Comsol's LiveLink for Matlab.

- `Darcy2PhaseRichards_CG.m`: solves the Richards’ equation with the CG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Cahn-Hilliard-Navier-Stokes equations with Comsol's LiveLink for Matlab.

- `Darcy2PhaseRichards_HDG.m`: solves the Richards’ equation with the HDG method for the macroscopic problem, informed by the effective properties derived from a unit cell problem solving the Cahn-Hilliard-Navier-Stokes equations with Comsol's LiveLink for Matlab.

## Multiphysics Coupling 🔄

Two key coupling strategies are implemented to facilitate multiphysics simulations:

- Volume-based coupling: Fully integrated physics within a single formulation (e.g., MHD, Plasma).

- Surface-based coupling: Interaction between subsystems (e.g., Fluid-Structure Interaction).

For instance, fluid-structure interaction (FSI) problems are solved by coupling:

- HDG-based formulations for weakly compressible flow.

- CG-based formulations for (non)linear elasticity.

- CG-based formulations for the fluid mesh motion according to the Arbitrary Lagrangian-Eulerian (ALE) description of motion.

## Simulation Types 📌

A variety of simulation modes are supported:

- `SingleSimulation`: Runs a single test case.

- `ConvergenceSpace`: Conducts a spatial convergence study.

- `ConvergenceTime`: Executes a temporal convergence study.

- `ParametricStudy`: Evaluates different parametric variations.

- `ScalingStrong`: Analyzes strong scaling performance.

- `ScalingWeak`: Assesses weak scaling performance.

## Numerical Methods & Solvers 🧮

### Spatial Discretization

- Element types: Triangles & Tetrahedra.

- Polynomial orders: Up to 8th order.

- Node distributions: Uniform & Fekete points.

### Time Integration

- Backward Differentiation Formulas (BDF): Supports up to 6th order.

- Predictors: Allows for high-order predictor schemes.

### Solvers & Preconditioners

Solvers:

- `backslash`: MATLAB’s \ operator to solve general linear systems of equations.

- `pcg`: Preconditioned conjugate gradient method for symmetric positive definite matrices.

- `minres`: Minimum residual method for symmetric and non-positive definite matrices.

- `gmres`: Generalized minimum residual method for non-symmetric and non-positive definite matrices.

Preconditioners:

- `ichol`: Incomplete Cholesky factorization.

- `ilu`: Incomplete LU factorization.

## Parallel & High-Performance Computing ⚡

Parallelized element computations utilizing `parfor`.

Up to ~10x performance improvement on large-scale problems.

Cluster execution support via test.sh for batch processing.

## Error Analysis 📈

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

- **GMSH**: High-order mesh generation.

- **ParaView**: Advanced visualization.

Paid tools:

- **Comsol Multiphysics**: Used for multiscale modeling (e.g., unit cell problems).

- **Advanpix**: Arbitrary precision arithmetic (applied to thermal problems but extendable).

## Contributions & Feedback 🙌

I encourage contributions and feedback! If this project proves helpful, consider starring ⭐ the repository, reporting issues, or submitting pull requests.

Happy coding! 🚀
