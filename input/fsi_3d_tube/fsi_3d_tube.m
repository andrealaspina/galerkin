close all; clear; clc;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='FluidStructureInteraction';  % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='WeaklyCompressibleFlowVP_HDG';% Formulation
Parameters(1).Problem='Fluid';                   % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).ConvectiveFlow='yes';              % Convective flow
Parameters(1).ArbitraryLagrangianEulerian='yes'; % Arbitrary Lagrangian-Eulerian description
Parameters(1).Degree=4;                          % Degree
Parameters(1).StabVelocity=100;                  % Stabilization for velocity
Parameters(1).StabPressure=100;                  % Stabilization for pressure
Parameters(1).ReferenceDensity=1;                % Reference density
Parameters(1).ReferencePressure=0;               % Reference pressure
Parameters(1).CompressibilityCoeff=0;            % Compressibility coefficient
Parameters(1).DynamicViscosity=0.03;             % Dynamic viscosity
Parameters(1).ScaledStrainRate=...               % Scaled strain rate
  @(x,y,z,t) [0*x,0*x,0*x,0*x,0*x,0*x];
Parameters(1).Velocity=@(x,y,z,t) [0*x,0*x,0*x]; % Velocity
Parameters(1).Pressure=@(x,y,z,t) 0*x;           % Pressure
Parameters(1).Traction=...                       % Traction
  @(x,y,z,t) [(abs(x)<1e-3)*1.3332e4*(t<=3e-3),...
              0*x,...
              0*x];
Parameters(1).ResidualContinuity=@(x,y,z,t) 0*x; % Residual of continuity equation
Parameters(1).Force=@(x,y,z,t) [0*x,0*x,0*x];    % Force
Parameters(1).Density=@(Pressure)...             % Equation of state
  Parameters(1).ReferenceDensity+...
 (Pressure-Parameters(1).ReferencePressure)*Parameters(1).CompressibilityCoeff;
Parameters(1).DDensityDPressure=@(Pressure)...   % dDensity/dPressure
  Parameters(1).CompressibilityCoeff;
Parameters(1).Displacement=...                   % Displacement
  @(x,y,z,t) [0*x,0*x,0*x];

Parameters(2).Formulation='Elasticity_CG';       % Formulation
Parameters(2).Problem='Mesh';                    % Problem
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Degree=4;                          % Degree
Parameters(2).NitschePenalty=1e6;                % Nitsche's penalty parameter
Parameters(2).Density=0;                         % Density
Parameters(2).YoungsModulus=1;                   % Young's modulus
Parameters(2).PoissonsRatio=0;                   % Poisson's ratio
Parameters(2).Displacement=...                   % Displacement
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(2).Traction=...                       % Traction
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x,0*x];
Parameters(2).Force=@(x,y,z,t) [0*x,0*x,0*x];    % Force

Parameters(3).Formulation='Elasticity_CG';       % Formulation
Parameters(3).Problem='Structural';              % Problem
Parameters(3).Model='LinearElasticity';          % Model
Parameters(3).Degree=4;                          % Degree
Parameters(3).NitschePenalty=3e9;                % Nitsche's penalty parameter
Parameters(3).Density=1.2;                       % Density
Parameters(3).YoungsModulus=3e6;                 % Young's modulus
Parameters(3).PoissonsRatio=0.3;                 % Poisson's ratio
Parameters(3).Displacement=...                   % Displacement
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(3).Traction=...                       % Traction
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x,0*x];
Parameters(3).Force=@(x,y,z,t) [0*x,0*x,0*x];    % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_fsi_3d_tube'};                  % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
Solver.Equilibrate='yes';                        % Equilibrate matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-3;                           % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=0.015;                            % Final time
Time.TimeStepSize=0.015/40;                      % Time step size
Time.BDFOrder=2;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries(1).Dirichlet_v_x=3;                   % Dirichlet portion for x-velocity
Boundaries(1).Dirichlet_v_y=3;                   % Dirichlet portion for y-velocity
Boundaries(1).Dirichlet_v_z=3;                   % Dirichlet portion for z-velocity
Boundaries(1).Dirichlet_p=[];                    % Dirichlet portion for pressure
Boundaries(1).Interface=[2,4,5,6];               % Interface portion
Boundaries(1).Neumann_t_x=1;                     % Neumann portion for x-traction
Boundaries(1).Neumann_t_y=1;                     % Neumann portion for y-traction
Boundaries(1).Neumann_t_z=1;                     % Neumann portion for z-traction
Boundaries(2).Dirichlet=[1,3];                   % Dirichlet portion
Boundaries(2).Interface=[2,4,5,6];               % Interface portion
Boundaries(2).Neumann=[];                        % Neumann portion
Boundaries(3).Dirichlet=[];                      % Dirichlet portion
Boundaries(3).Interface=[3,6,8,10];              % Interface portion
Boundaries(3).Neumann=[1,5,7,9];                 % Neumann portion
Boundaries(3).Fixed=[2,4];                       % Fixed portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='yes';                      % Plot geometry
Options.PlotMesh='yes';                          % Plot mesh
Options.PlotMeshDistortion='yes';                % Plot mesh distortion
Options.Export2Paraview=...                      % Export to Paraview
  {'Velocity';
   'Pressure';
   'Displacement'};
Options.Export2ParaviewTimeSteps='Time.TimeStep';% Export to Paraview time steps
Options.StoreTimeSteps='Time.TimeStep';          % Store time steps
Options.SaveResults='yes';                       % Save results
Options.SaveResultsFolder='fsi_3d_tube/';        % Save results folder
% --------------------------------------------------------------------------------------------------

%% Main

main