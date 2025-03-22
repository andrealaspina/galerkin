close all; clear; clc;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='FluidStructureInteraction';  % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='WeaklyCompressibleFlowDM_HDG';% Formulation
Parameters(1).Problem='Fluid';                   % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).ConvectiveFlow='yes';              % Convective flow
Parameters(1).ArbitraryLagrangianEulerian='yes'; % Arbitrary Lagrangian-Eulerian description
Parameters(1).Degree=4;                          % Degree
Parameters(1).StabDensity=100/1e-5;              % Stabilization for density
Parameters(1).StabMomentum=10;                   % Stabilization for momentum
Parameters(1).ReferenceDensity=1e3;              % Reference density
Parameters(1).ReferencePressure=0;               % Reference pressure
Parameters(1).CompressibilityCoeff=1e-5;         % Compressibility coefficient
Parameters(1).DynamicViscosity=1;                % Dynamic viscosity
Parameters(1).ScaledStrainRate=...               % Scaled strain rate
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(1).Density=@(x,y,z,t) 1e3*(x==x);     % Density
Parameters(1).Momentum=...                       % Momentum
  @(x,y,z,t) [(abs(x)<1e-3)*3/2.*y.*(0.41-y)/(0.41/2)^2*1e3*1*((1-cos(pi*t/2))/2*(t<=2)+(t>2)),...
              0*x];
Parameters(1).Traction=@(x,y,z,t) [0*x, 0*x];    % Traction
Parameters(1).ResidualContinuity=@(x,y,z,t) 0*x; % Residual of continuity equation
Parameters(1).Force=@(x,y,z,t) [0*x,0*x];        % Force
Parameters(1).Pressure=@(Density)...             % Equation of state
  Parameters(1).ReferencePressure+...
 (Density-Parameters(1).ReferenceDensity)/Parameters(1).CompressibilityCoeff;
Parameters(1).DPressureDDensity=@(Density)...    % dPressure/dDensity
  1/Parameters(1).CompressibilityCoeff;

Parameters(2).Formulation='Elasticity_CG';       % Formulation
Parameters(2).Problem='Mesh';                    % Problem
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Assumption='PlaneStrain';          % Assumption
Parameters(2).Degree=2;                          % Degree
Parameters(2).NitschePenalty=1e5;                % Nitsche's penalty parameter
Parameters(2).Density=0;                         % Density
Parameters(2).YoungsModulus=...                  % Young's modulus
  @(x,y,z) 1+(5-1)*max((1-abs(y-0.2)/0.2).^2,0)+(20-1)*((x-0.6).^2+(y-0.2).^2<0.03^2);
Parameters(2).PoissonsRatio=@(x,y,z) 0;          % Poisson's ratio
Parameters(2).Displacement=@(x,y,z,t) [0*x,0*x]; % Displacement
Parameters(2).Traction=...                       % Traction
  @(x,y,z,t,nx,ny,nz) [0*x,0*x];
Parameters(2).Force=@(x,y,z,t) [0*x,0*x];        % Force

Parameters(3).Formulation='Elasticity_CG';       % Formulation
Parameters(3).Problem='Structural';              % Problem
Parameters(3).Model='StVenantKirchhoff';         % Model
Parameters(3).Assumption='PlaneStrain';          % Assumption
Parameters(3).Degree=4;                          % Degree
Parameters(3).NitschePenalty=1e7;                % Nitsche's penalty parameter
Parameters(3).Density=10e3;                      % Density
Parameters(3).YoungsModulus=@(x,y,z) 1.4e6;      % Young's modulus
Parameters(3).PoissonsRatio=@(x,y,z) 0.4;        % Poisson's ratio
Parameters(3).Displacement=@(x,y,z,t) [0*x,0*x]; % Displacement
Parameters(3).Traction=...                       % Traction
  @(x,y,z,t,nx,ny,nz) [0*x,0*x];
Parameters(3).Force=@(x,y,z,t) [0*x,0*x];        % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_turek_fsi_unstructured'};       % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-2;                           % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=15;                               % Final time
Time.TimeStepSize=0.01;                          % Time step size
Time.BDFOrder=2;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries(1).Dirichlet_r=15;                    % Dirichlet portion for density
Boundaries(1).Dirichlet_w_x=[1:9,11,14];         % Dirichlet portion for x-momentum
Boundaries(1).Dirichlet_w_y=[1:9,11,14];         % Dirichlet portion for y-momentum
Boundaries(1).Interface=[10,12,13];              % Interface portion
Boundaries(1).Neumann_t_x=[];                    % Neumann portion for x-traction
Boundaries(1).Neumann_t_y=[];                    % Neumann portion for y-traction
Boundaries(2).Dirichlet=[1:9,11,14,15];          % Dirichlet portion
Boundaries(2).Interface=[10,12,13];              % Interface portion
Boundaries(2).Neumann=[];                        % Neumann portion
Boundaries(3).Dirichlet=[];                      % Dirichlet portion
Boundaries(3).Interface=[2,3,4];                 % Interface portion
Boundaries(3).Neumann=[];                        % Neumann portion
Boundaries(3).Fixed=1;                           % Fixed portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='yes';                      % Plot geometry
Options.PlotMesh='yes';                          % Plot mesh
Options.PlotMeshDistortion='yes';                % Plot mesh distortion
Options.PlotSolution=...                         % Plot solution
  {'ScaledStrainRate';
   'Density';
   'Momentum';
   'Displacement'};
Options.ComputeQuantity=...                      % Compute quantity of interest
  ['[~,Node]=ismembertol([0.6,0.2,0.0],Mesh(3).Nodes'',1e-6,''ByRows'',true);',...
   'Results(3).TipDisplacement(Time.TimeStep,:)=Block(3,3).SolutionGlobal(Node,:);',...
   'fprintf(''\n\nTipDisplacement = [%.4f,%.4f]*1e-3\n'',',...
   'Results(3).TipDisplacement(Time.TimeStep,1)*1e3,',...
   'Results(3).TipDisplacement(Time.TimeStep,2)*1e3)'];
Options.ComputeQuantityEnd=...                   % Compute quantity of interest (end of simulation)
  ['plot(Time.TimeStepSize:Time.TimeStepSize:Time.FinalTime,Results(3).TipDisplacement(:,2));',...
   'xlabel(''Time [s]'');',...
   'ylabel(''Tip vertical displacement [m]'');'];
Options.Export2Paraview=...                      % Export to Paraview
  {'ScaledStrainRate';
   'Density';
   'Momentum';
   'Displacement'};
Options.Export2ParaviewTimeSteps='Time.TimeStep';% Export to Paraview time steps
Options.StoreTimeSteps='Time.TimeStep';          % Store time steps
Options.SaveResults='yes';                       % Save results
Options.SaveResultsFolder='fsi_benchmark_dm/';   % Save results folder
% --------------------------------------------------------------------------------------------------

%% Main

main