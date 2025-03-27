close all; clear; clc;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Plasma1FluidElectromagnetic';% Problem
Simulation.NumProcessors=1;                      % Number of processors
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Plasma1FluidElectromagneticAdvanced_HDG';% Formulation
Parameters.Problem='Plasma1FluidElectromagnetic';% Problem
Parameters.Axisymmetric='no';                    % Axisymmetric
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.Electromagnetic2FluidCoupling='yes';  % Electromagnetic->Fluid coupling
Parameters.Fluid2ElectromagneticCoupling='yes';  % Fluid->Electromagnetic coupling
Parameters.SolveEnergyEquation='no';             % Solve energy equation
Parameters.EquilibrateLocalProblems='yes';       % Equilibrate local problems
Parameters.Degree=1; % =8 in the paper           % Degree
Parameters.MeshRefinement=1; % =3 in the paper   % Mesh refinement
Parameters.PredictorDegree=2;                    % Predictor degree
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.SpecificHeatRatio=1;                  % Ratio of specific heats
Parameters.MagneticPermeability=1e-6;            % Magnetic permeability
Parameters.ElectricPermittivity=1e-2;            % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.ParticleCharge=-1;                    % Particle charge
Parameters.IonCharge=+1;                         % Ion charge
Parameters.ParticleMass=1e-6/(1+1e-6);           % Particle mass
Parameters.IonMass=1/(1+1e-6);                   % Ion mass
Parameters.DampingFunction=@(x,y,z) [0*x,0*x];   % Damping function for PML
Parameters.Density=...                           % Density
  @(x,y,z,t) Parameters.ParticleMass*(1*1e-1*(sqrt(x.^2+y.^2)< 6)+...
                                      1*1e-0*(sqrt(x.^2+y.^2)>=6 & sqrt(x.^2+y.^2)<=8)+...
                                      1*1e-4*(sqrt(x.^2+y.^2)> 8)).*...
                                     (1+1e-3*cos(5*atan2(y,x)))*(t==0);
Parameters.IonDensity=(1/(1+1e-6))*(1*1e-4);     % Ion density
Parameters.Momentum=@(x,y,z,t) [0*x,0*x];        % Momentum
Parameters.Energy=@(x,y,z,t) 0*x;                % Energy
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,t) 5/Parameters.MagneticPermeability*(x==x);
Parameters.ElectricField=@(x,y,z,t) [0*x,0*x];   % Electric field
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier for hyperbolic div. cleaning
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.ForceDensity=@(x,y,z,t) 0*x;          % Force in density equation
Parameters.ForceMomentum=@(x,y,z,t) [0*x,0*x];   % Force in momentum equation
Parameters.ForceEnergy=@(x,y,z,t) 0*x;           % Force in energy equation
Parameters.CurrentDensity=@(x,y,z,t) [0*x,0*x];  % Current density
Parameters.IncidentField=@(x,y,z,t) [0*x,0*x];   % Incident field for absorbing boundary conditions
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={...                                   % Mesh file
  ['Mesh_diocotron_instability_structured_full',...
   '_k',num2str(Parameters.Degree),...
   '_r',num2str(Parameters.MeshRefinement)]};
Parameters.Mesh.File=Mesh.File;
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-6;                           % Tolerance
System.MaxIterations=5; % =10 in the paper       % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=4;                                % Final time
Time.TimeStepSize=1/2^3; % =1/2^8 in the paper   % Time step size
Time.BDFOrder=2;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Electrostatic initialization and projection correction -------------------------------------------
Parameters.ElectrostaticInitialization=...       % Electrostatic initialization
  @(Parameters) electrostaticInitialization(Parameters);
Parameters.ProjectionCorrection=...              % Projection correction
  @(Parameters,Block) projectionCorrection(Parameters,Block);
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
Solver.Equilibrate='yes';                        % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_r=[];                       % Dirichlet portion for density
Boundaries.Dirichlet_w=[];                       % Dirichlet portion for momentum
Boundaries.Dirichlet_g=[];                       % Dirichlet portion for energy
Boundaries.Dirichlet_e=1;                        % Dirichlet portion for electric field
Boundaries.FarField=[];                          % Far-field portion
Boundaries.InviscidWall=1;                       % Inviscid wall portion
Boundaries.Absorbing=[];                         % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='yes';                      % Plot geometry
Options.PlotMesh='yes';                          % Plot mesh
Options.PlotMeshDistortion='yes';                % Plot mesh distortion
Options.ComputeQuantityIterationEnd=...          % Compute quantity of interest at the iteration end
  ['[ResultsProjectionCorrection]=Parameters.ProjectionCorrection(Parameters,Block);',...
   'Block.SolutionLocal(:,6:7)=ResultsProjectionCorrection.ElectricFieldCorrected;'];
Options.Export2Paraview=...                      % Export to Paraview
  {'Density';
   'Momentum';
   'Energy';
   'MagneticField';
   'ElectricField'};
Options.Export2ParaviewTimeSteps='Time.TimeStep';% Export to Paraview time steps
Options.StoreTimeSteps='Time.TimeStep';          % Store time steps
Options.SaveResults='yes';                       % Save results
Options.SaveResultsFolder=...                    % Save results folder
  'diocotron_instability/';
Options.SaveResultsMatVersion='-v7.3';           % Save results (MAT-file version)
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Electrostatic initialization

function [Results] = electrostaticInitialization(ParametersAux)

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='ElectrostaticInitialization';% Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_HDG';            % Formulation
Parameters.Problem='ElectrostaticInitialization';% Problem
Parameters.Axisymmetric=...                      % Axisymmetric
  ParametersAux.Axisymmetric;
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.Degree=...                            % Degree
  ParametersAux.Degree;
Parameters.StabTemperature=1;                    % Stabilization for temperature
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [0*x,0*x];
Parameters.Temperature=@(x,y,z,t) 0*x;           % Temperature
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) 1/ParametersAux.ElectricPermittivity*...
         (ParametersAux.ParticleCharge/ParametersAux.ParticleMass*ParametersAux.Density(x,y,z,t)+...
          ParametersAux.IonCharge/ParametersAux.IonMass*ParametersAux.IonDensity);
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File=ParametersAux.Mesh.File;               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
Solver.Equilibrate='yes';                        % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=1;                          % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=[];                             % Robin portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  'Results.ElectricFieldInitial=Results.ScaledTemperatureGradient;';
Options.SaveResults='no';                        % Save results
% --------------------------------------------------------------------------------------------------

%% Main

fprintf('\n')
main

end

%% Projection correction

function [Results] = projectionCorrection(ParametersAux,BlockAux)

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='ProjectionCorrection';       % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_HDG';            % Formulation
Parameters.Problem='ProjectionCorrection';       % Problem
Parameters.Axisymmetric=...                      % Axisymmetric
  ParametersAux.Axisymmetric;
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.Degree=...                            % Degree
  ParametersAux.Degree;
Parameters.StabTemperature=1;                    % Stabilization for temperature
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [0*x,0*x];
Parameters.Temperature=@(x,y,z,t) 0*x;           % Temperature
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=@(x,y,z,t) 0*x;            % Heat source
Parameters.ElectricPermittivity=...              % Electric permittivity
  ParametersAux.ElectricPermittivity;
Parameters.ParticleCharge=...                    % Particle charge
  ParametersAux.ParticleCharge;
Parameters.IonCharge=...                         % Ion charge
  ParametersAux.IonCharge;
Parameters.ParticleMass=...                      % Particle mass
  ParametersAux.ParticleMass;
Parameters.IonMass=...                           % Ion mass
  ParametersAux.IonMass;
Parameters.IonDensity=...                        % Ion density
  ParametersAux.IonDensity;
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File=ParametersAux.Mesh.File;               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
Solver.Equilibrate='yes';                        % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=1;                          % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=[];                             % Robin portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeQuantityStart=...                 % Compute quantity of interest (start of simulation)
  ['Density=BlockAux.SolutionLocal(:,1);',...
   'ElectricField=BlockAux.SolutionLocal(:,6:7);',...
   'Elements(2).SolutionLocal=mat2cell([Density(Mesh.Elements(:),:),',...
   'ElectricField(Mesh.Elements(:),:)],ones(Sizes.NumElements,1)*Sizes.NumElementNodes,3);'];
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[~,~,DG2CG]=uniquetol(Mesh.Nodes'',1e-6,"ByRows",true);',...
   'ElementsCG=changem(Mesh.Elements,DG2CG,1:Sizes.NumLocalNodes);',...
   'CorrectionDG=Results.ScaledTemperatureGradient;',...
   'CorrectionCG(:,1)=accumarray(DG2CG,CorrectionDG(:,1),[],@mean);',...
   'CorrectionCG(:,2)=accumarray(DG2CG,CorrectionDG(:,2),[],@mean);',...
   'CorrectionCG2DG=CorrectionCG(ElementsCG(1:end),:);',...
   'Results.Correction=CorrectionCG2DG;',...
   'ElectricFieldCorrectedDG=ElectricField-Results.Correction;',...
   'ElectricFieldCorrectedCG(:,1)=accumarray(DG2CG,ElectricFieldCorrectedDG(:,1),[],@mean);',...
   'ElectricFieldCorrectedCG(:,2)=accumarray(DG2CG,ElectricFieldCorrectedDG(:,2),[],@mean);',...
   'ElectricFieldCorrectedCG2DG=ElectricFieldCorrectedCG(ElementsCG(1:end),:);',...
   'Results.ElectricFieldCorrected=ElectricFieldCorrectedCG2DG;'];
Options.SaveResults='no';                        % Save results
% --------------------------------------------------------------------------------------------------

%% Main

fprintf('\n')
main

end