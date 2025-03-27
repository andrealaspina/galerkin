close all; clear; clc;

%% Input (thermal problem)

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=1e2;                   % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=11.5;             % Thermal conductivity
Parameters.ConvectionCoefficient=...             % Convection coefficient
  @(x,y,z,b) 15*((b== 2)||(b== 6)||(b== 7)||(b== 8)||(b== 9))+...
             20*((b==13))+...
             30*((b==12)||(b==14)||(b==15))+...
             40*((b== 1)||(b==10))+...
             50*((b==11))+...
           1000*((b== 3)||(b== 4)||(b== 5));
Parameters.AmbientTemperature=...                % Ambient temperature
  @(x,y,z,b)150*((b==12)||(b==14)||(b==15))+...
            300*((b== 3)||(b== 4)||(b== 5))+...
            400*((b== 2)||(b== 6)||(b== 7)||(b== 8)||(b== 9))+...
            800*((b== 1))+...
           1000*((b==10)||(b==11)||(b==13));
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) 0*x;
Parameters.TemperatureGradient=...               % Temperature gradient
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=@(x,y,z,t) 0*x;            % Heat source               
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=importGeometry(createpde,'Blade.stl');  % Geometry definition
Mesh.MinElementSize=0.01/2;                      % Minimum element size
Mesh.MaxElementSize=0.01/2;                      % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
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
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[];                         % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=1:15;                           % Robin portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='yes';                      % Plot geometry
Options.PlotMesh='yes';                          % Plot mesh
Options.PlotMeshDistortion='yes';                % Plot mesh distortion
Options.PlotSolution={'Temperature'};            % Plot solution
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Temperature interpolant

TemperatureInterpolant=scatteredInterpolant(Mesh.Nodes(1,:)',Mesh.Nodes(2,:)',Mesh.Nodes(3,:)',...
  Results.Temperature,'linear','linear');

clearvars -except Geometry TemperatureInterpolant;

%% Input (structural problem)

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='LinearElasticity';             % Model
Parameters.Degree=2;                             % Degree
Parameters.NitschePenalty=1e12;                  % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.YoungsModulus=227e9;                  % Young's modulus
Parameters.PoissonsRatio=0.27;                   % Poisson's ratio
Parameters.CoefficientThermalExpansion=12.7e-6;  % Coefficient of thermal expansion
Parameters.ReferenceTemperature=300;             % Reference temperature
Parameters.Displacement=@(x,y,z,t) [0*x,0*x,0*x];% Displacement
Parameters.Traction=...                          % Traction
  @(x,y,z,t,b,nx,ny,nz) [-5.0e5*nx.*(b==11)-4.5e5*nx.*(b==10),...
                         -5.0e5*ny.*(b==11)-4.5e5*ny.*(b==10),...
                         -5.0e5*nz.*(b==11)-4.5e5*nz.*(b==10)];
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
Parameters.Temperature=...                       % Temperature
  @(x,y,z) TemperatureInterpolant(x,y,z);
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.MinElementSize=0.01/1; %0.01/2 for figure   % Minimum element size
Mesh.MaxElementSize=0.01/1; %0.01/2 for figure   % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
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
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=3;                          % Dirichlet portion
Boundaries.Neumann=[10,11];                      % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='yes';                          % Plot mesh
Options.PlotMeshDistortion='yes';                % Plot mesh distortion
Options.PlotSolution=...                         % Plot solution
  {'Displacement';
   'VonMisesStress'};
Options.Export2Paraview=...                      % Export to Paraview
  {'Temperature';
   'Displacement';
   'Stress';
   'VonMisesStress'};
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[Results]=Formulation{1}.evaluateStress(1,Results,Elements,Parameters,Mesh,RefElement,',...
   '                                        Sizes); ',...
   'Results.VonMisesStress=sqrt(',...
   '  1/2*((Results.Stress(:,1)-Results.Stress(:,5)).^2+',...
   '       (Results.Stress(:,5)-Results.Stress(:,9)).^2+',...
   '       (Results.Stress(:,9)-Results.Stress(:,1)).^2+',...
   '     6*(Results.Stress(:,2).^2+Results.Stress(:,3).^2+Results.Stress(:,6).^2)));'];
Options.SaveResults='yes';                       % Save results
Options.SaveResultsFolder='jet_engine_blade/';   % Save results folder
% --------------------------------------------------------------------------------------------------

%% Main

main