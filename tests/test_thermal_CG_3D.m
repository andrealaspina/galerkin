clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) sin(pi*x).*sin(pi*y).*sin(pi*z);
Parameters.TemperatureGradient=...               % Temperature gradient
  @(x,y,z,t) [pi*cos(pi*x).*sin(pi*y).*sin(pi*z),...
              pi*sin(pi*x).*cos(pi*y).*sin(pi*z),...
              pi*sin(pi*x).*sin(pi*y).*cos(pi*z)];
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) 3*pi^2*sin(pi*x).*sin(pi*y).*sin(pi*z);
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_cube_1'};                        % Mesh file
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
Boundaries.Dirichlet=[1,2,3,4,5,6];              % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError={'Temperature'};            % Compute error
Options.Test=...                                 % Test
  'abs(Results.TemperatureErrorL2-2.744756647666941e-01)<1e-12';
% --------------------------------------------------------------------------------------------------

%% Main

main