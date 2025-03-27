clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.NodesDistribution='Fekete';           % Nodes distribution
Parameters.Degree=2;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.SpecificHeatCapacity=1;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.Temperature=@(x,y,z,t) x.^2+y.^2;     % Temperature
Parameters.ThermalFlux=@(x,y,z,t) 0*x;           % Thermal flux
Parameters.HeatSource=@(x,y,z,t) -4*(x==x);      % Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_1'};                     % Mesh file
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
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=[];                             % Robin portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError={'Temperature'};            % Compute error
Options.Test=...                                 % Test
  'Results.TemperatureErrorL2<1e-12';
% --------------------------------------------------------------------------------------------------

%% Main

main