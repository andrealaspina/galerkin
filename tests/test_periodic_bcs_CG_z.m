clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=5;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) z.*(z-1/2).*(z-1).*x.*y;
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) -3*(2*z-1).*x.*y;
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_cube_2'};                       % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,3,4,5];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=[];                             % Robin portion
Boundaries.PeriodicMaster=2;                     % Periodic portion (master)
Boundaries.PeriodicSlave=6;                      % Periodic portion (slave)
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