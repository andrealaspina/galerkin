clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_HDG';            % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.Degree=1;                             % Degree
Parameters.StabTemperature=1;                    % Stabilization for temperature
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [-1*(x==x),-1*(x==x)];
Parameters.Temperature=@(x,y,z,t) x+y;           % Temperature
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=@(x,y,z,t) 0*x;            % Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='pcg';                               % Type
Solver.Preconditioner='ichol';                   % Preconditioner
Solver.PreconditionerAllTimeSteps='no';          % Compute preconditioner at all time steps
Solver.PreconditionerAllIterations='no';         % Compute preconditioner at all iterations
Solver.Tolerance=1e-6;                           % Tolerance
Solver.MaxIterations=100;                        % Maximum number of iterations
Solver.Equilibrate='no';                         % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledTemperatureGradient';
   'Temperature'};
Options.Test=...                                 % Test
  ['Results.ScaledTemperatureGradientErrorL2<4e-8 && ',...
   'Results.TemperatureErrorL2              <5e-9'];
% --------------------------------------------------------------------------------------------------

%% Main

main