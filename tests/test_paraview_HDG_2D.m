clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_HDG';            % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.NodesDistribution='Uniform';          % Nodes distribution
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=4;                             % Degree
Parameters.StabTemperature=1;                    % Stabilization for temperature
Parameters.Density=1;                            % Density
Parameters.SpecificHeatCapacity=1;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [-4*x.^3,-4*y.^3];
Parameters.Temperature=@(x,y,z,t) x.^4+y.^4;     % Temperature
Parameters.ThermalFlux=@(x,y,z,t) 0*x;           % Thermal flux
Parameters.HeatSource=@(x,y,z,t)-12*x.^2-12*y.^2;% Heat source
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
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview=...                      % Export to Paraview
  {'ThermalConductivity';
   'HeatSource';
   'ScaledTemperatureGradient';
   'Temperature';
   'TemperaturePost'};
Options.ComputeError=...                         % Compute error
  {'ScaledTemperatureGradient';
   'Temperature';
   'TemperaturePost'};
Options.SaveResultsFolder='tests/';              % Save results folder
Options.Test=...                                 % Test
  ['Results.ScaledTemperatureGradientErrorL2<1e-12 && ',...
   'Results.TemperatureErrorL2              <1e-12 && ',...
   'Results.TemperaturePostErrorL2          <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main