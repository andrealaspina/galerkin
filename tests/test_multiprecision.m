clearvars -except Files Passed CPUTime Test; addpath ../advanpix;

fprintf('\n'); mptest

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
Simulation.Digits=32;                            % Number of digits
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='Thermal_HDG';         % Formulation
Parameters(1).Problem='Thermal';                 % Problem
Parameters(1).PostProcessingHDG='yes';           % Perform HDG postprocessing
Parameters(1).NodesDistribution='Uniform';       % Nodes distribution
Parameters(1).Degree=3;                          % Degree
Parameters(1).StabTemperature=1;                 % Stabilization for temperature
Parameters(1).Density=1;                         % Density
Parameters(1).SpecificHeatCapacity=1;            % Specific heat capacity
Parameters(1).ThermalConductivity=1;             % Thermal conductivity
Parameters(1).ScaledTemperatureGradient=...      % Scaled temperature gradient
  @(x,y,z,t) [-3*x.^2*t,0*x];
Parameters(1).Temperature=@(x,y,z,t) x.^3*t;     % Temperature
Parameters(1).ThermalFlux=@(x,y,z,t) -3*x.^2*t;  % Thermal flux
Parameters(1).HeatSource=@(x,y,z,t) x.^3-6*x*t;  % Heat source
Parameters(2).Formulation='Thermal_CG';          % Formulation
Parameters(2).Problem='Thermal';                 % Problem
Parameters(2).NodesDistribution='Uniform';       % Nodes distribution
Parameters(2).Degree=4;                          % Degree
Parameters(2).NitschePenalty=10;                 % Nitsche's penalty parameter
Parameters(2).Density=1;                         % Density
Parameters(2).SpecificHeatCapacity=1;            % Specific heat capacity
Parameters(2).ThermalConductivity=1;             % Thermal conductivity
Parameters(2).Temperature=@(x,y,z,t) x.^3*t;     % Temperature
Parameters(2).ThermalFlux=@(x,y,z,t) +3*x.^2*t;  % Thermal flux
Parameters(2).HeatSource=@(x,y,z,t) x.^3-6*x*t;  % Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_multiprecision_coupled_1'};      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=mp('1/3');                     % Time step size
Time.BDFOrder=1;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries(1).Dirichlet=[2,4];                   % Dirichlet portion
Boundaries(1).Interface=3;                       % Interface portion
Boundaries(1).Neumann=1;                         % Neumann portion
Boundaries(2).Dirichlet=[2,4];                   % Dirichlet portion
Boundaries(2).Interface=1;                       % Interface portion
Boundaries(2).Neumann=3;                         % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError=...                         % Compute error
  {'ScaledTemperatureGradient';
   'Temperature';
   'TemperaturePost'};
Options.Test=...                                 % Test
  ['Results(1).ScaledTemperatureGradientErrorL2<1e-29 && ',...
   'Results(1).TemperatureErrorL2              <1e-30 && ',...
   'Results(1).TemperaturePostErrorL2          <1e-30 && ',...
   'Results(2).TemperatureErrorL2              <1e-30'];
% --------------------------------------------------------------------------------------------------

%% Main

main