clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ConvergenceSpace';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='Thermal_HDG';         % Formulation
Parameters(1).Problem='Thermal';                 % Problem
Parameters(1).PostProcessingHDG='yes';           % Perform HDG postprocessing
Parameters(1).Degree=2;                          % Degree
Parameters(1).StabTemperature=10;                % Stabilization for temperature
Parameters(1).Density=1;                         % Density
Parameters(1).SpecificHeatCapacity=1;            % Specific heat capacity
Parameters(1).ThermalConductivity=1;             % Thermal conductivity
Parameters(1).ScaledTemperatureGradient=...      % Scaled temperature gradient
  @(x,y,z,t) [-2*x,-2*y];
Parameters(1).Temperature=@(x,y,z,t) x.^2+y.^2;  % Temperature
Parameters(1).ThermalFlux=@(x,y,z,t) -2*x;       % Thermal flux
Parameters(1).HeatSource=@(x,y,z,t) -4*(x==x);   % Heat source
Parameters(2).Formulation='Thermal_CG';          % Formulation
Parameters(2).Problem='Thermal';                 % Problem
Parameters(2).Degree=2;                          % Degree
Parameters(2).NitschePenalty=100;                % Nitsche's penalty parameter
Parameters(2).Density=1;                         % Density
Parameters(2).SpecificHeatCapacity=1;            % Specific heat capacity
Parameters(2).ThermalConductivity=1;             % Thermal conductivity
Parameters(2).Temperature=@(x,y,z,t) x.^2+y.^2;  % Temperature
Parameters(2).ThermalFlux=@(x,y,z,t) +2*x;       % Thermal flux
Parameters(2).HeatSource=@(x,y,z,t) -4*(x==x);   % Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_coupled_2parts_1',...   % Mesh file
           'Mesh_square_coupled_2parts_2',...
           'Mesh_square_coupled_2parts_3'};
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
Boundaries(1).Dirichlet=[2,4;1,4;1,4];           % Dirichlet portion
Boundaries(1).Interface=[1;2;2];                 % Interface portion
Boundaries(1).Neumann=  [3;3;3];                 % Neumann portion
Boundaries(2).Dirichlet=[1,4;2,4;1,4];           % Dirichlet portion
Boundaries(2).Interface=[3;1;3];                 % Interface portion
Boundaries(2).Neumann=  [2;3;2];                 % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.PlotConvergence='no';                    % Plot convergence
Options.CharacteristicElementSize='Min';         % Characteristic element size
Options.ComputeError=...                         % Compute error
  {'ScaledTemperatureGradient';
   'Temperature';
   'TemperaturePost'};
Options.Test=...                                 % Test
  ['max(Results(1).ScaledTemperatureGradientErrorL2)<1e-11 && ',...
   'max(Results(1).TemperatureErrorL2              )<1e-12 && ',...
   'max(Results(1).TemperaturePostErrorL2          )<1e-12 && ',...
   'max(Results(2).TemperatureErrorL2              )<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main