clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ConvergenceSpace';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=[1,2,3];                       % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.SpecificHeatCapacity=1;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.Temperature=@(x,y,z,t) x.^4+y.^4;     % Temperature
Parameters.ThermalFlux=@(x,y,z,t) 0*x;           % Thermal flux
Parameters.HeatSource=@(x,y,z,t)-12*x.^2-12*y.^2;% Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=decsg([3,4,...                          % Geometry definition
               -1,+1,+1,-1,...
               +1,+1,-1,-1]','R','R');
Mesh.MinElementSize=...                          % Minimum element size
  2./2.^[1,2,3,4;...
         1,2,NaN,NaN;...
         1,2,3,NaN];
Mesh.MaxElementSize=Mesh.MinElementSize;         % Maximum element size
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
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=[];                             % Robin portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.PlotConvergence='no';                    % Plot convergence
Options.CharacteristicElementSize='Min';         % Characteristic element size
Options.ComputeError={'Temperature'};            % Compute error
Options.Test=...                                 % Test
  'abs(Results.TemperatureErrorL2Rate(3,3)-(Parameters.Degree(end)+1))<2.3e-1';
% --------------------------------------------------------------------------------------------------

%% Main

main