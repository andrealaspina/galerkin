clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ConvergenceTime';               % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.SpecificHeatCapacity=1;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.Temperature=@(x,y,z,t) (x+y)*t^3;     % Temperature
Parameters.ThermalFlux=@(x,y,z,t) 0*x;           % Thermal flux
Parameters.HeatSource=@(x,y,z,t) (x+y)*3*t^2;    % Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File=...                                   % Mesh file
  {'Mesh_square_2','Mesh_square_2','Mesh_square_3';...
   'Mesh_square_2','Mesh_square_3','Mesh_square_4'};
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1./2.^[0,1,2];                 % Time step size
Time.BDFOrder=[1,2];                             % BDF order
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
Options.PlotConvergence='yes';                    % Plot convergence
Options.CharacteristicElementSize='Min';         % Characteristic element size
Options.ComputeError={'Temperature'};            % Compute error
Options.Test=...                                 % Test
  'abs(Results.TemperatureErrorL2Rate(end)-Time.BDFOrder(end))<6e-2';
% --------------------------------------------------------------------------------------------------

%% Main

main