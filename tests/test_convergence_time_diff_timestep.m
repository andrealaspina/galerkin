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
Parameters.Temperature=@(x,y,z,t) (x+y)*t^4;     % Temperature
Parameters.ThermalFlux=@(x,y,z,t) 0*x;           % Thermal flux
Parameters.HeatSource=@(x,y,z,t) (x+y)*4*t^3;    % Heat source
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_1'};                    % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1./2.^[1,2,3,4;...             % Time step size
                         1,2,NaN,NaN;
                         1,2,3,NaN];
Time.BDFOrder=[1,2,3];                             % BDF order
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
  'abs(Results.TemperatureErrorL2Rate(3,3)-Time.BDFOrder(end))<1e-3';
% --------------------------------------------------------------------------------------------------

%% Main

main