clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=4;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.Amplitude=1;                          % Amplitude
Parameters.Frequency=2;                          % Frequency
A=Parameters.Amplitude; f=Parameters.Frequency;  omega=2*pi*f;
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) A*(x+1).*(x-1).*(y+1).*(y-1)*cos(omega*t);
Parameters.ThermalFlux=@(x,y,z,t) 0*x;           % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) -2*A*(x.^2+y.^2-2)*cos(omega*t);
clear f A omega
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=2;                                % Final time
Time.TimeStepSize=2^(-4);                        % Time step size
Time.BDFOrder=2;                                 % BDF order
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
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError={'Temperature'};            % Compute error
Options.StoreTimeSteps='1:1:Time.NumTimeSteps';  % Store time steps
Options.PlotSolutionTimeSteps='NaN';             % Plot solution at stored time steps
Options.ComputeQuantity=...                      % Compute quantity of interest
  ['[~,Node]=ismember([0,0,0],Mesh.Nodes'',''rows'');',...
   'Results.CenterTemperature(Time.TimeStep)=',...
   'Block.SolutionGlobal(Node,1);'];
Options.ComputeFFT=...                           % Compute FFT
  {'CenterTemperature';
   'Temperature'};
Options.Test=...                                 % Test
  ['Results.TemperatureErrorL2                                       <1e-12 && ',...
   'abs(abs(Results.CenterTemperatureFFT(5))-Parameters.Amplitude)   <1e-12 && '...
   'abs(max(abs(Results.TemperatureFFT(:,:,5)))-Parameters.Amplitude)<1e-12 && '...
   'abs(Results.Frequency(5)-Parameters.Frequency)                   <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main