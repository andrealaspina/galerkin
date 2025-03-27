clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_HDG';            % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=1;                             % Degree
Parameters.StabTemperature=10;                   % Stabilization for temperature
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ConvectionCoefficient=@(x,y,z,b) 0;   % Convection coefficient
Parameters.AmbientTemperature=@(x,y,z,b) 0;      % Ambient temperature
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [pi/2*x./sqrt(x.^2+y.^2).*sin(pi/2*sqrt(y.^2+x.^2)),...
              pi/2*y./sqrt(x.^2+y.^2).*sin(pi/2*sqrt(y.^2+x.^2))];
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) cos(pi/2*sqrt(x.^2+y.^2));
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) pi/2*(1./sqrt(x.^2+y.^2).*sin(pi/2*sqrt(x.^2+y.^2))+...
                                  pi/2*cos(pi/2*sqrt(x.^2+y.^2)));
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
Options.ComputeError=...                         % Compute error
  {'ScaledTemperatureGradient';
   'Temperature';
   'TemperaturePost'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledTemperatureGradientErrorL2-5.661404636327417e-01)<1e-12 && ',...
   'abs(Results.TemperatureErrorL2              -5.461589675232700e-02)<1e-12 && ',...
   'abs(Results.TemperaturePostErrorL2          -9.677073290852840e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main