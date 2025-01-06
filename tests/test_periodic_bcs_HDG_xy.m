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
Parameters.Degree=7;                             % Degree
Parameters.StabTemperature=1;                    % Stabilization for temperature
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [-(y.*z.*(6.*x.^2-6.*x+1).*(2.*y.^2-3.*y+1))/4,...
              -(x.*z.*(2.*x.^2-3.*x+1).*(6.*y.^2-6.*y+1))/4,...
              -x.*y.*(x-1).*(x-1/2).*(y-1).*(y-1/2)];
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) x.*(x-1/2).*(x-1).*y.*(y-1/2).*(y-1).*z;
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) (3.*z.*(2.*x-1).*(2.*y-1).*(-x.^2+x-y.^2+y))/2;
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_cube_1'};                        % Mesh file
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
Boundaries.Dirichlet=[2,4];                      % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.PeriodicMaster=[1,3];                 % Periodic portion (master)
Boundaries.PeriodicSlave=[5,6];                  % Periodic portion (slave)
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
  ['Results.ScaledTemperatureGradientErrorL2<1e-12 && ',...
   'Results.TemperatureErrorL2              <1e-12 && ',...
   'Results.TemperaturePostErrorL2          <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main