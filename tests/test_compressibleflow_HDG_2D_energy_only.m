clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Fluid';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='CompressibleFlow_HDG';   % Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveFlow='yes';                 % Convective flow
Parameters.SolveEnergyEquation='yes';            % Solve energy equation
Parameters.SolveEnergyEquationOnly='yes';        % Solve energy equation only
Parameters.Degree=2;                             % Degree
Parameters.StabDensity=1;                        % Stabilization for density
Parameters.StabMomentum=1;                       % Stabilization for momentum
Parameters.StabEnergy=1;                         % Stabilization for energy
Parameters.SpecificHeatConstantPressure=1;       % Specific heat at constant pressure
Parameters.SpecificHeatConstantVolume=1;         % Specific heat at constant volume
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.ThermalConductivity=1;                % Thermal conductivity
cv=Parameters.SpecificHeatConstantVolume;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [-2*x/cv,-2*y/cv];
Parameters.Density=@(x,y,z,t) 1*(x==x);          % Density
Parameters.Momentum=@(x,y,z,t) [0*x,0*x];        % Momentum
Parameters.Energy=@(x,y,z,t) x.^2+y.^2;          % Energy
Parameters.Velocity=@(x,y,z,t) [0*x,0*x];        % Velocity
Parameters.Temperature=@(x,y,z,t) (x.^2+y.^2)/cv;% Temperature
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
Parameters.HeatSource=@(x,y,z,t) -4/cv*(x==x);   % Heat source
Parameters.ComputePressure=...                   % Compute pressure
  @(Density,Momentum,Energy)...
   zeros(size(Density));
Parameters.ComputePressureLinearization=...      % Compute pressure linearization
  @(Density,Momentum,Energy)...
   [zeros(size(Density)),...
    zeros(size(Momentum)),...
    zeros(size(Energy))];
Parameters.ComputeTemperature=...                % Compute temperature
  @(Density,Momentum,Energy)...
   1/Parameters.SpecificHeatConstantVolume*...
   (Energy./Density-1/2*vecnorm(Momentum,2,2).^2./Density.^2);
Parameters.ComputeTemperatureLinearization=...   % Compute temperature linearization
  @(Density,Momentum,Energy)...
   [+1/Parameters.SpecificHeatConstantVolume*...
    (-Energy./Density.^2+vecnorm(Momentum,2,2).^2./Density.^3),...
    -1/Parameters.SpecificHeatConstantVolume*Momentum./Density.^2,...
    +1/Parameters.SpecificHeatConstantVolume*1./Density];
clear cv
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-10;                          % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'ScaledTemperatureGradient';
   'Density';
   'Momentum';
   'Energy';
   'VelocityPost';
   'TemperaturePost'};
Options.Test=...                                 % Test
  ['Results.ScaledStrainRateErrorL2         <1e-12 && ',...
   'Results.ScaledTemperatureGradientErrorL2<1e-12 && ',...
   'Results.DensityErrorL2                  <1e-12 && ',...
   'Results.MomentumErrorL2                 <1e-12 && ',...
   'Results.EnergyErrorL2                   <1e-12 && ',...
   'Results.VelocityPostErrorL2             <1e-12 && ',...
   'Results.TemperaturePostErrorL2          <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main