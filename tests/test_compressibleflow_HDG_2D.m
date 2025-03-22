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
Parameters.Degree=7;                             % Degree
Parameters.StabDensity=10;                       % Stabilization for density
Parameters.StabMomentum=10;                      % Stabilization for momentum
Parameters.StabEnergy=10;                        % Stabilization for energy
Parameters.SpecificHeatConstantPressure=2;       % Specific heat at constant pressure
Parameters.SpecificHeatConstantVolume=1;         % Specific heat at constant volume
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.ThermalConductivity=1;                % Thermal conductivity
cp=Parameters.SpecificHeatConstantPressure; cv=Parameters.SpecificHeatConstantVolume;
mu=Parameters.DynamicViscosity;             kappa=Parameters.ThermalConductivity;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [0.*x,...
              0.*x,...
              -2.*mu.^(1./2).*t.*(x==x)];
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [(kappa.^(1./2).*t.^2.*(x-1))./cv,...
              (kappa.^(1./2).*t.^2.*(y-1))./cv];
Parameters.Density=...                           % Density
  @(x,y,z,t) t.*(x-1).*(x+1).*(y-1).*(y+1)+1;
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [t.*(t.*(x-1).*(x+1).*(y-1).*(y+1)+1).*(y-1),...
              t.*(t.*(x-1).*(x+1).*(y-1).*(y+1)+1).*(x-1)];
Parameters.Energy=...                            % Energy
  @(x,y,z,t) t.*(x-1).*(x+1).*(y-1).*(y+1)+1;
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [t.*(y-1),...
              t.*(x-1)];
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) -((t.^2.*x.^2)./2-t.^2.*x+(t.^2.*y.^2)./2-t.^2.*y+t.^2-1)./cv;
Parameters.ResidualContinuity=...                % Residual of continuity equation
  @(x,y,z,t) t.*(t.*(x-1).*(x+1).*(y-1)+t.*(x-1).*(x+1).*(y+1)).*(x-1)+t.*(t.*(x-1).*(y-1).*(y+1)+t.*(x+1).*(y-1).*(y+1)).*(y-1)+(x-1).*(x+1).*(y-1).*(y+1);
Parameters.Force=...                             % Force
  @(x,y,z,t) [(y-1).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)-((cp-cv).*((t.^2.*(2.*x-2).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1))./2-t.*(x-1).*(y-1).*(y+1)-t.*(x+1).*(y-1).*(y+1)+t.^3.*x.*(y.^2-1).*(x.^2-2.*x+y.^2-2.*y+2)))./cv+t.^2.*(x-1).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)+2.*t.^3.*x.*(y-1).^3.*(y+1)+t.*(x-1).*(x+1).*(y-1).^2.*(y+1)+2.*t.^3.*y.*(x.^2-1).*(x-1).*(y-1),...
              (x-1).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)-((cp-cv).*((t.^2.*(2.*y-2).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1))./2-t.*(x-1).*(x+1).*(y-1)-t.*(x-1).*(x+1).*(y+1)+t.^3.*y.*(x.^2-1).*(x.^2-2.*x+y.^2-2.*y+2)))./cv+t.^2.*(y-1).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)+2.*t.^3.*y.*(x-1).^3.*(x+1)+t.*(x-1).^2.*(x+1).*(y-1).*(y+1)+2.*t.^3.*x.*(y.^2-1).*(x-1).*(y-1)];
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) -(2.*cv.*t-cv+4.*cv.*t.^2+2.*cv.*t.^3+2.*cv.*t.^4+cv.*x.^2+cv.*y.^2-2.*kappa.*t.^2+4.*cv.*mu.*t.^2+cv.*t.*x.^2-6.*cv.*t.^2.*x-2.*cv.*t.^3.*x+2.*cv.*t.^4.*x+cv.*t.*y.^2-6.*cv.*t.^2.*y-2.*cv.*t.^3.*y+2.*cv.*t.^4.*y-2.*cv.*t.^2.*x.^2+4.*cv.*t.^2.*x.^3-2.*cv.*t.^2.*x.^4-6.*cv.*t.^4.*x.^2+4.*cv.*t.^4.*x.^3-2.*cv.*t.^2.*y.^2+4.*cv.*t.^2.*y.^3-2.*cv.*t.^2.*y.^4-6.*cv.*t.^4.*y.^2+4.*cv.*t.^4.*y.^3-cv.*x.^2.*y.^2-2.*cv.*t.*x-2.*cv.*t.*y-4.*cv.*t.^2.*x.^2.*y.^3-4.*cv.*t.^2.*x.^3.*y.^2+2.*cv.*t.^2.*x.^2.*y.^4+2.*cv.*t.^2.*x.^4.*y.^2+10.*cv.*t.^4.*x.^2.*y.^2-8.*cv.*t.^4.*x.^2.*y.^3-8.*cv.*t.^4.*x.^3.*y.^2+6.*cv.*t.^4.*x.^3.*y.^3+4.*cv.*t.^2.*x.*y+2.*cv.*t.^3.*x.*y-14.*cv.*t.^4.*x.*y+6.*cv.*t.^2.*x.*y.^2+6.*cv.*t.^2.*x.^2.*y-2.*cv.*t.^2.*x.*y.^3-2.*cv.*t.^2.*x.^3.*y+8.*cv.*t.^4.*x.*y.^2+8.*cv.*t.^4.*x.^2.*y+2.*cv.*t.^4.*x.*y.^3+2.*cv.*t.^4.*x.^3.*y-6.*cv.*t.^4.*x.*y.^4-6.*cv.*t.^4.*x.^4.*y+2.*cv.*t.^4.*x.*y.^5+2.*cv.*t.^4.*x.^5.*y)./cv;
Parameters.ComputePressure=...                   % Compute pressure
  @(Density,Momentum,Energy)...
   (Parameters.SpecificHeatConstantPressure/Parameters.SpecificHeatConstantVolume-1)*...
   (Energy-1/2*vecnorm(Momentum,2,2).^2./Density);
Parameters.ComputePressureLinearization=...      % Compute pressure linearization
  @(Density,Momentum,Energy)...
   [+(Parameters.SpecificHeatConstantPressure/Parameters.SpecificHeatConstantVolume-1)*...
    1/2*vecnorm(Momentum,2,2).^2./Density.^2,...
    -(Parameters.SpecificHeatConstantPressure/Parameters.SpecificHeatConstantVolume-1)*...
    Momentum./Density,...
    +(Parameters.SpecificHeatConstantPressure/Parameters.SpecificHeatConstantVolume-1)*...
    ones(size(Energy))];
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
clear cp cv mu kappa nu r0 vx0 vy0 e0
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_1'};                     % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-12;                          % Tolerance
System.MaxIterations=10;                         % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/2;                           % Time step size
Time.BDFOrder=3;                                 % BDF order
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
  ['Results.ScaledStrainRateErrorL2         <1e-09 && ',...
   'Results.ScaledTemperatureGradientErrorL2<1e-09 && ',...
   'Results.DensityErrorL2                  <1e-10 && ',...
   'Results.MomentumErrorL2                 <1e-10 && ',...
   'Results.EnergyErrorL2                   <1e-10 && ',...
   'Results.VelocityPostErrorL2             <1e-10 && ',...
   'Results.TemperaturePostErrorL2          <1e-10'];
% --------------------------------------------------------------------------------------------------

%% Main

main