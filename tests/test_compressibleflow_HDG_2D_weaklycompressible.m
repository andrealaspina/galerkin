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
Parameters.SolveEnergyEquation='no';             % Solve energy equation
Parameters.Degree=1;                             % Degree
Parameters.StabDensity=100;                      % Stabilization for density
Parameters.StabMomentum=1;                       % Stabilization for momentum
Parameters.StabEnergy=1;                         % Stabilization for energy
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=0.1;             % Compressibility coefficient
Parameters.DynamicViscosity=0.1;                 % Dynamic viscosity
Parameters.ThermalConductivity=1;                % Thermal conductivity
r0=Parameters.ReferenceDensity;      p0=Parameters.ReferencePressure;
eps=Parameters.CompressibilityCoeff; mu=Parameters.DynamicViscosity;
vx=@(x,y,t) sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0));
vy=@(x,y,t) cos(pi.*x).*cos(pi.*y).*sin(pi.*t)-eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0));
p= @(x,y,t) pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y);
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [(mu.^(1./2).*pi.*(6.^(1./2).*eps.*pi.*cos(pi.*x).^2.*sin(pi.*t).^2-6.^(1./2).*eps.*pi.*sin(pi.*t).^2+6.^(1./2).*eps.*pi.*cos(pi.*y).^2.*sin(pi.*t).^2-6.*2.^(1./2).*r0.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)+6.^(1./2).*eps.*pi.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)))./(6.*r0),...
              (mu.^(1./2).*pi.*(6.^(1./2).*eps.*pi.*cos(pi.*x).^2.*sin(pi.*t).^2-6.^(1./2).*eps.*pi.*sin(pi.*t).^2+6.^(1./2).*eps.*pi.*cos(pi.*y).^2.*sin(pi.*t).^2+6.*2.^(1./2).*r0.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)+6.^(1./2).*eps.*pi.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)))./(6.*r0),...
              -(eps.*mu.^(1./2).*pi.^2.*(x.*pi.*sin(pi.*t).^2.*sin(2.*pi.*y)-2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x)+y.*pi.*sin(pi.*t).^2.*sin(2.*pi.*x)))./(2.*r0)];
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [0*x,0*x];
Parameters.Density=...                           % Density
  @(x,y,z,t) r0+eps*(p(x,y,t)-p0);
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [Parameters.Density(x,y,z,t).*vx(x,y,t),...
              Parameters.Density(x,y,z,t).*vy(x,y,t)];
Parameters.Energy=@(x,y,z,t) 0*x;                % Energy
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [vx(x,y,t),vy(x,y,t)];
Parameters.Temperature=@(x,y,z,t) 0*x;           % Temperature
Parameters.ResidualContinuity=...                % Residual of continuity equation
  @(x,y,z,t) pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t));
Parameters.Force=...                             % Force
  @(x,y,z,t) [(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-mu.*(eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0)+(pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))+eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^4.*x.*cos(2.*pi.*y).*sin(pi.*t).^2)./r0)+(2.*eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^3.*cos(pi.*x).*sin(pi.*t).^2.*sin(pi.*x))./r0))./3-2.*pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-(eps.*((pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(4.*r0)-(pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*x.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0))-pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).^2+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)),...
              mu.*(eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^4.*y.*cos(2.*pi.*x).*sin(pi.*t).^2)./r0)-eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))+(2.*eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*y).*sin(pi.*t).^2.*sin(pi.*y))./r0))./3+2.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-(eps.*((pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))./(2.*r0)+(pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(4.*r0))-pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*y.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0))+pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t)+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)).^2-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))];
Parameters.HeatSource=@(x,y,z,t) 0*x;            % Heat source
Parameters.ComputePressure=...                   % Compute pressure
  @(Density,Momentum,Energy)...
   Parameters.ReferencePressure+...
   (Density-Parameters.ReferenceDensity)/Parameters.CompressibilityCoeff;
Parameters.ComputePressureLinearization=...      % Compute pressure linearization
  @(Density,Momentum,Energy)...
   [1/Parameters.CompressibilityCoeff*ones(size(Density)),...
    zeros(size(Momentum)),...
    zeros(size(Energy))];
Parameters.ComputeTemperature=...                % Compute temperature
  @(Density,Momentum,Energy)...
   zeros(size(Density));
Parameters.ComputeTemperatureLinearization=...   % Compute temperature linearization
  @(Density,Momentum,Energy)...
   [zeros(size(Density)),...
    zeros(size(Momentum)),...
    zeros(size(Energy))];
clear r0 p0 eps mu vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square01_1'};                    % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-3;                           % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=0.5;                              % Final time
Time.TimeStepSize=2^(-3);                        % Time step size
Time.BDFOrder=4;                                 % BDF order
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
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'ScaledTemperatureGradient';
   'Density';
   'Momentum';
   'Energy';
   'VelocityPost';
   'TemperaturePost'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2         -9.655484433551047e-01)<1e-10 && ',...
   'abs(Results.ScaledTemperatureGradientErrorL2-0.000000000000000e-00)<1e-10 && ',...
   'abs(Results.DensityErrorL2                  -5.689157679981924e-02)<1e-10 && ',...
   'abs(Results.MomentumErrorL2                 -4.574051777313588e-01)<1e-10 && ',...
   'abs(Results.EnergyErrorL2                   -0.000000000000000e-00)<1e-10 && ',...
   'abs(Results.VelocityPostErrorL2             -4.140188763306042e-01)<1e-10 && ',...
   'abs(Results.TemperaturePostErrorL2          -0.000000000000000e-00)<1e-10'];
% --------------------------------------------------------------------------------------------------

%% Main

main