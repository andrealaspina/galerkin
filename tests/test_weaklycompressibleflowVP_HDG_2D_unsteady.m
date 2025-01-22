clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Fluid';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='WeaklyCompressibleFlowVP_HDG';% Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveFlow='yes';                 % Convective flow
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=1;                             % Degree
Parameters.StabVelocity=1;                       % Stabilization for velocity
Parameters.StabPressure=10;                      % Stabilization for pressure
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=0.1;             % Compressibility coefficient
Parameters.DynamicViscosity=0.1;                 % Dynamic viscosity
r0=Parameters.ReferenceDensity;      p0=Parameters.ReferencePressure;
eps=Parameters.CompressibilityCoeff; mu=Parameters.DynamicViscosity;
vx=@(x,y,t) sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0));
vy=@(x,y,t) cos(pi.*x).*cos(pi.*y).*sin(pi.*t)-eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0));
p= @(x,y,t) pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y);
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [(mu.^(1./2).*pi.*(6.^(1./2).*eps.*pi.*cos(pi.*x).^2.*sin(pi.*t).^2-6.^(1./2).*eps.*pi.*sin(pi.*t).^2+6.^(1./2).*eps.*pi.*cos(pi.*y).^2.*sin(pi.*t).^2-6.*2.^(1./2).*r0.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)+6.^(1./2).*eps.*pi.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)))./(6.*r0),...
              (mu.^(1./2).*pi.*(6.^(1./2).*eps.*pi.*cos(pi.*x).^2.*sin(pi.*t).^2-6.^(1./2).*eps.*pi.*sin(pi.*t).^2+6.^(1./2).*eps.*pi.*cos(pi.*y).^2.*sin(pi.*t).^2+6.*2.^(1./2).*r0.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)+6.^(1./2).*eps.*pi.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)))./(6.*r0),...
              -(eps.*mu.^(1./2).*pi.^2.*(x.*pi.*sin(pi.*t).^2.*sin(2.*pi.*y)-2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x)+y.*pi.*sin(pi.*t).^2.*sin(2.*pi.*x)))./(2.*r0)];
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [vx(x,y,t),vy(x,y,t)];
Parameters.Pressure=@(x,y,z,t) p(x,y,t);         % Pressure
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [0*x, 0*x];
Parameters.ResidualContinuity=...                % Residual of continuity equation
  @(x,y,z,t) pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t));
Parameters.Force=...                             % Force
  @(x,y,z,t) [(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-mu.*(eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0)+(pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))+eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^4.*x.*cos(2.*pi.*y).*sin(pi.*t).^2)./r0)+(2.*eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^3.*cos(pi.*x).*sin(pi.*t).^2.*sin(pi.*x))./r0))./3-2.*pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-(eps.*((pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(4.*r0)-(pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*x.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0))-pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).^2+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)),...
              mu.*(eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^4.*y.*cos(2.*pi.*x).*sin(pi.*t).^2)./r0)-eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))+(2.*eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*y).*sin(pi.*t).^2.*sin(pi.*y))./r0))./3+2.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-(eps.*((pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))./(2.*r0)+(pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(4.*r0))-pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*y.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0))+pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t)+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)).^2-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))];
Parameters.Density=@(Pressure)...                % Equation of state
  Parameters.ReferenceDensity+...
 (Pressure-Parameters.ReferencePressure)*Parameters.CompressibilityCoeff;
Parameters.DDensityDPressure=@(Pressure)...      % dDensity/dPressure
  Parameters.CompressibilityCoeff;
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
Boundaries.Dirichlet_v_x=[1,2,3,4];              % Dirichlet portion for x-velocity
Boundaries.Dirichlet_v_y=[1,2,3,4];              % Dirichlet portion for y-velocity
Boundaries.Dirichlet_p=[1,2,3,4];                % Dirichlet portion for pressure
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2-9.613585201526805e-01)<1e-12 && ',...
   'abs(Results.VelocityErrorL2        -4.698814311828181e-01)<1e-12 && ',...
   'abs(Results.PressureErrorL2        -5.728574451786637e-01)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main