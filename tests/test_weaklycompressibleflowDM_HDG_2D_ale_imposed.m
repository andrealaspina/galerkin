clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Fluid';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='WeaklyCompressibleFlowDM_HDG';% Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveFlow='yes';                 % Convective flow
Parameters.ArbitraryLagrangianEulerian='yes';    % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=1;                             % Degree
Parameters.StabDensity=100;                      % Stabilization for density
Parameters.StabMomentum=1;                       % Stabilization for momentum
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
Parameters.Density=...                           % Density
  @(x,y,z,t) r0+eps*(p(x,y,t)-p0);
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [Parameters.Density(x,y,z,t).*vx(x,y,t),...
              Parameters.Density(x,y,z,t).*vy(x,y,t)];
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [0*x, 0*x];
Parameters.ResidualContinuity=...                % Residual of continuity equation
  @(x,y,z,t) pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t));
Parameters.Force=...                             % Force
  @(x,y,z,t) [(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-mu.*(eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0)+(pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))+eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^4.*x.*cos(2.*pi.*y).*sin(pi.*t).^2)./r0)+(2.*eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^3.*cos(pi.*x).*sin(pi.*t).^2.*sin(pi.*x))./r0))./3-2.*pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-(eps.*((pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(4.*r0)-(pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*x.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0))-pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).^2+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)),...
              mu.*(eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^4.*y.*cos(2.*pi.*x).*sin(pi.*t).^2)./r0)-eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))+(2.*eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*y).*sin(pi.*t).^2.*sin(pi.*y))./r0))./3+2.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-(eps.*((pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))./(2.*r0)+(pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(4.*r0))-pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*y.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0))+pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t)+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)).^2-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))];
Parameters.Pressure=@(Density)...                % Equation of state
  Parameters.ReferencePressure+...
 (Density-Parameters.ReferenceDensity)/Parameters.CompressibilityCoeff;
Parameters.DPressureDDensity=@(Density)...       % dPressure/dDensity
  1/Parameters.CompressibilityCoeff;
Parameters.Displacement=...                      % Displacement
  @(x,y,z,t) [1/4*sin(2*pi*x).*(1-cos(2*pi*y)).*(1-cos(2*pi*t))*1/8,...
              1/4*sin(2*pi*y).*(1-cos(2*pi*x)).*(1-cos(2*pi*t))*1/8];
Parameters.ScaledStrainRateCenter=...            % Scaled strain rate at the center
  @(t) Parameters.ScaledStrainRate(0.5,0.5,0,t);
Parameters.DensityCenter=...                     % Density at the center
  @(t) Parameters.Density(0.5,0.5,0,t);
Parameters.MomentumCenter=...                    % Momentum at the center
  @(t) Parameters.Momentum(0.5,0.5,0,t);
clear r0 p0 eps mu vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square01_2'};                    % Mesh file
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
Boundaries.Dirichlet_r=[1,2,3,4];                % Dirichlet portion for density
Boundaries.Dirichlet_w_x=[1,2,3,4];              % Dirichlet portion for x-momentum
Boundaries.Dirichlet_w_y=[1,2,3,4];              % Dirichlet portion for y-momentum
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRateCenter','Number';
   'DensityCenter'         ,'Number';
   'MomentumCenter'        ,'Number'};
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[~,Node]=ismember([0.5,0.5,0],Mesh.Nodes'',''rows''); ',...
   'Results.ScaledStrainRateCenter=Results.ScaledStrainRate(Node,:); ',...
   'Results.DensityCenter=Results.Density(Node,:); ',...
   'Results.MomentumCenter=Results.Momentum(Node,:);'];
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateCenterErrorNumber-9.275922687857809e-01)<1e-12 && ',...
   'abs(Results.DensityCenterErrorNumber         -7.788901240095458e-03)<1e-12 && ',...
   'abs(Results.MomentumCenterErrorNumber        -7.037614890904743e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main