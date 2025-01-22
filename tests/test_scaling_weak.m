clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ScalingWeak';                   % Simulation type
Simulation.Problem='Fluid';                      % Problem
PoolAux=gcp; if isempty(PoolAux); PoolAux=[]; PoolAux.NumWorkers=1; end
Simulation.NumProcessors=PoolAux.NumWorkers*[1,1];% Number of processors
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='WeaklyCompressibleFlowDM_HDG';% Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=5;                             % Degree
Parameters.StabDensity=1000;                     % Stabilization for density
Parameters.StabMomentum=1;                       % Stabilization for momentum
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=1/300;           % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.ChannelLength=10;                     % Channel length
Parameters.ChannelHalfHeight=1;                  % Channel half height
Parameters.ChannelExitMeanVelocity=1;            % Channel exit mean velocity
r0=Parameters.ReferenceDensity;      p0=Parameters.ReferencePressure;
eps=Parameters.CompressibilityCoeff; mu=Parameters.DynamicViscosity;
L=Parameters.ChannelLength;          R=Parameters.ChannelHalfHeight;
U=Parameters.ChannelExitMeanVelocity;
vx=@(x,y) 3/2*U*(1-(y/R).^2)...
          -9/2*(mu*L*U^2)/(r0*R^2)*(1-x/L).*(1-(y/R).^2)*eps;
vy=@(x,y) 0*x;
p= @(x,y) p0+3*(mu*L*U)/(R^2)*(1-x/L)...
          -3/2*(mu^2*U^2)/(r0*R^2)*(3*(L/R)^2*(1-x/L).^2-(1-(y/R).^2))*eps;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [-(3*U^2*eps*mu^(3/2)*(R^2-y.^2)*(3*2^(1/2)+6^(1/2)))/(4*R^4*r0),...
              +(3*U^2*eps*mu^(3/2)*(R^2-y.^2)*(3*2^(1/2)-6^(1/2)))/(4*R^4*r0),...
              +(3*U*mu^(1/2)*y.*(R^2*r0-3*L*U*eps*mu+3*U*eps*mu*x))/(R^4*r0)];
Parameters.Density=@(x,y,z,t) r0+eps*(p(x,y)-p0);% Density
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [Parameters.Density(x,y,t).*vx(x,y),...
              Parameters.Density(x,y,t).*vy(x,y)];
Parameters.Velocity=@(x,y,z,t) [vx(x,y),vy(x,y)];% Velocity
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [-((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*...
                (+p0-(3*L*U*mu*(x/L-1))/R^2-(3*U^2*eps*mu^2*(y.^2/R^2+(3*L^2*(x/L-1).^2)/R^2-1))/...
                (2*R^2)+(6*U^2*eps*mu^2*(y.^2/R^2-1))/R^2)-...
               ((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y)    .*((3*U*y)/R^2+...
                (9*L*U^2*eps*mu*y.*(x/L-1))/R^4)*mu,...
              -((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*((3*U*y)/R^2+...
                (9*L*U^2*eps*mu*y.*(x/L-1))/R^4)*mu+...
               ((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y)    .*...
                (-p0+(3*L*U*mu*(x/L-1))/R^2+(3*U^2*eps*mu^2*(y.^2/R^2+(3*L^2*(x/L-1).^2)/R^2-1))/...
                (2*R^2)+(3*U^2*eps*mu^2*(y.^2/R^2-1))/R^2)];
Parameters.ResidualContinuity=...                % Residual of continuity equation
  @(x,y,z,t) (27*mu^2*U^3)/(4*r0*R^8)*(R^2-y.^2).*...
             (6*R^2*(L-x)-mu*U/r0*(9*(L-x).^2-R^2+y.^2)*eps)*eps^2;
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
Parameters.Pressure=@(Density)...                % Equation of state
  Parameters.ReferencePressure+...
 (Density-Parameters.ReferenceDensity)/Parameters.CompressibilityCoeff;
Parameters.DPressureDDensity=@(Density)...       % dPressure/dDensity
  1/Parameters.CompressibilityCoeff;
clear r0 p0 eps mu L R U vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_scaling_weak_1',...              % Mesh file
          'Mesh_scaling_weak_1'};
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-6;                           % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
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
Options.PlotEfficiency='yes';                    % Plot efficiency
Options.CharacteristicTimer='Evaluation';        % Characteristic timer
Options.Test='Results.Efficiency(1)==1';         % Test
% --------------------------------------------------------------------------------------------------

%% Main

main