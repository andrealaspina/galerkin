clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetohydrodynamic';        % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagnetohydrodynamicsCURL_HDG';% Formulation
Parameters.Problem='Magnetohydrodynamic';        % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.CouplingTerms='no';                   % Coupling terms
Parameters.Degree=3;                             % Degree
Parameters.StabVelocity=10;                      % Stabilization for velocity
Parameters.StabPressure=10;                      % Stabilization for pressure
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=1/300;           % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ChannelLength=10;                     % Channel length
Parameters.ChannelHalfHeight=1;                  % Channel half height
Parameters.ChannelExitMeanVelocity=1;            % Channel exit mean velocity
r0=Parameters.ReferenceDensity;      p0=Parameters.ReferencePressure;
eps=Parameters.CompressibilityCoeff; mu=Parameters.DynamicViscosity;
L=Parameters.ChannelLength;           R=Parameters.ChannelHalfHeight;
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
Parameters.Velocity=@(x,y,z,t) [vx(x,y),vy(x,y)];% Velocity
Parameters.Pressure=@(x,y,z,t) p(x,y);           % Pressure
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
Parameters.Density=@(Pressure)...                % Equation of state
  Parameters.ReferenceDensity+...
 (Pressure-Parameters.ReferencePressure)*Parameters.CompressibilityCoeff;
Parameters.DDensityDPressure=@(Pressure)...      % dDensity/dPressure
  Parameters.CompressibilityCoeff;
Parameters.ScaledMagneticCurl=@(x,y,z,t) 0*x;    % Scaled magnetic curl
Parameters.ScaledMagneticGradient=...            % Scaled magnetic gradient
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [0*x,0*x];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.PseudoTraction=@(x,y,z,t) [0*x,0*x];  % Pseudo-traction
Parameters.Source=@(x,y,z,t) [0*x,0*x];          % Source
clear r0 p0 eps mu eta mu_m L R U vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_poiseuille_1'};                 % Mesh file
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
Boundaries.Dirichlet_v_x=[1,3,4];                % Dirichlet portion for x-velocity
Boundaries.Dirichlet_v_y=[1,2,3,4];              % Dirichlet portion for y-velocity
Boundaries.Dirichlet_p=[];                       % Dirichlet portion for pressure
Boundaries.Neumann_t_x=2;                        % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
Boundaries.Dirichlet_b_x=[1,2,3,4];              % Dirichlet portion for x-magnetic induction
Boundaries.Dirichlet_b_y=[1,2,3,4];              % Dirichlet portion for y-magnetic induction
Boundaries.Dirichlet_q=[1,2,3,4];                % Dirichlet portion for lagrange multiplier
Boundaries.Neumann_s_x=[];                       % Neumann portion for x-pseudo-traction
Boundaries.Neumann_s_y=[];                       % Neumann portion for y-pseudo-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure';
   'VelocityPost';
   'ScaledMagneticCurl';
   'ScaledMagneticGradient';
   'MagneticInduction';
   'LagrangeMultiplier';
   'MagneticInductionPost'};
Options.Test=...                                 % Test
  ['Results.ScaledStrainRateErrorL2      <1e-11 && ',...
   'Results.VelocityErrorL2              <1e-12 && ',...
   'Results.PressureErrorL2              <1e-10 && ',...
   'Results.VelocityPostErrorL2          <1e-11 && ',...
   'Results.ScaledMagneticCurlErrorL2    <1e-12 && ',...
   'Results.ScaledMagneticGradientErrorL2<1e-12 && ',...
   'Results.MagneticInductionErrorL2     <1e-12 && ',...
   'Results.LagrangeMultiplierErrorL2    <1e-12 && ',...
   'Results.MagneticInductionPostErrorL2 <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main