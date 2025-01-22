clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetic';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Magnetic_HDG';           % Formulation
Parameters.Problem='Magnetic';                   % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveTerm='yes';                 % Convective term
Parameters.Degree=1;                             % Degree
Parameters.StabMagneticInduction=10;             % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=10;            % Stabilization for Lagrange multiplier
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.Hartmann=1;                           % Hartmann number
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.ChannelLength=10;                     % Channel length
Parameters.ChannelHalfHeight=1;                  % Channel half height
Parameters.ChannelExitMeanVelocity=1;            % Channel exit mean velocity
eta=Parameters.MagneticDiffusivity;   mu_m=Parameters.MagneticPermeability;
mu=Parameters.DynamicViscosity;
L=Parameters.ChannelLength;           R=Parameters.ChannelHalfHeight;
U=Parameters.ChannelExitMeanVelocity; Ha=Parameters.Hartmann;
Parameters.ScaledMagneticGradient=...            % Scaled magnetic gradient
  @(x,y,z,t) [0*x,...
              0*x,...
              3*sqrt(mu*mu_m)*U/R*(sinh(Ha)-Ha*cosh(y/R*Ha))*(Ha*sinh(Ha))^(-1)];
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [3*sqrt(mu*mu_m/eta)*U*(sinh(y/R*Ha)/sinh(Ha)-y/R)*Ha^(-1),...
                sqrt(mu*mu_m*eta)/R*Ha*(x==x)];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.PseudoTraction=...                    % Pseudo-traction
  @(x,y,z,t) [-((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y)    .*3*sqrt(mu*mu_m*eta)/R*U.*...
                (sinh(Ha)-Ha*cosh(y/R*Ha))*(Ha*sinh(Ha))^(-1),...
              -((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*3*sqrt(mu*mu_m*eta)/R*U.*...
                (sinh(Ha)-Ha*cosh(y/R*Ha))*(Ha*sinh(Ha))^(-1)];
Parameters.Source=@(x,y,z,t) [0*x,0*x];          % Source
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [3/2*U*2*(1-cosh(y/R*Ha)/cosh(Ha))*(Ha*tanh(Ha))^(-1),...
              0*x];
clear eta mu_m Ha mu L R U
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_hartmann_3'};                    % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
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
  {'ScaledMagneticGradient';
   'MagneticInduction';
   'LagrangeMultiplier';
   'MagneticInductionPost'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledMagneticGradientErrorL2-9.427830133993462e-02)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorL2     -7.997581579478961e-03)<1e-12 && ',...
   'abs(Results.LagrangeMultiplierErrorL2    -1.723236638790498e-03)<1e-12 && ',...
   'abs(Results.MagneticInductionPostErrorL2 -1.100731614961559e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main