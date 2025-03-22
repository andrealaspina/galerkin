clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Fluid';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='WeaklyCompressibleFlowVP_HDG';% Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=1;                             % Degree
Parameters.StabVelocity=10;                      % Stabilization for velocity
Parameters.StabPressure=10;                      % Stabilization for pressure
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=0;               % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
Parameters.Hartmann=1;                           % Hartmann number
Parameters.ChannelLength=10;                     % Channel length
Parameters.ChannelHalfHeight=1;                  % Channel half height
Parameters.ChannelExitMeanVelocity=1;            % Channel exit mean velocity
p0=Parameters.ReferencePressure;    mu=Parameters.DynamicViscosity;
eta=Parameters.MagneticDiffusivity; Ha=Parameters.Hartmann;
L=Parameters.ChannelLength;         R=Parameters.ChannelHalfHeight;
U=Parameters.ChannelExitMeanVelocity;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [0*x,...
              0*x,...
              3*sqrt(mu)*U/R*sinh(y/R*Ha)*sinh(Ha)^(-1)];
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [3/2*U*2*(1-cosh(y/R*Ha)/cosh(Ha))*(Ha*tanh(Ha))^(-1),...
              0*x];
Parameters.Pressure=...                          % Pressure
  @(x,y,z,t) p0+3*mu*L*U/R^2*(1-x/L)-9/2*mu*U^2/eta*(sinh(y/R*Ha)/sinh(Ha)-y/R).^2*Ha^(-2);
Parameters.Traction=@(x,y,z,t) [0*x, 0*x];       % Traction
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
Parameters.Force=...                             % Force
  @(x,y,z,t) [-3*mu*U/R^2*(1-Ha*cosh(y/R*Ha)/sinh(Ha)),...
              -9*mu*U^2/(eta*R^2)*(y*sinh(Ha)-R*sinh(y/R*Ha)).*...
                                   (sinh(Ha)-Ha*cosh(y/R*Ha))*(Ha*sinh(Ha))^(-2)];
Parameters.Density=@(Pressure)...                % Equation of state
  Parameters.ReferenceDensity+...
 (Pressure-Parameters.ReferencePressure)*Parameters.CompressibilityCoeff;
Parameters.DDensityDPressure=@(Pressure)...      % dDensity/dPressure
  Parameters.CompressibilityCoeff;
clear r0 p0 eps mu L R U vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_hartmann_3'};                   % Mesh file
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
   'Pressure';
   'VelocityPost'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2-2.804100310651381e-01)<1e-12 && ',...
   'abs(Results.VelocityErrorL2        -1.996947872130190e-02)<1e-12 && ',...
   'abs(Results.PressureErrorL2        -2.630098602304956e-02)<1e-12 && ',...
   'abs(Results.VelocityPostErrorL2    -3.111002104775944e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main