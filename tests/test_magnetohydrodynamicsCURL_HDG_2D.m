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
Parameters.CouplingTerms='yes';                  % Coupling terms
Parameters.Degree=2;                             % Degree
Parameters.StabVelocity=10;                      % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.StabMagneticInduction=10;             % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=0;               % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ChannelLength=10;                     % Channel length
Parameters.ChannelHalfHeight=1;                  % Channel half height
Parameters.ChannelExitMeanVelocity=1;            % Channel exit mean velocity
Parameters.Hartmann=1;                           % Hartmann number
p0=Parameters.ReferencePressure;      mu=Parameters.DynamicViscosity;
eta=Parameters.MagneticDiffusivity;   mu_m=Parameters.MagneticPermeability;
L=Parameters.ChannelLength;           R=Parameters.ChannelHalfHeight;
U=Parameters.ChannelExitMeanVelocity; Ha=Parameters.Hartmann;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [0*x,...
              0*x,...
              3*sqrt(mu)*U/R*sinh(y/R*Ha)*sinh(Ha)^(-1)];
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [3/2*U*2*(1-cosh(y/R*Ha)/cosh(Ha))*(Ha*tanh(Ha))^(-1),...
              0*x];
Parameters.Pressure=...                          % Pressure
  @(x,y,z,t) p0+3*mu*L*U/R^2*(1-x/L)-9/2*mu*U^2/eta*(sinh(y/R*Ha)/sinh(Ha)-y/R).^2*Ha^(-2);
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [-((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*(p0+3*mu*L*U/R^2*(1-x/L)-...
                (9/2*mu*U^2/eta*(y*sinh(Ha)-R*sinh(y/R*Ha))^2)/(Ha^2*R^2*sinh(Ha)^2))-...
               ((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y).*3*mu*U/R*sinh(y/R*Ha)*sinh(Ha)^(-1),...
              -((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*3*mu*U/R*sinh(y/R*Ha)*sinh(Ha)^(-1)-...
               ((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y).*(p0+3*mu*L*U/R^2*(1-x/L)-...
                (9/2*mu*U^2/eta*(y*sinh(Ha)-R*sinh(y/R*Ha))^2)/(Ha^2*R^2*sinh(Ha)^2))];
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
Parameters.Density=@(Pressure)...                % Equation of state
  Parameters.ReferenceDensity+...
 (Pressure-Parameters.ReferencePressure)*Parameters.CompressibilityCoeff;
Parameters.DDensityDPressure=@(Pressure)...      % dDensity/dPressure
  Parameters.CompressibilityCoeff;
Parameters.ScaledMagneticCurl=...                % Scaled magnetic curl
  @(x,y,z,t) -3*sqrt(mu/eta)*U/R*(sinh(Ha)-Ha*cosh(y/R*Ha))*(Ha*sinh(Ha))^(-1);
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
clear p0 mu eta mu_m L R U Ha
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_hartmann_1'};                   % Mesh file
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
Boundaries.Dirichlet_v_x=[1,2,3,4];              % Dirichlet portion for x-velocity
Boundaries.Dirichlet_v_y=[1,2,3,4];              % Dirichlet portion for y-velocity
Boundaries.Dirichlet_p=[1,2,3,4];                % Dirichlet portion for pressure
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
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
  ['abs(Results.ScaledStrainRateErrorL2      -8.856351705708788e-02)<1e-12 && ',...
   'abs(Results.VelocityErrorL2              -1.842147791908669e-02)<1e-12 && ',...
   'abs(Results.PressureErrorL2              -1.431964579321530e-02)<1e-11 && ',...
   'abs(Results.VelocityPostErrorL2          -2.494448856388832e-02)<1e-12 && ',...
   'abs(Results.ScaledMagneticCurlErrorL2    -1.675249085068344e-01)<1e-12 && ',...
   'abs(Results.ScaledMagneticGradientErrorL2-1.681436687531524e-01)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorL2     -3.997925275479078e-02)<1e-12 && ',...
   'abs(Results.LagrangeMultiplierErrorL2    -1.471549360938696e-02)<1e-12 && ',...
   'abs(Results.MagneticInductionPostErrorL2 -6.245815199816839e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main