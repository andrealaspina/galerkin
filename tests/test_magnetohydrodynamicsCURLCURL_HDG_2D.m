clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetohydrodynamic';        % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagnetohydrodynamicsCURLCURL_HDG';% Formulation
Parameters.Problem='Magnetohydrodynamic';        % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.CouplingTerms='yes';                  % Coupling terms
Parameters.Degree=2;                             % Degree
Parameters.StabVelocity=10;                      % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
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
  @(x,y,z,t) 3*sqrt(eta)*sqrt(mu*mu_m/eta)*U*(sinh(Ha)-Ha*cosh(Ha*y/R))/(R*Ha*sinh(Ha));
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [3*sqrt(mu*mu_m/eta)*U*(sinh(y/R*Ha)/sinh(Ha)-y/R)*Ha^(-1),...
                sqrt(mu*mu_m*eta)/R*Ha*(x==x)];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.TangentialElectricField=...           % Tangential electric field
  @(x,y,z,t) [0*x,0*x];
Parameters.NormalMagneticInduction=...           % Normal magnetic induction
  @(x,y,z,t) 0*x;
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
Boundaries.Dirichlet_m=[1,2,3,4];                % Dirichlet portion for magnetic field
Boundaries.Natural_m=[];                         % Natural portion for magnetic field
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate'     ,'L2';
   'Velocity'             ,'L2';
   'Pressure'             ,'L2';
   'VelocityPost'         ,'L2';
   'ScaledMagneticCurl'   ,'L2';
   'MagneticInduction'    ,'L2';
   'LagrangeMultiplier'   ,'L2';
   'MagneticInduction'    ,'Hcurl';
   'MagneticInductionPost','Hcurl'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2        -8.482238455657563e-02)<1e-12 && ',...
   'abs(Results.VelocityErrorL2                -1.773669472381909e-02)<1e-12 && ',...
   'abs(Results.PressureErrorL2                -2.169872843438172e-02)<1e-12 && ',...
   'abs(Results.VelocityPostErrorL2            -2.442944420362691e-02)<1e-12 && ',...
   'abs(Results.ScaledMagneticCurlErrorL2      -2.954268464258091e-02)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorL2       -3.262227427953723e-02)<1e-12 && ',...
   'abs(Results.LagrangeMultiplierErrorL2      -4.214607945392719e-03)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorHcurl    -8.461152946735624e-02)<1e-12 && ',...
   'abs(Results.MagneticInductionPostErrorHcurl-3.024205556502943e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main