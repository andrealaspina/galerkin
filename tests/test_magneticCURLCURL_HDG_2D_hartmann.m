clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetic';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagneticCURLCURL_HDG';   % Formulation
Parameters.Problem='Magnetic';                   % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveTerm='yes';                 % Convective term
Parameters.Degree=1;                             % Degree
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
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
Parameters.ScaledMagneticCurl=...                % Scaled magnetic curl
  @(x,y,z,t) 3*sqrt(eta)*sqrt(mu*mu_m/eta)*U*(sinh(Ha)-Ha*cosh(Ha*y/R))/(R*Ha*sinh(Ha));
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [3*sqrt(mu*mu_m/eta)*U*(sinh(y/R*Ha)/sinh(Ha)-y/R)*Ha^(-1),...
                sqrt(mu*mu_m*eta)/R*Ha*(x==x)];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.TangentialElectricField=...           % Tangential electric field
  @(x,y,z,t) [0*x,...
              0*x];
Parameters.NormalMagneticInduction=...           % Normal magnetic induction
  @(x,y,z,t) 0*x;
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
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Natural=[];                           % Natural portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError=...                         % Compute error
  {'ScaledMagneticCurl'   ,'L2';
   'MagneticInduction'    ,'L2';
   'LagrangeMultiplier'   ,'L2';
   'MagneticInduction'    ,'Hcurl';
   'MagneticInductionPost','Hcurl'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledMagneticCurlErrorL2      -1.994120349867612e-02)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorL2       -4.010428990955935e-02)<1e-12 && ',...
   'abs(Results.LagrangeMultiplierErrorL2      -6.678000216524408e-03)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorHcurl    -5.205065529430997e-01)<1e-12 && ',...
   'abs(Results.MagneticInductionPostErrorHcurl-3.784669013863075e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main