clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetohydrodynamic';        % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagnetohydrodynamicsCURL_HDG';% Formulation
Parameters.Problem='Magnetohydrodynamic';        % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveFlow='yes';                 % Convective flow
Parameters.CouplingTerms='no';                   % Coupling terms
Parameters.Degree=1;                             % Degree
Parameters.StabVelocity=10;                      % Stabilization for velocity
Parameters.StabPressure=10;                      % Stabilization for pressure
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=1/10;            % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.CylinderHeight=4;                     % Cylinder height
Parameters.InnerRadius=1;                        % Inner radius
Parameters.OuterRadius=2;                        % Outer radius
Parameters.InnerAngularVelocity=0;               % Inner wall angular velocity
Parameters.OuterAngularVelocity=1/2;             % Outer wall angular velocity
r0=Parameters.ReferenceDensity;      p0=Parameters.ReferencePressure;
eps=Parameters.CompressibilityCoeff; mu=Parameters.DynamicViscosity;
r1=Parameters.InnerRadius;           r2=Parameters.OuterRadius;
w1=Parameters.InnerAngularVelocity;  w2=Parameters.OuterAngularVelocity;
vx=@(x,y) +(y.*(-w2*r1^2*r2^2+w1*r1^2*r2^2+w2*r2^2*y.^2-...
                 w1*r1^2*y.^2+w2*r2^2*x.^2-w1*r1^2*x.^2))./((y.^2+x.^2)*(-r2^2+r1^2));
vy=@(x,y) -(x.*(-w2*r1^2*r2^2+w1*r1^2*r2^2+w2*r2^2*y.^2-...
                 w1*r1^2*y.^2+w2*r2^2*x.^2-w1*r1^2*x.^2))./((y.^2+x.^2)*(-r2^2+r1^2));
p= @(x,y) p0+(r0*(exp(eps*(((-w2*r2^2+w1*r1^2)^2*(-r1^2+y.^2+x.^2))/...
          (2*(r1^2-r2^2)^2)+(r1^2*r2^4*(-w2+w1)^2*(-r1^2+y.^2+x.^2))./...
          (2*(y.^2+x.^2)*(r1^2-r2^2)^2)-(2*r1^2*r2^2*(log(y.^2+x.^2)/2-...
          log(r1))*(-w2+w1)*(-w2*r2^2+w1*r1^2))/(r1^2-r2^2)^2))-1))/eps;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [+(2*2^(1/2)*mu^(1/2)*r1^2*r2^2*x.*y*(-w2+w1))./((x.^2+y.^2).^2*(-r2^2+r1^2)),...
              -(2*2^(1/2)*mu^(1/2)*r1^2*r2^2*x.*y*(-w2+w1))./((x.^2+y.^2).^2*(-r2^2+r1^2)),...
              -(2*mu^(1/2)*r1^2*r2^2*(-y.^2+x.^2)*(-w2+w1))./((x.^2+y.^2).^2*(-r2^2+r1^2))];
Parameters.Velocity=@(x,y,z,t) [vx(x,y),vy(x,y)];% Velocity
Parameters.Pressure=@(x,y,z,t) p(x,y);           % Pressure
Parameters.Traction=@(x,y,z,t) [0*x,0*x];        % Traction
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
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
clear r0 p0 eps mu r1 r2 w1 w2 vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_taylorcouette_1'};              % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-3;                           % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_v_x=1:16;                   % Dirichlet portion for x-velocity
Boundaries.Dirichlet_v_y=1:16;                   % Dirichlet portion for y-velocity
Boundaries.Dirichlet_p=1:16;                     % Dirichlet portion for pressure
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
Boundaries.Dirichlet_b_x=1:16;                   % Dirichlet portion for x-magnetic induction
Boundaries.Dirichlet_b_y=1:16;                   % Dirichlet portion for y-magnetic induction
Boundaries.Dirichlet_q=1:16;                     % Dirichlet portion for lagrange multiplier
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
   'ScaledMagneticCurl';
   'ScaledMagneticGradient';
   'MagneticInduction';
   'LagrangeMultiplier'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2      -3.729485507228658e-01)<1e-12 && ',...
   'abs(Results.VelocityErrorL2              -3.595960620937351e-02)<1e-12 && ',...
   'abs(Results.PressureErrorL2              -4.082658880469952e-02)<1e-12 && ',...
   'abs(Results.ScaledMagneticCurlErrorL2    -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ScaledMagneticGradientErrorL2-0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.MagneticInductionErrorL2     -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.LagrangeMultiplierErrorL2    -0.000000000000000e-00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main