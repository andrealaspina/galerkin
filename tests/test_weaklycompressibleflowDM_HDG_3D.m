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
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=1;                             % Degree
Parameters.StabDensity=100;                      % Stabilization for density
Parameters.StabMomentum=10;                      % Stabilization for momentum
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=1/10;            % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.CylinderHeight=4;                     % Cylinder height
Parameters.InnerRadius=1;                        % Inner radius
Parameters.OuterRadius=2;                        % Outer radius
Parameters.InnerAngularVelocity=0;               % Inner wall angular velocity
Parameters.OuterAngularVelocity=1/2;             % Outer wall angular velocity
r0=Parameters.ReferenceDensity;      p0=Parameters.ReferencePressure;
eps=Parameters.CompressibilityCoeff; mu=Parameters.DynamicViscosity;
H=Parameters.CylinderHeight;
r1=Parameters.InnerRadius;          r2=Parameters.OuterRadius;
w1=Parameters.InnerAngularVelocity; w2=Parameters.OuterAngularVelocity;
vx=@(x,y,z) +(y.*(-w2*r1^2*r2^2+w1*r1^2*r2^2+w2*r2^2*y.^2-...
                   w1*r1^2*y.^2+w2*r2^2*x.^2-w1*r1^2*x.^2))./((y.^2+x.^2)*(-r2^2+r1^2));
vy=@(x,y,z) -(x.*(-w2*r1^2*r2^2+w1*r1^2*r2^2+w2*r2^2*y.^2-...
                   w1*r1^2*y.^2+w2*r2^2*x.^2-w1*r1^2*x.^2))./((y.^2+x.^2)*(-r2^2+r1^2));
vz=@(x,y,z) 0*x;
p= @(x,y,z) p0+(r0*(exp(eps*(((-w2*r2^2+w1*r1^2)^2*(-r1^2+y.^2+x.^2))/...
            (2*(r1^2-r2^2)^2)+(r1^2*r2^4*(-w2+w1)^2*(-r1^2+y.^2+x.^2))./...
            (2*(y.^2+x.^2)*(r1^2-r2^2)^2)-(2*r1^2*r2^2*(log(y.^2+x.^2)/2-...
            log(r1))*(-w2+w1)*(-w2*r2^2+w1*r1^2))/(r1^2-r2^2)^2))-1))/eps;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [+(2*2^(1/2)*mu^(1/2)*r1^2*r2^2*x.*y*(-w2+w1))./((x.^2+y.^2).^2*(-r2^2+r1^2)),...
              -(2*2^(1/2)*mu^(1/2)*r1^2*r2^2*x.*y*(-w2+w1))./((x.^2+y.^2).^2*(-r2^2+r1^2)),...
              0*x,...
              -(2*mu^(1/2)*r1^2*r2^2*(-y.^2+x.^2)*(-w2+w1))./((x.^2+y.^2).^2*(-r2^2+r1^2)),...
              0*x,...
              0*x];
Parameters.Density=...                           % Density
  @(x,y,z,t) r0+eps*(p(x,y,z)-p0);
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [Parameters.Density(x,y,z,t).*vx(x,y,z),...
              Parameters.Density(x,y,z,t).*vy(x,y,z),...
              Parameters.Density(x,y,z,t).*vz(x,y,z)];
Parameters.Traction=@(x,y,z,t) [0*x,0*x,0*x];    % Traction
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
Parameters.Pressure=@(Density)...                % Equation of state
  Parameters.ReferencePressure+...
 (Density-Parameters.ReferenceDensity)/Parameters.CompressibilityCoeff;
Parameters.DPressureDDensity=@(Density)...       % dPressure/dDensity
  1/Parameters.CompressibilityCoeff;
clear r0 p0 eps mu H r1 r2 w1 w2 vx vy vz p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_taylorcouette3d_test'};         % Mesh file
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
Boundaries.Dirichlet_r=1:18;                     % Dirichlet portion for density
Boundaries.Dirichlet_w_x=1:18;                   % Dirichlet portion for x-momentum
Boundaries.Dirichlet_w_y=1:18;                   % Dirichlet portion for y-momentum
Boundaries.Dirichlet_w_z=1:18;                   % Dirichlet portion for z-momentum
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
Boundaries.Neumann_t_z=[];                       % Neumann portion for z-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Density';
   'Momentum'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2-7.711068927404321e-01)<1e-12 && ',...
   'abs(Results.DensityErrorL2         -1.590378056745913e-02)<1e-12 && ',...
   'abs(Results.MomentumErrorL2        -1.279396228741930e-01)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main