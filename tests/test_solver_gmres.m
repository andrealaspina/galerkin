clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Fluid';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='WeaklyCompressibleFlowVP_HDG';% Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=2;                             % Degree
Parameters.StabVelocity=1;                       % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=0;               % Compressibility coefficient
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
Parameters.Velocity=@(x,y,z,t) [vx(x,y),vy(x,y)];% Velocity
Parameters.Pressure=@(x,y,z,t) p(x,y);           % Pressure
Parameters.Traction=@(x,y,z,t) [0*x, 0*x];       % Traction
Parameters.ResidualContinuity=...                % Residual of continuity equation
  @(x,y,z,t) (27*mu^2*U^3)/(4*r0*R^8)*(R^2-y.^2).*...
             (6*R^2*(L-x)-mu*U/r0*(9*(L-x).^2-R^2+y.^2)*eps)*eps^2;
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
Parameters.Density=@(Pressure)...                % Equation of state
  Parameters.ReferenceDensity+...
 (Pressure-Parameters.ReferencePressure)*Parameters.CompressibilityCoeff;
Parameters.DDensityDPressure=@(Pressure)...      % dDensity/dPressure
  Parameters.CompressibilityCoeff;
clear r0 p0 eps mu L R U vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_poiseuille_1'};                    % Mesh file
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
Solver.Type='gmres';                             % Type
Solver.Preconditioner='ilu';                     % Preconditioner
Solver.PreconditionerAllTimeSteps='no';          % Compute preconditioner at all time steps
Solver.PreconditionerAllIterations='no';         % Compute preconditioner at all iterations
Solver.Restart=[];                               % Restart
Solver.Tolerance=1e-3;                           % Tolerance
Solver.MaxIterations=100;                        % Maximum number of iterations
Solver.Equilibrate='no';                         % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_v_x=[1,3,4,4];              % Dirichlet portion for x-velocity
Boundaries.Dirichlet_v_y=[1,2,3,4];              % Dirichlet portion for y-velocity
Boundaries.Dirichlet_p=[1,2,3,4];                % Dirichlet portion for pressure
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure'};
Options.Test=...                                 % Test
  ['Results.ScaledStrainRateErrorL2<1e-12 && ',...
   'Results.VelocityErrorL2        <1e-12 && ',...
   'Results.PressureErrorL2        <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main