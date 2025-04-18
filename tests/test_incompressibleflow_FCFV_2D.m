clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Fluid';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='IncompressibleFlow_FCFV';% Formulation
Parameters.Problem='Fluid';                      % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=0;                             % Degree
Parameters.StabVelocity=1;                       % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.Density=1;                            % Density
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.ChannelLength=10;                     % Channel length
Parameters.ChannelHalfHeight=1;                  % Channel half height
Parameters.ChannelExitMeanVelocity=1;            % Channel exit mean velocity
mu=Parameters.DynamicViscosity; L=Parameters.ChannelLength;
R=Parameters.ChannelHalfHeight; U=Parameters.ChannelExitMeanVelocity;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [0*x,...
              0*x,...
              3*U*mu^(1/2)*y/R^2];
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [3/2*U*(1-(y/R).^2),...
              0*x];
Parameters.Pressure=...                          % Pressure
  @(x,y,z,t) 3*(mu*L*U)/(R^2)*(1-x/L);
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [-((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*(-(3*L*U*mu*(x/L-1))/R^2)-...
               ((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y)    .*((3*U*y)/R^2)*mu,...
              -((x<(0+1e-6)) |(x>(L-1e-6))) .*sign(x-L/2).*((3*U*y)/R^2)*mu+...
               ((y<(-R+1e-6))|(y>(+R-1e-6))).*sign(y)    .*(+(3*L*U*mu*(x/L-1))/R^2)];
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
clear r mu L R U
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
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2-6.164713540975509e-00)<1e-12 && ',...
   'abs(Results.VelocityErrorL2        -3.876876460280624e-00)<1e-12 && ',...
   'abs(Results.PressureErrorL2        -6.299486269227092e+01)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main