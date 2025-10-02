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
Parameters.ConvectiveFlow='yes';                 % Convective flow
Parameters.ArbitraryLagrangianEulerian='no';     % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=0;                             % Degree
Parameters.StabVelocity=1;                       % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.Density=1;                            % Density
Parameters.DynamicViscosity=0.1;                 % Dynamic viscosity
r=Parameters.Density; mu=Parameters.DynamicViscosity;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [-2*pi*cos(pi*x)*sin(pi*y)*sin(pi*t),...
              +2*pi*cos(pi*x)*sin(pi*y)*sin(pi*t),...
              0*x];
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [sin(pi.*x).*sin(pi.*y).*sin(pi.*t),...
              cos(pi.*x).*cos(pi.*y).*sin(pi.*t)];
Parameters.Pressure=...                          % Pressure
  @(x,y,z,t) pi.*cos(pi.*x).*sin(pi.*y).*sin(pi.*t);
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [0*x, 0*x];
Parameters.Force=...                             % Force
  @(x,y,z,t) [pi.*sin(pi.*x).*(+r.*cos(pi.*x)+r.*cos(pi.*t).*sin(pi.*y)-pi.*sin(pi.*t).*sin(pi.*y)-r.*cos(pi.*t)^2.*cos(pi.*x)+2.*mu.*pi.*sin(pi.*t).*sin(pi.*y)),...
              pi.*cos(pi.*y).*(-r.*sin(pi.*y)+r.*cos(pi.*t).*cos(pi.*x)+pi.*cos(pi.*x).*sin(pi.*t)+r.*cos(pi.*t)^2.*sin(pi.*y)+2.*mu.*pi.*cos(pi.*x).*sin(pi.*t))];
clear r mu
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square01_1'};                   % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-10;                          % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=0.5;                              % Final time
Time.TimeStepSize=2^(-3);                        % Time step size
Time.BDFOrder=4;                                 % BDF order
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
   'Pressure'};
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateErrorL2-2.423462777444449e-00)<1e-12 && ',...
   'abs(Results.VelocityErrorL2        -2.532299856794472e-01)<1e-12 && ',...
   'abs(Results.PressureErrorL2        -5.193011528487051e-01)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main