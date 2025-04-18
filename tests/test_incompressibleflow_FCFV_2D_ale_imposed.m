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
Parameters.ArbitraryLagrangianEulerian='yes';    % Arbitrary Lagrangian-Eulerian description
Parameters.Degree=0;                             % Degree
Parameters.StabVelocity=1;                       % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.Density=1;                            % Density
Parameters.DynamicViscosity=0.1;                 % Dynamic viscosity
r=Parameters.Density; mu=Parameters.DynamicViscosity;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [-2^(1/2)*mu^(1/2)*pi*cos(pi*x)*sin(pi*t)*sin(pi*y),...
              +2^(1/2)*mu^(1/2)*pi*cos(pi*x)*sin(pi*t)*sin(pi*y),...
              0*x];
Parameters.Velocity=...                          % Velocity
  @(x,y,z,t) [sin(pi.*t).*sin(pi.*x).*sin(pi.*y),...
              cos(pi.*x).*cos(pi.*y).*sin(pi.*t)];
Parameters.Pressure=...                          % Pressure
  @(x,y,z,t) pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y);
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [0*x, 0*x];
Parameters.Force=...                             % Force
  @(x,y,z,t) [pi.*sin(pi.*x).*(+r.*cos(pi.*x)+r.*cos(pi.*t).*sin(pi.*y)-pi.*sin(pi.*t).*sin(pi.*y)-r.*cos(pi.*t)^2.*cos(pi.*x)+2.*mu.*pi.*sin(pi.*t).*sin(pi.*y)),...
              pi.*cos(pi.*y).*(-r.*sin(pi.*y)+r.*cos(pi.*t).*cos(pi.*x)+pi.*cos(pi.*x).*sin(pi.*t)+r.*cos(pi.*t)^2.*sin(pi.*y)+2.*mu.*pi.*cos(pi.*x).*sin(pi.*t))];
Parameters.Displacement=...                      % Displacement
  @(x,y,z,t) [1/4*sin(2*pi*x).*(1-cos(2*pi*y)).*(1-cos(2*pi*t))*1/8,...
              1/4*sin(2*pi*y).*(1-cos(2*pi*x)).*(1-cos(2*pi*t))*1/8];
Parameters.ScaledStrainRateCenter=...            % Scaled strain rate at the center
  @(t) Parameters.ScaledStrainRate(0.5,0.5,0,t);
Parameters.VelocityCenter=...                    % Velocity at the center
  @(t) Parameters.Velocity(0.5,0.5,0,t);
Parameters.PressureCenter=...                    % Pressure at the center
  @(t) Parameters.Pressure(0.5,0.5,0,t);
clear r mu
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square01_2'};                   % Mesh file
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
  {'ScaledStrainRateCenter','Number';
   'VelocityCenter'        ,'Number';
   'PressureCenter'        ,'Number'};
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[~,Node]=ismember([0.5,0.5,0],Mesh.Nodes'',''rows''); ',...
   'Results.ScaledStrainRateCenter=Results.ScaledStrainRate(Node,:); ',...
   'Results.VelocityCenter=Results.Velocity(Node,:); ',...
   'Results.PressureCenter=Results.Pressure(Node,:);'];
Options.Test=...                                 % Test
  ['abs(Results.ScaledStrainRateCenterErrorNumber-6.082459475196617e-01)<1e-12 && ',...
   'abs(Results.VelocityCenterErrorNumber        -5.692053686573615e-01)<1e-12 && ',...
   'abs(Results.PressureCenterErrorNumber        -1.080293813337339e-00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main