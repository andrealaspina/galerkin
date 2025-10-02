clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='FluidALE';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='IncompressibleFlow_FCFV';% Formulation
Parameters(1).Problem='Fluid';                   % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).ConvectiveFlow='yes';              % Convective flow
Parameters(1).ArbitraryLagrangianEulerian='yes'; % Arbitrary Lagrangian-Eulerian description
Parameters(1).Degree=0;                          % Degree
Parameters(1).StabVelocity=1;                    % Stabilization for velocity
Parameters(1).StabPressure=1;                    % Stabilization for pressure
Parameters(1).Density=1;                         % Density
Parameters(1).DynamicViscosity=0.1;              % Dynamic viscosity
r=Parameters(1).Density; mu=Parameters(1).DynamicViscosity;
Parameters(1).ScaledStrainRate=...               % Scaled strain rate
  @(x,y,z,t) [-2*pi*cos(pi*x)*sin(pi*y)*sin(pi*t),...
              +2*pi*cos(pi*x)*sin(pi*y)*sin(pi*t),...
              0*x];
Parameters(1).Velocity=...                       % Velocity
  @(x,y,z,t) [sin(pi.*x).*sin(pi.*y).*sin(pi.*t),...
              cos(pi.*x).*cos(pi.*y).*sin(pi.*t)];
Parameters(1).Pressure=...                       % Pressure
  @(x,y,z,t) pi.*cos(pi.*x).*sin(pi.*y).*sin(pi.*t);
Parameters(1).Traction=...                       % Traction
  @(x,y,z,t) [0*x, 0*x];
Parameters(1).Force=...                          % Force
  @(x,y,z,t) [pi.*sin(pi.*x).*(+r.*cos(pi.*x)+r.*cos(pi.*t).*sin(pi.*y)-pi.*sin(pi.*t).*sin(pi.*y)-r.*cos(pi.*t)^2.*cos(pi.*x)+2.*mu.*pi.*sin(pi.*t).*sin(pi.*y)),...
              pi.*cos(pi.*y).*(-r.*sin(pi.*y)+r.*cos(pi.*t).*cos(pi.*x)+pi.*cos(pi.*x).*sin(pi.*t)+r.*cos(pi.*t)^2.*sin(pi.*y)+2.*mu.*pi.*cos(pi.*x).*sin(pi.*t))];
Parameters(1).ScaledStrainRateCenter=...         % Scaled strain rate at the center
  @(t) Parameters(1).ScaledStrainRate(0.5,0.5,0,t);
Parameters(1).VelocityCenter=...                 % Velocity at the center
  @(t) Parameters(1).Velocity(0.5,0.5,0,t);
Parameters(1).PressureCenter=...                 % Pressure at the center
  @(t) Parameters(1).Pressure(0.5,0.5,0,t);
Parameters(2).Formulation='Elasticity_CG';       % Formulation
Parameters(2).Problem='Mesh';                    % Problem
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Assumption='PlaneStrain';          % Assumption
Parameters(2).Degree=1;                          % Degree
Parameters(2).NitschePenalty=1000;               % Nitsche's penalty parameter
Parameters(2).Density=1;                         % Density
Parameters(2).YoungsModulus=@(x,y,z) 5/4;        % Young's modulus
Parameters(2).PoissonsRatio=@(x,y,z) 1/4;        % Poisson's ratio
Parameters(2).Displacement=...                   % Displacement
  @(x,y,z,t) [1/4*sin(2*pi*x).*(1-cos(2*pi*y)).*(1-cos(2*pi*t))*1/8,...
              1/4*sin(2*pi*y).*(1-cos(2*pi*x)).*(1-cos(2*pi*t))*1/8];
Parameters(2).Traction=@(x,y,z,t) [0*x,0*x];     % Traction
Parameters(2).Force=...                          % Force
  @(x,y,z,t) [-pi^2/16*sin(2*pi*x).*(6*cos(2*pi*y)+cos(2*pi*t)-4*cos(2*pi*y)*cos(2*pi*t)-3),...
              -pi^2/16*sin(2*pi*y).*(6*cos(2*pi*x)+cos(2*pi*t)-4*cos(2*pi*x)*cos(2*pi*t)-3)];
Parameters(2).DisplacementCenter=...             % Displacement at the center
  @(t) Parameters(2).Displacement(0.5,0.5,0,t);
clear r mu
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square01_ale_2'};               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-9;                           % Tolerance
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
Boundaries(1).Dirichlet_v_x=[1,2,3,4];           % Dirichlet portion for x-velocity
Boundaries(1).Dirichlet_v_y=[1,2,3,4];           % Dirichlet portion for y-velocity
Boundaries(1).Dirichlet_p=[1,2,3,4];             % Dirichlet portion for pressure
Boundaries(1).Neumann_t_x=[];                    % Neumann portion for x-traction
Boundaries(1).Neumann_t_y=[];                    % Neumann portion for y-traction
Boundaries(2).Dirichlet=[1,2,3,4];               % Dirichlet portion
Boundaries(2).Neumann=[];                        % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRateCenter','Number';
   'VelocityCenter'        ,'Number';
   'PressureCenter'        ,'Number';
   'DisplacementCenter'    ,'Number'};
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[~,Node1]=ismember([0.5,0.5,0],Mesh(1).Nodes'',''rows''); ',...
   '[~,Node2]=ismember([0.5,0.5,0],Mesh(2).Nodes'',''rows''); ',...
   'Results(1).ScaledStrainRateCenter=Results(1).ScaledStrainRate(Node1,:); ',...
   'Results(1).VelocityCenter=Results(1).Velocity(Node1,:); ',...
   'Results(1).PressureCenter=Results(1).Pressure(Node1,:); ',...
   'Results(2).DisplacementCenter=Results(2).Displacement(Node2,:);'];
Options.Test=...                                 % Test
  ['abs(Results(1).ScaledStrainRateCenterErrorNumber-2.683934577110453e-00)<1e-12 && ',...
   'abs(Results(1).VelocityCenterErrorNumber        -5.383729373024970e-01)<1e-12 && ',...
   'abs(Results(1).PressureCenterErrorNumber        -1.058816192322839e-00)<1e-12 && ',...
   'abs(Results(2).DisplacementCenterErrorNumber    -3.510466001877625e-17)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main