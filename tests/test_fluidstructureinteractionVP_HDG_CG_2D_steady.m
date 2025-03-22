clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='FluidStructureInteraction';  % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='WeaklyCompressibleFlowVP_HDG';% Formulation
Parameters(1).Problem='Fluid';                   % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).ConvectiveFlow='yes';              % Convective flow
Parameters(1).ArbitraryLagrangianEulerian='yes'; % Arbitrary Lagrangian-Eulerian description
Parameters(1).Degree=4;                          % Degree
Parameters(1).StabVelocity=100;                  % Stabilization for velocity
Parameters(1).StabPressure=100;                  % Stabilization for pressure
Parameters(1).ReferenceDensity=1e3;              % Reference density
Parameters(1).ReferencePressure=0;               % Reference pressure
Parameters(1).CompressibilityCoeff=1e-5;         % Compressibility coefficient
Parameters(1).DynamicViscosity=1;                % Dynamic viscosity
Parameters(1).ScaledStrainRate=...               % Scaled strain rate
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(1).Velocity=...                       % Velocity
  @(x,y,z,t) [(abs(x)<1e-3)*3/2.*y.*(0.41-y)/(0.41/2)^2*0.2,...
              0*x];
Parameters(1).Pressure=@(x,y,z,t) 0*x;           % Pressure
Parameters(1).Traction=@(x,y,z,t) [0*x, 0*x];    % Traction
Parameters(1).ResidualContinuity=@(x,y,z,t) 0*x; % Residual of continuity equation
Parameters(1).Force=@(x,y,z,t) [0*x,0*x];        % Force
Parameters(1).Density=@(Pressure)...             % Equation of state
  Parameters(1).ReferenceDensity+...
 (Pressure-Parameters(1).ReferencePressure)*Parameters(1).CompressibilityCoeff;
Parameters(1).DDensityDPressure=@(Pressure)...   % dDensity/dPressure
  Parameters(1).CompressibilityCoeff;
Parameters(1).Displacement=@(x,y,z,t) [0*x,0*x]; % Displacement
Parameters(2).Formulation='Elasticity_CG';       % Formulation
Parameters(2).Problem='Mesh';                    % Problem
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Assumption='PlaneStrain';          % Assumption
Parameters(2).Degree=2;                          % Degree
Parameters(2).NitschePenalty=1e5;                % Nitsche's penalty parameter
Parameters(2).Density=0;                         % Density
Parameters(2).YoungsModulus=...                  % Young's modulus
  @(x,y,z) 1+(5-1)*max((1-abs(y-0.2)/0.2).^2,0)+(20-1)*((x-0.6).^2+(y-0.2).^2<0.03^2);
Parameters(2).PoissonsRatio=@(x,y,z) 0;          % Poisson's ratio
Parameters(2).Displacement=@(x,y,z,t) [0*x,0*x]; % Displacement
Parameters(2).Traction=...                       % Traction
  @(x,y,z,t,nx,ny,nz) [0*x,0*x];
Parameters(2).Force=@(x,y,z,t) [0*x,0*x];        % Force
Parameters(2).RelaxationParameter=0.8;           % Relaxation parameter
Parameters(3).Formulation='Elasticity_CG';       % Formulation
Parameters(3).Problem='Structural';              % Problem
Parameters(3).Model='StVenantKirchhoff';         % Model
Parameters(3).Assumption='PlaneStrain';          % Assumption
Parameters(3).Degree=4;                          % Degree
Parameters(3).NitschePenalty=1e3;                % Nitsche's penalty parameter
Parameters(3).Density=1e3;                       % Density
Parameters(3).YoungsModulus=@(x,y,z) 1.4e6;      % Young's modulus
Parameters(3).PoissonsRatio=@(x,y,z) 0.4;        % Poisson's ratio
Parameters(3).Displacement=@(x,y,z,t) [0*x,0*x]; % Displacement
Parameters(3).Traction=...                       % Traction
  @(x,y,z,t,nx,ny,nz) [0*x,0*x];
Parameters(3).Force=@(x,y,z,t) [0*x,0*x];        % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_turek_fsi_test'};               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-9;                           % Tolerance
System.MaxIterations=4;                          % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries(1).Dirichlet_v_x=[1:9,11,14];         % Dirichlet portion for x-velocity
Boundaries(1).Dirichlet_v_y=[1:9,11,14];         % Dirichlet portion for y-velocity
Boundaries(1).Dirichlet_p=15;                    % Dirichlet portion for pressure
Boundaries(1).Interface=[10,12,13];              % Interface portion
Boundaries(1).Neumann_t_x=[];                    % Neumann portion for x-traction
Boundaries(1).Neumann_t_y=[];                    % Neumann portion for y-traction
Boundaries(2).Dirichlet=[1:9,11,14,15];          % Dirichlet portion
Boundaries(2).Interface=[10,12,13];              % Interface portion
Boundaries(2).Neumann=[];                        % Neumann portion
Boundaries(3).Dirichlet=[];                      % Dirichlet portion
Boundaries(3).Interface=[2,3,4];                 % Interface portion
Boundaries(3).Neumann=[];                        % Neumann portion
Boundaries(3).Fixed=1;                           % Fixed portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure';
   'Displacement'};
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[~,Node]=ismembertol([0.6,0.2,0.0],Mesh(3).Nodes'',1e-6,''ByRows'',true); ',...
   'Results(3).DisplacementTip=Results(3).Displacement(Node,2);'];
Options.Test=...                                 % Test
  ['abs(Results(1).ScaledStrainRateErrorL2-2.677624263695627e-00)<1e-12 && ',...
   'abs(Results(1).VelocityErrorL2        -2.234157037107025e-01)<1e-12 && ',...
   'abs(Results(1).PressureErrorL2        -3.152898004748829e+01)<1e-11 && ',...
   'abs(Results(1).DisplacementErrorL2    -1.304474758703843e-04)<1e-12 && ',...
   'abs(Results(2).DisplacementErrorL2    -1.304475290827737e-04)<1e-12 && ',...
   'abs(Results(3).DisplacementErrorL2    -4.145099893799288e-05)<1e-12 && ',...
   'abs(Results(3).DisplacementTip        -9.538071176018161e-04)<1e-11'];
% --------------------------------------------------------------------------------------------------

%% Main

main