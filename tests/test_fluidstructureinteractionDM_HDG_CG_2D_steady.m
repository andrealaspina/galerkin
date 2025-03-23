clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='FluidStructureInteraction';  % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='WeaklyCompressibleFlowDM_HDG';% Formulation
Parameters(1).Problem='Fluid';                   % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).ConvectiveFlow='yes';              % Convective flow
Parameters(1).ArbitraryLagrangianEulerian='yes'; % Arbitrary Lagrangian-Eulerian description
Parameters(1).Degree=4;                          % Degree
Parameters(1).StabDensity=100/1e-5;              % Stabilization for density
Parameters(1).StabMomentum=1;                    % Stabilization for momentum
Parameters(1).ReferenceDensity=1e3;              % Reference density
Parameters(1).ReferencePressure=0;               % Reference pressure
Parameters(1).CompressibilityCoeff=1e-5;         % Compressibility coefficient
Parameters(1).DynamicViscosity=1;                % Dynamic viscosity
Parameters(1).ScaledStrainRate=...               % Scaled strain rate
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(1).Density=@(x,y,z,t) 1e3*(x==x);     % Density
Parameters(1).Momentum=...                       % Momentum
  @(x,y,z,t) [(abs(x)<1e-3)*3/2.*y.*(0.41-y)/(0.41/2)^2*1e3*0.2,...
              0*x];
Parameters(1).Traction=@(x,y,z,t) [0*x,0*x];     % Traction
Parameters(1).ResidualContinuity=@(x,y,z,t) 0*x; % Residual of continuity equation
Parameters(1).Force=@(x,y,z,t) [0*x,0*x];        % Force
Parameters(1).Pressure=@(Density)...             % Equation of state
  Parameters(1).ReferencePressure+...
 (Density-Parameters(1).ReferenceDensity)/Parameters(1).CompressibilityCoeff;
Parameters(1).DPressureDDensity=@(Density)...    % dPressure/dDensity
  1/Parameters(1).CompressibilityCoeff;
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
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x];
Parameters(2).Force=@(x,y,z,t) [0*x,0*x];        % Force
Parameters(2).RelaxationParameter=0.7;           % Relaxation parameter
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
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x];
Parameters(3).Force=@(x,y,z,t) [0*x,0*x];        % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_turek_fsi_test'};               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-4;                           % Tolerance
System.MaxIterations=4;                          % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries(1).Dirichlet_r=15;                    % Dirichlet portion for density
Boundaries(1).Dirichlet_w_x=[1:9,11,14];         % Dirichlet portion for x-momentum
Boundaries(1).Dirichlet_w_y=[1:9,11,14];         % Dirichlet portion for y-momentum
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
   'Density';
   'Momentum';
   'Displacement'};
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['[~,Node]=ismembertol([0.6,0.2,0.0],Mesh(3).Nodes'',1e-6,''ByRows'',true); ',...
   'Results(3).DisplacementTip=Results(3).Displacement(Node,2);'];
Options.Test=...                                 % Test
  ['abs(Results(1).ScaledStrainRateErrorL2-2.678190046210910e-00)<1e-07 && ',...
   'abs(Results(1).DensityErrorL2         -3.149155884565349e-04)<1e-10 && ',...
   'abs(Results(1).MomentumErrorL2        -2.234262193053535e+02)<1e-05 && ',...
   'abs(Results(1).DisplacementErrorL2    -1.082396567008379e-04)<1e-09 && ',...
   'abs(Results(2).DisplacementErrorL2    -1.082397013381004e-04)<1e-09 && ',...
   'abs(Results(3).DisplacementErrorL2    -2.924802482496271e-05)<1e-09 && ',...
   'abs(Results(3).DisplacementTip        -6.710736922334208e-04)<1e-08'];
% --------------------------------------------------------------------------------------------------

%% Main

main