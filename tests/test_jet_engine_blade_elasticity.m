clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='StVenantKirchhoff';            % Model
Parameters.Degree=2;                             % Degree
Parameters.NitschePenalty=1e12;                  % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.YoungsModulus=227e9;                  % Young's modulus
Parameters.PoissonsRatio=0.27;                   % Poisson's ratio
Parameters.Displacement=@(x,y,z,t) [0*x,0*x,0*x];% Displacement
Parameters.Traction=...                          % Traction
  @(x,y,z,t,b,nx,ny,nz) [-5.0e5*nx.*(b==11)-4.5e5*nx.*(b==10),...
                         -5.0e5*ny.*(b==11)-4.5e5*ny.*(b==10),...
                         -5.0e5*nz.*(b==11)-4.5e5*nz.*(b==10)];
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=importGeometry(createpde,'Blade.stl');  % Geometry definition
Geometry.translate(-min(Geometry.Vertices));     % Translate geometry
Mesh.MinElementSize=0.02;                        % Minimum element size
Mesh.MaxElementSize=0.02;                        % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e+12;                          % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=3;                          % Dirichlet portion
Boundaries.Neumann=[10,11];                      % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Test=...                                 % Test
  ['abs(max(Results.Displacement(:,1))-6.658784409406569e-10)<1e-12 && ',...
   'abs(max(Results.Displacement(:,2))-1.036833692703720e-05)<1e-12 && ',...
   'abs(max(Results.Displacement(:,3))-1.508248922949512e-04)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main
