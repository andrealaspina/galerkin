clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='NeoHook';                      % Model
Parameters.Assumption='PlaneStrain';             % Assumption
Parameters.Degree=2;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=1e3;                          % Density
Parameters.YoungsModulus=@(x,y,z) 1.4e6;         % Young's modulus
Parameters.PoissonsRatio=@(x,y,z) 0.4;           % Poisson's ratio
Parameters.Displacement=@(x,y,z,t) [0*x,0*x];    % Displacement
Parameters.Velocity=@(x,y,z,t)[0*x,0*x];         % Velocity
Parameters.Traction=@(x,y,z,t)[0*x,0*x];         % Traction
Parameters.Force=...                             % Force
  @(x,y,z,t) [0*x,-Parameters.Density*2*(x==x)];
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_turek_structure'};               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
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
Boundaries.Dirichlet=2;                          % Dirichlet portion
Boundaries.Neumann=[1,3,4];                      % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeQuantity=...                      % Compute quantity of interest
  ['[~,Node]=ismember([0.6,0.2,0.0],Mesh.Nodes'',''rows'');',...
   'Results.TipDisplacement(Time.TimeStep,:)=',...
   'Block.SolutionGlobal(Node,:);'];
Options.Test=...                                 % Test
  'abs(Results.TipDisplacement(:,2)-(-6.560967263936134e-02))<1e-12';
% --------------------------------------------------------------------------------------------------

%% Main

main