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
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.YoungsModulus=...                     % Young's modulus
  @(x,y,z) 1*      ((y>=-1.0)&(y<=-0.5)|(y>=+0.0)&(y<=+0.5))...
          +10*     ((y>=-0.5)&(y<=+0.0)|(y>=+0.5)&(y<=+1.0));
Parameters.PoissonsRatio=...                     % Poisson's ratio
  @(x,y,z) 0.49999*((y>=-1.0)&(y<=-0.5)|(y>=+0.0)&(y<=+0.5))...
          +0.3*    ((y>=-0.5)&(y<=+0.0)|(y>=+0.5)&(y<=+1.0));
Parameters.Displacement=@(x,y,z,t) [0*x,0*x,0*x];% Displacement
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [0*x,...
              -13/6000*(y>=(+1-1e-6)),...
              0*x];
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_beam_test'};                     % Mesh file
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
Boundaries.Neumann=[1,3,4,5,6];                  % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeQuantity=...                      % Compute quantity of interest
  ['[~,Node]=ismembertol([0,0,10],Mesh.Nodes'',1e-6,''ByRows'',true);',...
   'Results.TipDisplacement(Time.TimeStep,:)=',...
   'Block.SolutionGlobal(Node,:);'];
Options.Test=...                                 % Test
  'abs(Results.TipDisplacement(:,2)-(-2.195728817356443e-01))<1e-12';
% --------------------------------------------------------------------------------------------------

%% Main

main