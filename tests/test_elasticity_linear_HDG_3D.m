clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='ElasticityLinear_HDG';   % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.Model='LinearElasticity';             % Model
Parameters.Degree=2;                             % Degree
Parameters.StabDisplacement=10;                  % Stabilization for displacement
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
System.Nonlinear='no';                           % Nonlinear
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
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeQuantity=...                      % Compute quantity of interest
  ['[~,Node]=ismembertol([0,0,10],Mesh.Nodes'',1e-6,''ByRows'',true);',...
   'Results.TipDisplacement(Time.TimeStep,:)=',...
   'Block.SolutionLocal(Node,Sizes.NumVoigtComp+(1:Sizes.NumSpaceDim));'];
Options.Test=...                                 % Test
  'abs(Results.TipDisplacement(:,2)-(-9.309384053988653e-01))<1e-08';
% --------------------------------------------------------------------------------------------------

%% Main

main