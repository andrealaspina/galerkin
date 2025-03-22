clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='ElasticityLinear_HDG';% Formulation
Parameters(1).Problem='Structural';              % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).Model='LinearElasticity';          % Model
Parameters(1).Degree=2;                          % Degree
Parameters(1).StabDisplacement=10;               % Stabilization for displacement
Parameters(1).Density=1;                         % Density
Parameters(1).YoungsModulus=@(x,y,z) 1;          % Young's modulus
Parameters(1).PoissonsRatio=@(x,y,z) 0.49999;    % Poisson's ratio
Parameters(1).Displacement=...                   % Displacement
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(1).Traction=...                       % Traction
  @(x,y,z,t,nx,ny,nz) [0*x,-13/6000*(y>=(+1-1e-6)),0*x];
Parameters(1).Force=@(x,y,z,t) [0*x,0*x,0*x];    % Force
Parameters(2).Formulation='Elasticity_CG';       % Formulation
Parameters(2).Problem='Structural';              % Problem
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Degree=2;                          % Degree
Parameters(2).NitschePenalty=100;                % Nitsche's penalty parameter
Parameters(2).Density=1;                         % Density
Parameters(2).YoungsModulus=@(x,y,z) 10;         % Young's modulus
Parameters(2).PoissonsRatio=@(x,y,z) 0.3;        % Poisson's ratio
Parameters(2).Displacement=...                   % Displacement
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters(2).Traction=...                       % Traction
  @(x,y,z,t,nx,ny,nz) [0*x,-13/6000*(y>=(+1-1e-6)),0*x];
Parameters(2).Force=@(x,y,z,t) [0*x,0*x,0*x];    % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_beam_test_coupled'};            % Mesh file
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
Boundaries(1).Dirichlet=[2,8];                   % Dirichlet portion
Boundaries(1).Interface=[5,7,11];                % Interface portion
Boundaries(1).Neumann=[1,3,4,6,9,10,12];         % Neumann portion
Boundaries(2).Dirichlet=[2,8];                   % Dirichlet portion
Boundaries(2).Interface=[1,5,7];                 % Interface portion
Boundaries(2).Neumann=[3,4,6,9,10,11,12];        % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeQuantity=...                      % Compute quantity of interest
  ['[~,Node]=ismembertol([0,0,10],Mesh(1).Nodes'',1e-6,''ByRows'',true);',...
   'Results(1).TipDisplacement(Time.TimeStep,:)=',...
   'Block(1,1).SolutionLocal(Node,Sizes(1).NumVoigtComp+(1:Sizes(1).NumSpaceDim));'];
Options.Test=...                                 % Test
  'abs(Results(1).TipDisplacement(:,2)-(-9.098889620715099e-01))<1e-07';
% --------------------------------------------------------------------------------------------------

%% Main

main