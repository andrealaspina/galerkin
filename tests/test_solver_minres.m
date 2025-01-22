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
Parameters.Assumption='PlaneStrain';             % Assumption
Parameters.Degree=1;                             % Degree
Parameters.StabDisplacement=1;                   % Stabilization for displacement
Parameters.Density=0;                            % Density
Parameters.YoungsModulus=@(x,y,z) 5/4;           % Young's modulus
Parameters.PoissonsRatio=@(x,y,z) 1/4;           % Poisson's ratio
Parameters.Displacement=@(x,y,z,t) [x+y,x+y];    % Displacement
Parameters.Traction=@(x,y,z,t) [0*x,0*x];        % Traction
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='minres';                            % Type
Solver.Preconditioner='ichol';                   % Preconditioner
Solver.PreconditionerAllTimeSteps='no';          % Compute preconditioner at all time steps
Solver.PreconditionerAllIterations='no';         % Compute preconditioner at all iterations
Solver.Tolerance=1e-6;                           % Tolerance
Solver.MaxIterations=100;                        % Maximum number of iterations
Solver.Equilibrate='no';                         % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError={'Displacement'};           % Compute error
Options.Test=...                                 % Test
  'Results.DisplacementErrorL2<5e-7';
% --------------------------------------------------------------------------------------------------

%% Main

main