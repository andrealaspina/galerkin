clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ConvergenceTime';               % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='LinearElasticity';             % Model
Parameters.Assumption='PlaneStrain';             % Assumption
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=10^5;                  % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.YoungsModulus=@(x,y,z) 5/4;           % Young's modulus
Parameters.PoissonsRatio=@(x,y,z) 1/4;           % Poisson's ratio
Parameters.Displacement=...                      % Displacement
  @(x,y,z,t) [(x+y)*t^6,(x+y)*t^6];
Parameters.Traction=...                          % Traction
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x];
Parameters.Force=...                             % Force
  @(x,y,z,t) [(x+y)*6*5*t^4,(x+y)*6*5*t^4];
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_1'};                     % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1./[3,4];                      % Time step size
Time.BDFOrder=4;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.PlotConvergence='no';                    % Plot convergence
Options.CharacteristicElementSize='Min';         % Characteristic element size
Options.ComputeError={'Displacement'};           % Compute error
Options.Test=...                                 % Test
  'abs(Results.DisplacementErrorL2Rate(end)-Time.BDFOrder)<1e-2';
% --------------------------------------------------------------------------------------------------

%% Main

main