clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ConvergenceTime';               % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='StVenantKirchhoff';            % Model
Parameters.Assumption='PlaneStrain';             % Assumption
Parameters.PredictorDegree=3;                    % Predictor degree
Parameters.Degree=2;                             % Degree
Parameters.NitschePenalty=10^5;                  % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.YoungsModulus=@(x,y,z) 5/4;           % Young's modulus
Parameters.PoissonsRatio=@(x,y,z) 1/4;           % Poisson's ratio
Parameters.Displacement=...                      % Displacement
  @(x,y,z,t) [(x.^2+y.^2)*t^3,...
              (x.^2+y.^2)*t^3];
Parameters.Traction=...                          % Traction
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x];
Parameters.Force=...                             % Force
  @(x,y,z,t) [-2*t*(24*t^8*x.^2+24*t^8*y.^2+16*t^5*x+8*t^5*y+2*t^2-3*x.^2-3*y.^2),...
              -2*t*(24*t^8*x.^2+24*t^8*y.^2+16*t^5*y+8*t^5*x+2*t^2-3*x.^2-3*y.^2)];
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_1'};                     % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-10;                          % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/4;                           % Time step size
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
  ['Results.DisplacementErrorL2<1e-12 && ',...
   'System.Iteration==1'];
% --------------------------------------------------------------------------------------------------

%% Main

main