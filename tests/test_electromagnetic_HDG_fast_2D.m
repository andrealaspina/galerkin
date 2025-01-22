clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Electromagnetic';            % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Electromagnetic_HDG_fast';% Formulation
Parameters.Problem='Electromagnetic';            % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=3;                             % Degree
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,t) 1/2*(x-y).*(x.*y-t^2);
Parameters.ElectricField=...                     % Electric field
  @(x,y,z,t) [(1/2*x.^2-x.*y)*t,...
              (1/2*y.^2-x.*y)*t];
Parameters.CurrentDensity=...                    % Current density
  @(x,y,z,t) [1/2*t^2*(x==x),...
              1/2*t^2*(x==x)];
tol=1e-9;
nx=@(x) (abs(mean(x)-(+1))<tol)-(abs(mean(x)-(-1))<tol);
ny=@(y) (abs(mean(y)-(+1))<tol)-(abs(mean(y)-(-1))<tol);
Parameters.IncidentField=...                     % Incident field for absorbing boundary conditions
 @(x,y,z,t) [+ny(y)*(ny(y)*t*(-x.^2/2+y.*x)-nx(x)*t*(-y.^2/2+x.*y))-ny(y)*(x/2-y/2).*(-t^2+x.*y),...
             -nx(x)*(ny(y)*t*(-x.^2/2+y.*x)-nx(x)*t*(-y.^2/2+x.*y))+nx(x)*(x/2-y/2).*(-t^2+x.*y)];
clear tol nx ny
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.SymmetrizeMatrixOnlyOnce='yes';           % Symmetrize matrix only once
System.Nonlinear='no';                           % Nonlinear
System.IncrementalForm='no';                     % Incremental form
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/2;                           % Time step size
Time.BDFOrder=3;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='pcg';                               % Type
Solver.Preconditioner='ichol';                   % Preconditioner
Solver.PreconditionerAllTimeSteps='no';          % Compute preconditioner at all time steps
Solver.PreconditionerAllIterations='no';         % Compute preconditioner at all iterations
Solver.Tolerance=1e-12;                          % Tolerance
Solver.MaxIterations=100;                        % Maximum number of iterations
Solver.Equilibrate='no';                         % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2];                      % Dirichlet portion
Boundaries.Absorbing=[3,4];                      % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'MagneticField','L2';
   'ElectricField','L2';
   'ElectricFieldPost','L2';
   'ElectricFieldPost','Hcurl'};
Options.Test=...                                 % Test
  ['Results.MagneticFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldPostErrorL2   <1e-12 && ',...
   'Results.ElectricFieldPostErrorHcurl<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main