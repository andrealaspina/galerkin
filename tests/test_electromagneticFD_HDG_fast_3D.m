clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Electromagnetic';            % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='ElectromagneticFD_HDG_fast';% Formulation
Parameters.Problem='Electromagnetic';            % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=2;                             % Degree
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,w) [(1-1i)*(-y/2+z/2),...
              (1-1i)*(-x/1-z/2),...
              (1-1i)*(+x/1+y/2)];
Parameters.ElectricField=...                     % Electric field
  @(x,y,z,w) [(1+1i)*(+x.*x/2+x.*y/1+x.*z/1)*w,...
              (1+1i)*(-x.*y/2-y.*y/4-y.*z/2)*w,...
              (1+1i)*(-x.*z/2-y.*z/2-z.*z/4)*w];
Parameters.CurrentDensity=...                    % Current density
  @(x,y,z,w) [-(-1+1i)/1*((x.*y+x.*z+x.^2/2)*w^2+1),...
              +(-1+1i)/2*((x.*y+y.*z+y.^2/2)*w^2+1),...
              +(-1+1i)/2*((x.*z+y.*z+z.^2/2)*w^2+1)];
Parameters.IncidentField=...                     % Incident field for absorbing boundary conditions
 @(x,y,z,w) [0*x,...
             0*x,...
             0*x];
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_cube_1'};                        % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
System.IncrementalForm='no';                     % Incremental form
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
Time.Frequency=1/(2*pi);                         % Frequency of time-harmonic problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='gmres';                             % Type
Solver.Preconditioner='ilu';                     % Preconditioner
Solver.PreconditionerAllTimeSteps='no';          % Compute preconditioner at all time steps
Solver.PreconditionerAllIterations='no';         % Compute preconditioner at all iterations
Solver.Restart=[];                               % Restart
Solver.Tolerance=1e-12;                          % Tolerance
Solver.MaxIterations=100;                        % Maximum number of iterations
Solver.Equilibrate='no';                         % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3,4,5,6];              % Dirichlet portion
Boundaries.Absorbing=[];                         % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'MagneticField'    ,'L2';
   'ElectricField'    ,'L2';
   'ElectricField'    ,'Hcurl';
   'ElectricFieldPost','Hcurl'};
Options.Test=...                                 % Test
  ['Results.MagneticFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldErrorHcurl    <1e-12 && ',...
   'Results.ElectricFieldPostErrorHcurl<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main