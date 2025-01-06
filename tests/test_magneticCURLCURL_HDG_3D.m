clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetic';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagneticCURLCURL_HDG';   % Formulation
Parameters.Problem='Magnetic';                   % Problem
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.ConvectiveTerm='yes';                 % Convective term
Parameters.Degree=2;                             % Degree
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
eta=Parameters.MagneticDiffusivity;
Parameters.ScaledMagneticCurl=...                % Scaled magnetic curl
  @(x,y,z,t) [y/2-z/2,...
              x+z/2,...
              -x-y/2];
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [+x.*x/2+x.*y/1+x.*z/1,...
              -x.*y/2-y.*y/4-y.*z/2,...
              -x.*z/2-y.*z/2-z.*z/4];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.Source=...                            % Source
  @(x,y,z,t) [-eta/1+2*x.*x+4*x.*y+4*x.*z,...
              +eta/2-2*y.*x-1*y.*y-2*y.*z,...
              +eta/2-2*z.*x-2*z.*y-1*z.*z];
Parameters.Velocity=@(x,y,z,t) [x,y,z];          % Velocity
Parameters.TangentialElectricField=...           % Tangential electric field
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters.NormalMagneticInduction=...           % Normal magnetic induction
  @(x,y,z,t) 0*x;
clear eta
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=multicuboid(1,1,1);                     % Geometry definition
Mesh.MinElementSize=1/1;                         % Minimum element size
Mesh.MaxElementSize=1/1;                         % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3,4,5,6];              % Dirichlet portion
Boundaries.Natural=[];                           % Natural portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError=...                         % Compute error
  {'ScaledMagneticCurl';
   'MagneticInduction';
   'LagrangeMultiplier'};
Options.Test=...                                 % Test
  ['Results.ScaledMagneticCurlErrorL2<1e-12 && ',...
   'Results.MagneticInductionErrorL2 <1e-12 && ',...
   'Results.LagrangeMultiplierErrorL2<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main