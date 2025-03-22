clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetic';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagneticCURLCURL_HDG';   % Formulation
Parameters.Problem='Magnetic';                   % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveTerm='yes';                 % Convective term
Parameters.Degree=2;                             % Degree
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
eta=Parameters.MagneticDiffusivity;
Parameters.ScaledMagneticCurl=...                % Scaled magnetic curl
  @(x,y,z,t) sqrt(eta)*(x-y);
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [x.^2/2-x.*y,...
              y.^2/2-x.*y];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.Source=...                            % Source
  @(x,y,z,t) [-eta+(3*x.^2)/2-3*x.*y,...
              -eta+(3*y.^2)/2-3*x.*y];
Parameters.Velocity=@(x,y,z,t) [x,y];            % Velocity
tol=1e-9;
nx=@(x) (abs(mean(x)-(+1))<tol)-(abs(mean(x)-(-1))<tol);
ny=@(y) (abs(mean(y)-(+1))<tol)-(abs(mean(y)-(-1))<tol);
Parameters.TangentialElectricField=...           % Tangential electric field
  @(x,y,z,t) [-3/2*x.*y.*(x-y).*ny(y)-eta*(x-y).*ny(y),...
              +3/2*x.*y.*(x-y).*nx(x)+eta*(x-y).*nx(x)];
Parameters.NormalMagneticInduction=...           % Normal magnetic induction
  @(x,y,z,t) (x.^2/2-x.*y).*nx(x)+(y.^2/2-x.*y).*ny(y);
clear eta tol nx ny
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_1'};                     % Mesh file
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
Boundaries.Dirichlet=[1,2];                      % Dirichlet portion
Boundaries.Natural=[3,4];                        % Natural portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledMagneticCurl','L2';
   'MagneticInduction','L2';
   'LagrangeMultiplier','L2';
   'MagneticInductionPost','L2';
   'MagneticInductionPost','Hcurl'};
Options.Test=...                                 % Test
  ['Results.ScaledMagneticCurlErrorL2      <1e-12 && ',...
   'Results.MagneticInductionErrorL2       <1e-12 && ',...
   'Results.LagrangeMultiplierErrorL2      <1e-12 && ',...
   'Results.MagneticInductionPostErrorL2   <1e-12 && ',...
   'Results.MagneticInductionPostErrorHcurl<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main