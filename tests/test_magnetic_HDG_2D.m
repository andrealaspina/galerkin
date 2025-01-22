clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetic';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Magnetic_HDG';           % Formulation
Parameters.Problem='Magnetic';                   % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveTerm='yes';                 % Convective term
Parameters.Degree=2;                             % Degree
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
eta=Parameters.MagneticDiffusivity;
Parameters.ScaledMagneticGradient=...            % Scaled magnetic gradient
  @(x,y,z,t) [-sqrt(2*eta)*(x-y),...
              +sqrt(2*eta)*(x-y),...
              +sqrt(eta)*(x+y)];
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [x.^2/2-x.*y,...
              y.^2/2-x.*y];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.PseudoTraction=...                    % Pseudo-traction
  @(x,y,z,t) [+((x<(-1+1e-6))|(x>(+1-1e-6))).*sign(x)*2*eta.*(x-y)-...
               ((y<(-1+1e-6))|(y>(+1-1e-6))).*sign(y)*1*eta.*(x+y),...
              -((x<(-1+1e-6))|(x>(+1-1e-6))).*sign(x)*1*eta.*(x+y)-...
               ((y<(-1+1e-6))|(y>(+1-1e-6))).*sign(y)*2*eta.*(x-y)];
Parameters.Source=...                            % Source
  @(x,y,z,t) [-eta+(3*x.^2)/2-3*x.*y,...
              -eta+(3*y.^2)/2-3*x.*y];
Parameters.Velocity=@(x,y,z,t) [x,y];            % Velocity
clear eta
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
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
Boundaries.Dirichlet_b_x=[1,4];                  % Dirichlet portion for x-magnetic induction
Boundaries.Dirichlet_b_y=[2,3];                  % Dirichlet portion for y-magnetic induction
Boundaries.Dirichlet_q=[];                       % Dirichlet portion for lagrange multiplier
Boundaries.Neumann_s_x=[2,3];                    % Neumann portion for x-pseudo-traction
Boundaries.Neumann_s_y=[1,4];                    % Neumann portion for y-pseudo-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'ScaledMagneticGradient';
   'MagneticInduction';
   'LagrangeMultiplier';
   'MagneticInductionPost'};
Options.Test=...                                 % Test
  ['Results.ScaledMagneticGradientErrorL2<1e-12 && ',...
   'Results.MagneticInductionErrorL2     <1e-12 && ',...
   'Results.LagrangeMultiplierErrorL2    <1e-12 && ',...
   'Results.MagneticInductionPostErrorL2 <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main