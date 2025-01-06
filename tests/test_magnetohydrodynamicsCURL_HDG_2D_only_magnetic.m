clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Magnetohydrodynamic';        % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='MagnetohydrodynamicsCURL_HDG';% Formulation
Parameters.Problem='Magnetohydrodynamic';        % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.ConvectiveFlow='no';                  % Convective flow
Parameters.CouplingTerms='no';                   % Coupling terms
Parameters.Degree=2;                             % Degree
Parameters.StabVelocity=1;                       % Stabilization for velocity
Parameters.StabPressure=1;                       % Stabilization for pressure
Parameters.StabMagneticInduction=1;              % Stabilization for magnetic induction
Parameters.StabLagrangeMultiplier=1;             % Stabilization for Lagrange multiplier
Parameters.ReferenceDensity=1;                   % Reference density
Parameters.ReferencePressure=0;                  % Reference pressure
Parameters.CompressibilityCoeff=0;               % Compressibility coefficient
Parameters.DynamicViscosity=1;                   % Dynamic viscosity
Parameters.MagneticDiffusivity=1;                % Magnetic diffusivity
Parameters.MagneticPermeability=1;               % Magnetic permeability
eta=Parameters.MagneticDiffusivity; mu_m=Parameters.MagneticPermeability;
Parameters.ScaledStrainRate=...                  % Scaled strain rate
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters.Velocity=@(x,y,z,t) [0*x,0*x];        % Velocity
Parameters.Pressure=@(x,y,z,t) 0*x;              % Pressure
Parameters.Traction=@(x,y,z,t) [0*x,0*x];        % Traction
Parameters.ResidualContinuity=@(x,y,z,t) 0*x;    % Residual of continuity equation
Parameters.Force=@(x,y,z,t) [0*x,0*x];           % Force
Parameters.Density=@(Pressure)...                % Equation of state
  Parameters.ReferenceDensity+...
 (Pressure-Parameters.ReferencePressure)*Parameters.CompressibilityCoeff;
Parameters.DDensityDPressure=@(Pressure)...      % dDensity/dPressure
  Parameters.CompressibilityCoeff;
Parameters.ScaledMagneticCurl=...                % Scaled magnetic curl
  @(x,y,z,t) -sqrt(mu_m)*(x-y);
Parameters.ScaledMagneticGradient=...            % Scaled magnetic gradient
  @(x,y,z,t) [-sqrt(2*eta)*(x-y),...
              +sqrt(2*eta)*(x-y),...
              +sqrt(eta)*(x+y)];
Parameters.MagneticInduction=...                 % Magnetic induction
  @(x,y,z,t) [x.^2/2-x.*y,...
              y.^2/2-x.*y];
Parameters.LagrangeMultiplier=@(x,y,z,t) 0*x;    % Lagrange multiplier
Parameters.PseudoTraction=...                    % Pseudo-traction
  @(x,y,z,t) [+((x<(-1+1e-6))|(x>(+1-1e-6))).*sign(x)*2*eta.*(x.^2-y.^2)-...
               ((y<(-1+1e-6))|(y>(+1-1e-6))).*sign(y)*4*eta.*x.*y,...
              -((x<(-1+1e-6))|(x>(+1-1e-6))).*sign(x)*4*eta.*x.*y-...
               ((y<(-1+1e-6))|(y>(+1-1e-6))).*sign(y)*2*eta.*(x.^2-y.^2)];
Parameters.Source=...                            % Source
  @(x,y,z,t) [-eta*(x==x),...
              -eta*(x==x)];
clear eta mu_m
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
Boundaries.Dirichlet_v_x=[1,2,3,4];              % Dirichlet portion for x-velocity
Boundaries.Dirichlet_v_y=[1,2,3,4];              % Dirichlet portion for y-velocity
Boundaries.Dirichlet_p=[1,2,3,4];                % Dirichlet portion for pressure
Boundaries.Neumann_t_x=[];                       % Neumann portion for x-traction
Boundaries.Neumann_t_y=[];                       % Neumann portion for y-traction
Boundaries.Dirichlet_b_x=[1,4];                  % Dirichlet portion for x-magnetic induction
Boundaries.Dirichlet_b_y=[2,3];                  % Dirichlet portion for y-magnetic induction
Boundaries.Dirichlet_q=[1,2,3,4];                % Dirichlet portion for lagrange multiplier
Boundaries.Neumann_s_x=[];                       % Neumann portion for x-pseudo-traction
Boundaries.Neumann_s_y=[];                       % Neumann portion for y-pseudo-traction
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure';
   'VelocityPost';
   'ScaledMagneticCurl';
   'ScaledMagneticGradient';
   'MagneticInduction';
   'LagrangeMultiplier';
   'MagneticInductionPost'};
Options.Test=...                                 % Test
  ['Results.ScaledStrainRateErrorL2      <1e-12 && ',...
   'Results.VelocityErrorL2              <1e-12 && ',...
   'Results.PressureErrorL2              <1e-11 && ',...
   'Results.VelocityPostErrorL2          <1e-12 && ',...
   'Results.ScaledMagneticCurlErrorL2    <1e-12 && ',...
   'Results.ScaledMagneticGradientErrorL2<1e-12 && ',...
   'Results.MagneticInductionErrorL2     <1e-12 && ',...
   'Results.LagrangeMultiplierErrorL2    <1e-12 && ',...
   'Results.MagneticInductionPostErrorL2 <1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main