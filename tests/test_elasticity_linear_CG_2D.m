clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='LinearElasticity';             % Model
Parameters.Assumption='PlaneStrain';             % Assumption
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=2.5e3;                 % Nitsche's penalty parameter
Parameters.Density=1;                            % Density
Parameters.YoungsModulus=...                     % Young's modulus
  @(x,y,z) 25*     ((x>0)&(y<0)|(x<0)&(y>0))...
          +250*    ((x<0)&(y<0)|(x>0)&(y>0));
Parameters.PoissonsRatio=...                     % Poisson's ratio
  @(x,y,z) 0.49999*((x>0)&(y<0)|(x<0)&(y>0))...
          +0.3*    ((x<0)&(y<0)|(x>0)&(y>0));
E=Parameters.YoungsModulus; nu=Parameters.PoissonsRatio;
Parameters.Displacement=...                      % Displacement
  @(x,y,z,t) [1./E(x,y,z).*(sin(2.*pi.*y).*(-1+cos(2.*pi.*x)).*(2+2.*nu(x,y,z)))+x.*y.*sin(pi.*x).*sin(pi.*y).*((1+nu(x,y,z)).*(1-2.*nu(x,y,z)))./(1-nu(x,y,z)-2.*nu(x,y,z).^2+E(x,y,z).*nu(x,y,z)),...
              1./E(x,y,z).*(sin(2.*pi.*x).*(+1-cos(2.*pi.*y)).*(2+2.*nu(x,y,z)))+x.*y.*sin(pi.*x).*sin(pi.*y).*((1+nu(x,y,z)).*(1-2.*nu(x,y,z)))./(1-nu(x,y,z)-2.*nu(x,y,z).^2+E(x,y,z).*nu(x,y,z))];
Parameters.Traction=...                          % Traction
  @(x,y,z,t) [0*x,0*x];
Parameters.Force=...                             % Force
  @(x,y,z,t) [(24.*pi.^2.*sin(2.*pi.*y)-24.*pi.^2.*nu(x,y,z).*sin(2.*pi.*y)-48.*pi.^2.*nu(x,y,z).^2.*sin(2.*pi.*y)+E(x,y,z).*sin(pi.*x).*sin(pi.*y)-64.*pi.^2.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+24.*E(x,y,z).*pi.^2.*nu(x,y,z).*sin(2.*pi.*y)+pi.*E(x,y,z).*x.*cos(pi.*x).*sin(pi.*y)+2.*pi.*E(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)+4.*pi.*E(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+pi.*E(x,y,z).*y.*cos(pi.*y).*sin(pi.*x)+64.*pi.^2.*nu(x,y,z).*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+128.*pi.^2.*nu(x,y,z).^2.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)-4.*pi.*E(x,y,z).*nu(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+E(x,y,z).*pi.^2.*x.*y.*cos(pi.*x).*cos(pi.*y)-3.*E(x,y,z).*pi.^2.*x.*y.*sin(pi.*x).*sin(pi.*y)-64.*E(x,y,z).*pi.^2.*nu(x,y,z).*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+4.*E(x,y,z).*pi.^2.*nu(x,y,z).*x.*y.*sin(pi.*x).*sin(pi.*y))./(2.*nu(x,y,z)-2.*E(x,y,z).*nu(x,y,z)+4.*nu(x,y,z).^2-2),...
              (24.*pi.^2.*nu(x,y,z).*sin(2.*pi.*x)-24.*pi.^2.*sin(2.*pi.*x)+48.*pi.^2.*nu(x,y,z).^2.*sin(2.*pi.*x)+E(x,y,z).*sin(pi.*x).*sin(pi.*y)+64.*pi.^2.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-24.*E(x,y,z).*pi.^2.*nu(x,y,z).*sin(2.*pi.*x)+pi.*E(x,y,z).*x.*cos(pi.*x).*sin(pi.*y)+4.*pi.*E(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)+2.*pi.*E(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+pi.*E(x,y,z).*y.*cos(pi.*y).*sin(pi.*x)-64.*pi.^2.*nu(x,y,z).*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-128.*pi.^2.*nu(x,y,z).^2.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+E(x,y,z).*pi.^2.*x.*y.*cos(pi.*x).*cos(pi.*y)-3.*E(x,y,z).*pi.^2.*x.*y.*sin(pi.*x).*sin(pi.*y)+64.*E(x,y,z).*pi.^2.*nu(x,y,z).*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)+4.*E(x,y,z).*pi.^2.*nu(x,y,z).*x.*y.*sin(pi.*x).*sin(pi.*y))./(2.*nu(x,y,z)-2.*E(x,y,z).*nu(x,y,z)+4.*nu(x,y,z).^2-2)];
clear E nu
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_2'};                      % Mesh file
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
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError={'Displacement'};           % Compute error
Options.Test=...                                 % Test
  'abs(Results.DisplacementErrorL2-1.979147177254772e-01)<1e-12';
% --------------------------------------------------------------------------------------------------

%% Main

main