clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='Elasticity_CG';       % Formulation
Parameters(1).Problem='Structural';              % Problem
Parameters(1).Model='LinearElasticity';          % Model
Parameters(1).Assumption='PlaneStrain';          % Assumption
Parameters(1).Degree=2;                          % Degree
Parameters(1).NitschePenalty=2.5e3;              % Nitsche's penalty parameter
Parameters(1).Density=1;                         % Density
Parameters(1).YoungsModulus=@(x,y,z) 250;        % Young's modulus
Parameters(1).PoissonsRatio=@(x,y,z) 0.3;        % Poisson's ratio
E=Parameters(1).YoungsModulus; nu=Parameters(1).PoissonsRatio;
Parameters(1).Displacement=...                   % Displacement
  @(x,y,z,t) [1./E(x,y,z).*(sin(2.*pi.*y).*(-1+cos(2.*pi.*x)).*(2+2.*nu(x,y,z)))+x.*y.*sin(pi.*x).*sin(pi.*y).*((1+nu(x,y,z)).*(1-2.*nu(x,y,z)))./(1-nu(x,y,z)-2.*nu(x,y,z).^2+E(x,y,z).*nu(x,y,z)),...
              1./E(x,y,z).*(sin(2.*pi.*x).*(+1-cos(2.*pi.*y)).*(2+2.*nu(x,y,z)))+x.*y.*sin(pi.*x).*sin(pi.*y).*((1+nu(x,y,z)).*(1-2.*nu(x,y,z)))./(1-nu(x,y,z)-2.*nu(x,y,z).^2+E(x,y,z).*nu(x,y,z))];
Parameters(1).Traction=...                        % Traction
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x];
Parameters(1).Force=...                          % Force
  @(x,y,z,t) [(24.*pi.^2.*sin(2.*pi.*y)-24.*pi.^2.*nu(x,y,z).*sin(2.*pi.*y)-48.*pi.^2.*nu(x,y,z).^2.*sin(2.*pi.*y)+E(x,y,z).*sin(pi.*x).*sin(pi.*y)-64.*pi.^2.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+24.*E(x,y,z).*pi.^2.*nu(x,y,z).*sin(2.*pi.*y)+pi.*E(x,y,z).*x.*cos(pi.*x).*sin(pi.*y)+2.*pi.*E(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)+4.*pi.*E(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+pi.*E(x,y,z).*y.*cos(pi.*y).*sin(pi.*x)+64.*pi.^2.*nu(x,y,z).*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+128.*pi.^2.*nu(x,y,z).^2.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)-4.*pi.*E(x,y,z).*nu(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+E(x,y,z).*pi.^2.*x.*y.*cos(pi.*x).*cos(pi.*y)-3.*E(x,y,z).*pi.^2.*x.*y.*sin(pi.*x).*sin(pi.*y)-64.*E(x,y,z).*pi.^2.*nu(x,y,z).*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+4.*E(x,y,z).*pi.^2.*nu(x,y,z).*x.*y.*sin(pi.*x).*sin(pi.*y))./(2.*nu(x,y,z)-2.*E(x,y,z).*nu(x,y,z)+4.*nu(x,y,z).^2-2),...
              (24.*pi.^2.*nu(x,y,z).*sin(2.*pi.*x)-24.*pi.^2.*sin(2.*pi.*x)+48.*pi.^2.*nu(x,y,z).^2.*sin(2.*pi.*x)+E(x,y,z).*sin(pi.*x).*sin(pi.*y)+64.*pi.^2.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-24.*E(x,y,z).*pi.^2.*nu(x,y,z).*sin(2.*pi.*x)+pi.*E(x,y,z).*x.*cos(pi.*x).*sin(pi.*y)+4.*pi.*E(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)+2.*pi.*E(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+pi.*E(x,y,z).*y.*cos(pi.*y).*sin(pi.*x)-64.*pi.^2.*nu(x,y,z).*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-128.*pi.^2.*nu(x,y,z).^2.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+E(x,y,z).*pi.^2.*x.*y.*cos(pi.*x).*cos(pi.*y)-3.*E(x,y,z).*pi.^2.*x.*y.*sin(pi.*x).*sin(pi.*y)+64.*E(x,y,z).*pi.^2.*nu(x,y,z).*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)+4.*E(x,y,z).*pi.^2.*nu(x,y,z).*x.*y.*sin(pi.*x).*sin(pi.*y))./(2.*nu(x,y,z)-2.*E(x,y,z).*nu(x,y,z)+4.*nu(x,y,z).^2-2)];
Parameters(2).Formulation='ElasticityLinear_HDG';% Formulation
Parameters(2).Problem='Structural';              % Problem
Parameters(2).PostProcessingHDG='yes';           % Perform HDG postprocessing
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Assumption='PlaneStrain';          % Assumption
Parameters(2).Degree=1;                          % Degree
Parameters(2).StabDisplacement=2.5e2;            % Stabilization for displacement
Parameters(2).Density=1;                         % Density
Parameters(2).YoungsModulus=@(x,y,z) 25;         % Young's modulus
Parameters(2).PoissonsRatio=@(x,y,z) 0.49999;    % Poisson's ratio
E=Parameters(2).YoungsModulus; nu=Parameters(2).PoissonsRatio;
Parameters(2).Displacement=...                   % Displacement
  @(x,y,z,t) [1./E(x,y,z).*(sin(2.*pi.*y).*(-1+cos(2.*pi.*x)).*(2+2.*nu(x,y,z)))+x.*y.*sin(pi.*x).*sin(pi.*y).*((1+nu(x,y,z)).*(1-2.*nu(x,y,z)))./(1-nu(x,y,z)-2.*nu(x,y,z).^2+E(x,y,z).*nu(x,y,z)),...
              1./E(x,y,z).*(sin(2.*pi.*x).*(+1-cos(2.*pi.*y)).*(2+2.*nu(x,y,z)))+x.*y.*sin(pi.*x).*sin(pi.*y).*((1+nu(x,y,z)).*(1-2.*nu(x,y,z)))./(1-nu(x,y,z)-2.*nu(x,y,z).^2+E(x,y,z).*nu(x,y,z))];
Parameters(2).Traction=...                       % Traction
  @(x,y,z,t,b,nx,ny,nz) [0*x,0*x];
Parameters(2).Force=...                          % Force
  @(x,y,z,t) [(24.*pi.^2.*sin(2.*pi.*y)-24.*pi.^2.*nu(x,y,z).*sin(2.*pi.*y)-48.*pi.^2.*nu(x,y,z).^2.*sin(2.*pi.*y)+E(x,y,z).*sin(pi.*x).*sin(pi.*y)-64.*pi.^2.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+24.*E(x,y,z).*pi.^2.*nu(x,y,z).*sin(2.*pi.*y)+pi.*E(x,y,z).*x.*cos(pi.*x).*sin(pi.*y)+2.*pi.*E(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)+4.*pi.*E(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+pi.*E(x,y,z).*y.*cos(pi.*y).*sin(pi.*x)+64.*pi.^2.*nu(x,y,z).*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+128.*pi.^2.*nu(x,y,z).^2.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)-4.*pi.*E(x,y,z).*nu(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+E(x,y,z).*pi.^2.*x.*y.*cos(pi.*x).*cos(pi.*y)-3.*E(x,y,z).*pi.^2.*x.*y.*sin(pi.*x).*sin(pi.*y)-64.*E(x,y,z).*pi.^2.*nu(x,y,z).*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y)+4.*E(x,y,z).*pi.^2.*nu(x,y,z).*x.*y.*sin(pi.*x).*sin(pi.*y))./(2.*nu(x,y,z)-2.*E(x,y,z).*nu(x,y,z)+4.*nu(x,y,z).^2-2),...
              (24.*pi.^2.*nu(x,y,z).*sin(2.*pi.*x)-24.*pi.^2.*sin(2.*pi.*x)+48.*pi.^2.*nu(x,y,z).^2.*sin(2.*pi.*x)+E(x,y,z).*sin(pi.*x).*sin(pi.*y)+64.*pi.^2.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-24.*E(x,y,z).*pi.^2.*nu(x,y,z).*sin(2.*pi.*x)+pi.*E(x,y,z).*x.*cos(pi.*x).*sin(pi.*y)+4.*pi.*E(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)+2.*pi.*E(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+pi.*E(x,y,z).*y.*cos(pi.*y).*sin(pi.*x)-64.*pi.^2.*nu(x,y,z).*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-128.*pi.^2.*nu(x,y,z).^2.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*x.*cos(pi.*y).*sin(pi.*x)-4.*pi.*E(x,y,z).*nu(x,y,z).*y.*cos(pi.*x).*sin(pi.*y)+E(x,y,z).*pi.^2.*x.*y.*cos(pi.*x).*cos(pi.*y)-3.*E(x,y,z).*pi.^2.*x.*y.*sin(pi.*x).*sin(pi.*y)+64.*E(x,y,z).*pi.^2.*nu(x,y,z).*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x)+4.*E(x,y,z).*pi.^2.*nu(x,y,z).*x.*y.*sin(pi.*x).*sin(pi.*y))./(2.*nu(x,y,z)-2.*E(x,y,z).*nu(x,y,z)+4.*nu(x,y,z).^2-2)];
clear E nu
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_square_coupled_4parts_2'};      % Mesh file
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
Boundaries(1).Dirichlet=[1,3,6,8];               % Dirichlet portion
Boundaries(1).Interface=[2,4,5,7];               % Interface portion
Boundaries(1).Neumann=[];                        % Neumann portion
Boundaries(2).Dirichlet=[1,4,6,7];               % Dirichlet portion
Boundaries(2).Interface=[2,3,5,8];               % Interface portion
Boundaries(2).Neumann=[];                        % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'Displacement';
   'DisplacementPost'};
Options.Test=...                                 % Test
  ['abs(Results(1).DisplacementErrorL2    -6.272951140645987e-03)<1e-12 && ',...
   'abs(Results(2).DisplacementErrorL2    -1.268994341859454e-01)<1e-12 && ',...
   'abs(Results(2).DisplacementPostErrorL2-1.365773359217927e-01)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main