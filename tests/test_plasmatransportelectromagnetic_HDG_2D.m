clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='PlasmaTransportElectromagnetic';% Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='PlasmaTransportElectromagnetic_HDG';% Formulation
Parameters.Problem='PlasmaTransportElectromagnetic';% Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=1;                             % Degree
Parameters.StabDensity=10;                       % Stabilization for density
Parameters.StabElectricField=10;                 % Stabilization for electric field
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.ParticleCharge=1;                     % Particle charge
Parameters.ParticleMass=1;                       % Particle mass
Parameters.DomainCoordinates=[-1/2,1/2;-1/2,1/2];% Domain coordinates
Parameters.PMLWidth=[1/6,1/6];                   % Width of PML region in x,y,z
Parameters.PolynomialScaling=2;                  % Polynomial scaling for PML
Parameters.DampingConstants=[100,100];           % Damping constants for PML

% Derived parameters
xi=Parameters.DomainCoordinates(1,1);
xf=Parameters.DomainCoordinates(1,2);
yi=Parameters.DomainCoordinates(2,1);
yf=Parameters.DomainCoordinates(2,2);
dx=Parameters.PMLWidth(1);
dy=Parameters.PMLWidth(2);
m=Parameters.PolynomialScaling;
sigmax=Parameters.DampingConstants(1);
sigmay=Parameters.DampingConstants(2);

Parameters.DampingFunction=...                   % Damping function for PML
  @(x,y,z) [abs(x-xi).^m/dx^m*sigmax.*(x<xi)+abs(x-xf).^m/dx^m*sigmax.*(x>xf),...
            abs(y-yi).^m/dy^m*sigmay.*(y<yi)+abs(y-yf).^m/dy^m*sigmay.*(y>yf)];
Parameters.Density=...                           % Density
  @(x,y,z,t) 1/4*(1+cos(pi*(x-1/6)/0.2)).*(1+cos(pi*(y-1/6)/0.2)).*...
             ((((x-1/6)/0.2).^2+((y-1/6)/0.2).^2)<=1);
Parameters.MagneticField=@(x,y,z,t) 0*x;         % Magnetic field
Parameters.ElectricField=@(x,y,z,t) [0*x,0*x];   % Electric field
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.Velocity=@(x,y,z,t) [-y,+x];          % Velocity
Parameters.IncidentField=@(x,y,z,t) [0*x,0*x];   % Incident field for absorbing boundary conditions
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
load('Mesh_square01_3.mat');                     % Mesh file
Mesh.Nodes(1,:)=(Mesh.Nodes(1,:)-1/2)*(1+2*dx);
Mesh.Nodes(2,:)=(Mesh.Nodes(2,:)-1/2)*(1+2*dy);
clear xi xf yi yf dx dy m sigmax sigmay
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=2*pi;                             % Final time
Time.TimeStepSize=2*pi/4;                        % Time step size
Time.BDFOrder=2;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_r=[1,2,3,4];                % Dirichlet portion for density
Boundaries.Dirichlet_e=[];                       % Dirichlet portion for electric field
Boundaries.Absorbing=[1,2,3,4];                  % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'Density'          ,'L2';
   'MagneticField'    ,'L2';
   'ElectricField'    ,'L2';
   'MagneticFieldAux' ,'L2';
   'ElectricFieldAux' ,'L2';
   'ElectricField'    ,'Hcurl';
   'ElectricFieldPost','Hcurl'};
Options.Test=...                                 % Test
  ['abs(Results.DensityErrorL2             -1.435747461239719e-01)<1e-12 && ',...
   'abs(Results.MagneticFieldErrorL2       -2.321210511829748e-03)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorL2       -2.137550026583197e-02)<1e-12 && ',...
   'abs(Results.MagneticFieldAuxErrorL2    -2.386601820555109e-02)<1e-12 && ',...
   'abs(Results.ElectricFieldAuxErrorL2    -3.770096649303220e-01)<1e-11 && ',...
   'abs(Results.ElectricFieldErrorHcurl    -3.306307660434313e-02)<1e-12 && ',...
   'abs(Results.ElectricFieldPostErrorHcurl-3.302490555045126e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main