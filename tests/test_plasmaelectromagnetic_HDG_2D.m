clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='PlasmaElectromagnetic';      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='PlasmaElectromagnetic_HDG';% Formulation
Parameters.Problem='PlasmaElectromagnetic';      % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=1;                             % Degree
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
Parameters.MagneticField=@(x,y,z,t) 0*x;         % Magnetic field
Parameters.ElectricField=@(x,y,z,t) [0*x,0*x];   % Electric field
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.Density=...                           % Density
  @(x,y,z,t) 1/4*(1+cos(pi*(x-sqrt(2)/6*cos(pi/4+t))/0.2)).*...
                 (1+cos(pi*(y-sqrt(2)/6*sin(pi/4+t))/0.2)).*...
             ((((x-sqrt(2)/6*cos(pi/4+t))/0.2).^2+...
               ((y-sqrt(2)/6*sin(pi/4+t))/0.2).^2)<=1);
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
Boundaries.Dirichlet_e=[];                       % Dirichlet portion for electric field
Boundaries.Absorbing=[1,2,3,4];                  % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'MagneticField'    ,'L2';
   'ElectricField'    ,'L2';
   'MagneticFieldAux' ,'L2';
   'ElectricFieldAux' ,'L2';
   'ElectricField'    ,'Hcurl';
   'ElectricFieldPost','Hcurl'};
Options.Test=...                                 % Test
  ['abs(Results.MagneticFieldErrorL2       -3.310721370997904e-03)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorL2       -6.368599921664606e-02)<1e-12 && ',...
   'abs(Results.MagneticFieldAuxErrorL2    -1.478371922608692e-02)<1e-12 && ',...
   'abs(Results.ElectricFieldAuxErrorL2    -1.885105947222024e-01)<1e-11 && ',...
   'abs(Results.ElectricFieldErrorHcurl    -6.619126535608760e-02)<1e-12 && ',...
   'abs(Results.ElectricFieldPostErrorHcurl-6.620275665127154e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main