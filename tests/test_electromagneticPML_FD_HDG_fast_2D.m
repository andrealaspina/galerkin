clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Electromagnetic';            % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='ElectromagneticPML_FD_HDG_fast';% Formulation
Parameters.Problem='Electromagnetic';            % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=4;                             % Degree
Parameters.StabElectricField=10;                 % Stabilization for electric field
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.Frequency=2;                          % Frequency
Parameters.DomainCoordinates=[0,1;0,1];          % Domain coordinates
Parameters.PMLWidth=[1/3,0];                     % Width of PML region in x,y,z
Parameters.PolynomialScaling=1;                  % Polynomial scaling for PML
Parameters.DampingConstants=[100,0];             % Damping constants for PML

% Derived parameters
xi=Parameters.DomainCoordinates(1,1);
xf=Parameters.DomainCoordinates(1,2);
yi=Parameters.DomainCoordinates(2,1);
yf=Parameters.DomainCoordinates(2,2);
dx=Parameters.PMLWidth(1);
dy=Parameters.PMLWidth(2);
m=Parameters.PolynomialScaling;
sigmax=Parameters.DampingConstants(1);

% Additional parameters
mu0=Parameters.MagneticPermeability;
epsilon0=Parameters.ElectricPermittivity;
f=Parameters.Frequency;
omega=2*pi*f;
Hz_inc=@(x,y,w) sqrt(epsilon0/mu0)*cos(sqrt(epsilon0*mu0)*omega*x)...
            -1i*sqrt(epsilon0/mu0)*sin(sqrt(epsilon0*mu0)*omega*x);
Ex_inc=@(x,y,w) 0*x;
Ey_inc=@(x,y,w) cos(sqrt(epsilon0*mu0)*omega*x)...
            -1i*sin(sqrt(epsilon0*mu0)*omega*x);
tol=1e-9;
nx=@(x) (abs(mean(x)-xf)<tol)-(abs(mean(x)-xi)<tol);
ny=@(y) (abs(mean(y)-yf)<tol)-(abs(mean(y)-yi)<tol);

Parameters.DampingFunction=...                   % Damping function for PML
  @(x,y,z) [abs(x-xi).^m/dx^m*sigmax.*(x<xi)+abs(x-xf).^m/dx^m*sigmax.*(x>xf),...
            0*x];
Parameters.MagneticField=@(x,y,z,t) 0*x;         % Magnetic field
Parameters.ElectricField=@(x,y,z,t) [0*x,0*x];   % Electric field
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.CurrentDensity=@(x,y,z,t) [0*x,0*x];  % Current density
Parameters.IncidentField=...                     % Incident field for absorbing boundary conditions
  @(x,y,z,w) [-sqrt(epsilon0/mu0)*ny(y)*(Ex_inc(x,y,w)*ny(y)-Ey_inc(x,y,w)*nx(x))-Hz_inc(x,y,w)*ny(y),...
              +sqrt(epsilon0/mu0)*nx(x)*(Ex_inc(x,y,w)*ny(y)-Ey_inc(x,y,w)*nx(x))+Hz_inc(x,y,w)*nx(x)];
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
load('Mesh_square01_2.mat');                    % Mesh file
Mesh.Nodes(1,:)=Mesh.Nodes(1,:)*(1+dx);
Mesh.Nodes(2,:)=Mesh.Nodes(2,:)*(1+dy);
clear xi xf yi yf dx dy m sigmax mu0 epsilon0 f omega Hz_inc Ex_inc Ey_inc tol nx ny
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
System.IncrementalForm='no';                     % Incremental form
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
Time.Frequency=Parameters.Frequency;             % Frequency of time-harmonic problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,3,4];                    % Dirichlet portion
Boundaries.Absorbing=2;                          % Absorbing portion
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
  ['abs(Results.MagneticFieldErrorL2       -5.167724958057509e-01)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorL2       -5.232052938104808e-01)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorHcurl    -6.178301799660463e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldPostErrorHcurl-6.140721105406652e-00)<1e-12 && ',...
   'abs(abs(max(Results.ElectricField(Mesh.Nodes(1,:)<1,2)))-1)   <5e-02'];
% --------------------------------------------------------------------------------------------------

%% Main

main