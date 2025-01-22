clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Plasma1FluidElectromagnetic';% Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Plasma1FluidElectromagnetic_HDG';% Formulation
Parameters.Problem='Plasma1FluidElectromagnetic';% Problem
Parameters.Axisymmetric='yes';                   % Axisymmetric
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Electromagnetic2FluidCoupling='yes';  % Electromagnetic->Fluid coupling
Parameters.Fluid2ElectromagneticCoupling='yes';  % Fluid->Electromagnetic coupling
Parameters.Degree=2;                             % Degree
Parameters.StabElectricField=10;                 % Stabilization for electric field
Parameters.SpecificHeatRatio=7/5;                % Ratio of specific heats
Parameters.InnerRadius=3/4;                      % Inner radius
Parameters.WaveNumber=2*pi;                      % Wave number
Parameters.MagneticPermeability=@(x,y,z)1*(x==x);% Magnetic permeability
Parameters.ElectricPermittivity=...              % Electric permittivity
  @(x,y,z) 1.5^2*(abs(y)<=Parameters.InnerRadius)+...
           1.0^2*(abs(y)> Parameters.InnerRadius);
Parameters.ElectricConductivity=@(x,y,z) 0*x;    % Electric conductivity
Parameters.ParticleCharge=0;                     % Particle charge
Parameters.ParticleMass=1;                       % Particle mass
Parameters.DampingFunction=@(x,y,z) [0*x,0*x];   % Damping function for PML
Parameters.Eigenmode=solveEigenmode(Parameters); % Eigenmode
Parameters.Density=@(x,y,z,t) 1*(x==x);          % Density
Parameters.Momentum=@(x,y,z,t) [0*x,0*x];        % Momentum
Parameters.Energy=@(x,y,z,t) 1*(x==x);           % Energy
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,t)  computeSolution(x,y,t,Parameters,'Hz');
Parameters.ElectricField=...                     % Electric field
  @(x,y,z,t) [computeSolution(x,y,t,Parameters,'Ex'),...
              computeSolution(x,y,t,Parameters,'Ey')];
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.ForceDensity=@(x,y,z,t) 0*x;          % Force in density equation
Parameters.ForceMomentum=@(x,y,z,t) [0*x,0*x];   % Force in momentum equation
Parameters.ForceEnergy=@(x,y,z,t) 0*x;           % Force in energy equation
Parameters.CurrentDensity=@(x,y,z,t) [0*x,0*x];  % Current density
Parameters.IncidentField=@(x,y,z,t) [0*x,0*x];   % Incident field for absorbing boundary conditions
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_axisymmetric_maxwell_test'};     % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/2^3;                         % Time step size
Time.BDFOrder=2;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_r=[1,2,3,4];                % Dirichlet portion for density
Boundaries.Dirichlet_w=[1,2,3,4];                % Dirichlet portion for momentum
Boundaries.Dirichlet_g=[1,2,3,4];                % Dirichlet portion for energy
Boundaries.Dirichlet_e=[1,2,3,4];                % Dirichlet portion for electric field
Boundaries.FarField=[];                          % Far-field portion
Boundaries.InviscidWall=[];                      % Inviscid wall portion
Boundaries.Absorbing=[];                         % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'Density'          ,'L2';
   'Momentum'         ,'L2';
   'Energy'           ,'L2';
   'MagneticField'    ,'L2';
   'ElectricField'    ,'L2';
   'MagneticFieldAux' ,'L2';
   'ElectricFieldAux' ,'L2';
   'ElectricField'    ,'Hcurl';
   'ElectricFieldPost','Hcurl'};
Options.Test=...                                 % Test
  ['abs(Results.DensityErrorL2             -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.MomentumErrorL2            -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.EnergyErrorL2              -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.MagneticFieldErrorL2       -3.435810475774289e-01)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorL2       -1.273922024834438e-01)<1e-12 && ',...
   'abs(Results.MagneticFieldAuxErrorL2    -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldAuxErrorL2    -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorHcurl    -1.176233027969729e+00)<1e-12 && ',...
   'abs(Results.ElectricFieldPostErrorHcurl-1.468440342698116e+00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Eigenmode

function [beta]=solveEigenmode(Parameters)

% Parameters
a=Parameters.InnerRadius;
k=Parameters.WaveNumber;
eps1=Parameters.ElectricPermittivity(0,a-1e-3);
eps2=Parameters.ElectricPermittivity(0,a+1e-3);

% Variables
u=@(beta) a*sqrt(+(eps1*k^2-beta^2));
w=@(beta) a*sqrt(-(eps2*k^2-beta^2));

% Function
fun=@(beta) eps1*besselj(1,u(beta))/(u(beta)*besselj(0,u(beta)))+...
            eps2*besselk(1,w(beta))/(w(beta)*besselk(0,w(beta)));

% Initial guess
beta0=sqrt(2)/2*sqrt(eps1+eps2)*k;

% Solve
beta=fzero(fun,beta0);

end

%% Solution

function [sol]=computeSolution(x,y,t,Parameters,var)

% Parameters
a=Parameters.InnerRadius;
k=Parameters.WaveNumber;
eps1=Parameters.ElectricPermittivity(0,a-1e-3);
eps2=Parameters.ElectricPermittivity(0,a+1e-3);
beta=Parameters.Eigenmode;
omega=k;

% Variables
u=a*sqrt(+(eps1*k^2-beta^2));
w=a*sqrt(-(eps2*k^2-beta^2));

% Radial coordinate
r=abs(y);

% Return solution
switch var
  case 'Ex'
    Ex=                                      besselj(0,u/a*r).*cos(omega*t-beta*x).*(r<=a)...
                  +besselj(0,u)/besselk(0,w)*besselk(0,w/a*r).*cos(omega*t-beta*x).*(r> a);
    sol=Ex;
  case 'Ey'
    Er=  -beta*a/u                          *besselj(1,u/a*r).*sin(omega*t-beta*x).*(r<=a)...
         +beta*a/w*besselj(0,u)/besselk(0,w)*besselk(1,w/a*r).*sin(omega*t-beta*x).*(r> a);
    sol=sign(y).*Er;
  case 'Hz'
    Hp=-k*eps1*a/u                          *besselj(1,u/a*r).*sin(omega*t-beta*x).*(r<=a)...
       +k*eps2*a/w*besselj(0,u)/besselk(0,w)*besselk(1,w/a*r).*sin(omega*t-beta*x).*(r> a);
    sol=sign(y).*Hp;
end

end