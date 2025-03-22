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
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.SpecificHeatRatio=7/5;                % Ratio of specific heats
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.ParticleCharge=0;                     % Particle charge
Parameters.ParticleMass=1;                       % Particle mass
Parameters.DampingFunction=@(x,y,z) [0*x,0*x];   % Damping function for PML
Parameters=referenceSolution(Parameters);        % Reference solution
Parameters.Density=...                           % Density
  @(x,y,z,t) (t==0)*Parameters.RefDensity(1)+...
             (t >0)*interp1(Parameters.RefDomain,Parameters.RefDensity,y,'linear','extrap');
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [0*x,...
              (t==0)*Parameters.RefMomentum(1)+...
              (t >0)*interp1(Parameters.RefDomain,Parameters.RefMomentum,y,'linear','extrap')];
Parameters.Energy=...                            % Energy
  @(x,y,z,t) (t==0)*Parameters.RefEnergy(1)+...
             (t >0)*interp1(Parameters.RefDomain,Parameters.RefEnergy,y,'linear','extrap');
Parameters.MagneticField=@(x,y,z,t) 0*x;         % Magnetic field
Parameters.ElectricField=@(x,y,z,t) [0*x,0*x];   % Electric field
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.ForceDensity=@(x,y,z,t) 0*x;          % Force in density equation
Parameters.ForceMomentum=@(x,y,z,t) [0*x,0*x];   % Force in momentum equation
Parameters.ForceEnergy=@(x,y,z,t) 0*x;           % Force in energy equation
Parameters.CurrentDensity=@(x,y,z,t) [0*x,0*x];  % Current density
Parameters.IncidentField=@(x,y,z,t) [0*x,0*x];   % Incident field for absorbing boundary conditions
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_axisymmetric_euler_1'};         % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-8;                           % Tolerance
System.MaxIterations=20;                         % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/2^3;                         % Time step size
Time.BDFOrder=1;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_r=[];                       % Dirichlet portion for density
Boundaries.Dirichlet_w=[];                       % Dirichlet portion for momentum
Boundaries.Dirichlet_g=[];                       % Dirichlet portion for energy
Boundaries.Dirichlet_e=[1,4];                    % Dirichlet portion for electric field
Boundaries.FarField=[1,4];                       % Far-field portion
Boundaries.InviscidWall=[];                      % Inviscid wall portion
Boundaries.Absorbing=[];                         % Absorbing portion
Boundaries.PeriodicMaster=2;                     % Periodic portion (master)
Boundaries.PeriodicSlave=3;                      % Periodic portion (slave)
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
  ['abs(Results.DensityErrorL2             -1.731436168759418e-04)<1e-12 && ',...
   'abs(Results.MomentumErrorL2            -9.591023302537317e-02)<1e-12 && ',...
   'abs(Results.EnergyErrorL2              -7.887462710752101e+01)<1e-10 && ',...
   'abs(Results.MagneticFieldErrorL2       -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorL2       -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.MagneticFieldAuxErrorL2    -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldAuxErrorL2    -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorHcurl    -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldPostErrorHcurl-0.000000000000000e-00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Reference solution

function [Parameters]=referenceSolution(Parameters)

% Reference domain
load('Mesh_axisymmetric_euler_1','Mesh');
rspan=minmax(Mesh.Nodes(2,:));

% Parameters
gamma=Parameters.SpecificHeatRatio;

% Inlet values
rhoI=2;
urI=200;
pI=2e5;

% Set ODE
alpha1=rspan(1)*rhoI*urI;
alpha3=rspan(1)*urI*(gamma/(gamma-1)*pI+1/2*rhoI*urI^2);
drho_dr=@(r,rho) rho/((alpha3/alpha1^3*rho^2*r^2-(gamma+1)/(2*(gamma-1)))*(gamma-1)*r);

% Solution
[r,rho]=ode23(drho_dr,rspan,rhoI,odeset('RelTol',3e-14,'AbsTol',eps));
ur=alpha1./(r.*rho);
p=(gamma-1)/gamma*(alpha3./(r.*ur)-1/2*rho.*ur.^2);

% Reference data
Parameters.RefDomain=r;
Parameters.RefDensity=rho;
Parameters.RefMomentum=rho.*ur;
Parameters.RefEnergy=p/(gamma-1)+1/2*rho.*ur.^2;

end