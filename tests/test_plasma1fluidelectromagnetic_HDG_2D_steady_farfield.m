clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Plasma1FluidElectromagnetic';% Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Plasma1FluidElectromagnetic_HDG';% Formulation
Parameters.Problem='Plasma1FluidElectromagnetic';% Problem
Parameters.Axisymmetric='no';                    % Axisymmetric
Parameters.PostProcessingHDG='no';               % Perform HDG postprocessing
Parameters.Electromagnetic2FluidCoupling='yes';  % Electromagnetic->Fluid coupling
Parameters.Fluid2ElectromagneticCoupling='yes';  % Fluid->Electromagnetic coupling
Parameters.Degree=1;                             % Degree
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.SpecificHeatRatio=7/5;                % Ratio of specific heats
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.ParticleCharge=0;                     % Particle charge
Parameters.ParticleMass=1;                       % Particle mass
Parameters.DampingFunction=@(x,y,z) [0*x,0*x];   % Damping function for PML
Parameters.Density=...                           % Density
  @(x,y,z,t) computeSolution(x,y,Parameters.SpecificHeatRatio,'r');
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [computeSolution(x,y,Parameters.SpecificHeatRatio,'wx'),...
              computeSolution(x,y,Parameters.SpecificHeatRatio,'wy')];
Parameters.Energy=...                            % Energy
  @(x,y,z,t) computeSolution(x,y,Parameters.SpecificHeatRatio,'rE');
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
Mesh.File={'Mesh_ringleb_flow_1'};               % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-10;                          % Tolerance
System.MaxIterations=2;                          % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.SteadySolution='yes';                       % Steady solution from time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=5;                                % Final time
Time.TimeStepSize=1;                             % Time step size
Time.BDFOrder=1;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet_r=[];                       % Dirichlet portion for density
Boundaries.Dirichlet_w=[];                       % Dirichlet portion for momentum
Boundaries.Dirichlet_g=[];                       % Dirichlet portion for energy
Boundaries.Dirichlet_e=[1,2,3,4];                % Dirichlet portion for electric field
Boundaries.FarField=[1,2,3,4];                   % Far-field portion
Boundaries.InviscidWall=[];                      % Inviscid wall portion
Boundaries.Absorbing=[];                         % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'Density'         ,'L2';
   'Momentum'        ,'L2';
   'Energy'          ,'L2';
   'MagneticField'   ,'L2';
   'ElectricField'   ,'L2';
   'MagneticFieldAux','L2';
   'ElectricFieldAux','L2';};
Options.Test=...                                 % Test
  ['abs(Results.DensityErrorL2         -1.822770305625944e-02)<1e-12 && ',...
   'abs(Results.MomentumErrorL2        -2.964960109875233e-02)<1e-12 && ',...
   'abs(Results.EnergyErrorL2          -3.290783500125130e-02)<1e-12 && ',...
   'abs(Results.MagneticFieldErrorL2   -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldErrorL2   -0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.MagneticFieldAuxErrorL2-0.000000000000000e-00)<1e-12 && ',...
   'abs(Results.ElectricFieldAuxErrorL2-0.000000000000000e-00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Analytical solution

function [sol]=computeSolution(x,y,gamma,var)
r=@(c) c.^(2/(gamma-1));
V=@(c) sqrt(2*(1-c.^2)/(gamma-1));
J=@(c) 1./c+1./(3*c.^3)+1./(5*c.^5)-1/2*log((1+c)./(1-c));
c=zeros(size(x));
for i=1:numel(x)
  f=@(c) (x(i)+J(c)/2)^2+y(i)^2-1/(4*r(c)^2*V(c)^4);
  c(i)=fzero(f,1/2);
end
r=r(c);
V=V(c);
theta=1/2*asin(2*r.*V.^2.*y);
switch var
  case 'c'
    sol=c;
  case 'vx'
    sol=-sign(y).*V.*sin(theta);
  case 'vy'
    sol=V.*cos(theta);
  case 'p'
    sol=1/gamma*c.^(2*gamma/(gamma-1));
  case 'M'
    sol=sqrt((-sign(y).*V.*sin(theta)).^2+(V.*cos(theta)).^2)./c;
  case 'r'
    sol=r;
  case 'wx'
    sol=r.*(-sign(y).*V.*sin(theta));
  case 'wy'
    sol=r.*(V.*cos(theta));
  case 'rE'
    sol=(1/gamma*c.^(2*gamma/(gamma-1)))/(gamma-1)...
      +1/2*r.*((-sign(y).*V.*sin(theta)).^2+(V.*cos(theta)).^2);
end
end