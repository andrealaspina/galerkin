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
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Electromagnetic2FluidCoupling='yes';  % Electromagnetic->Fluid coupling
Parameters.Fluid2ElectromagneticCoupling='yes';  % Fluid->Electromagnetic coupling
Parameters.Degree=6;                             % Degree
Parameters.StabElectricField=10;                 % Stabilization for electric field
Parameters.SpecificHeatRatio=7/5;                % Ratio of specific heats
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=1;               % Electric conductivity
Parameters.ParticleCharge=1;                     % Particle charge
Parameters.ParticleMass=1;                       % Particle mass
Parameters.DampingFunction=@(x,y,z) [0*x,0*x];   % Damping function for PML
Parameters.Density=...                           % Density
  @(x,y,z,t) 1+(x-1).*(x+1).*(y-1).*(y+1)*t;
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [+(1+(x-1).*(x+1).*(y-1).*(y+1)*t).*y*t,...
              -(1+(x-1).*(x+1).*(y-1).*(y+1)*t).*x*t];
Parameters.Energy=...                            % Energy
  @(x,y,z,t) 1+(x-1).*(x+1).*(y-1).*(y+1)*t;
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,t) (x-y).*(x.*y-t^2);
Parameters.ElectricField=...                     % Electric field
  @(x,y,z,t) [x.*(x-2*y)*t,...
              y.*(y-2*x)*t];
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.ForceDensity=...                      % Force in density equation
  @(x,y,z,t) -2.*t.^2.*x.^3.*y+2.*t.^2.*x.*y.^3+x.^2.*y.^2-x.^2-y.^2+1;
Parameters.ForceMomentum=...                     % Force in momentum equation
  @(x,y,z,t) [y.*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)-(7.*t.^2.*x.*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1))./5+(2.*t.*(x-1).*(y-1).*(y+1))./5+(2.*t.*(x+1).*(y-1).*(y+1))./5-t.*x.*(x-2.*y).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)-2.*t.^3.*x.*y.^2.*(x.^2-1)+2.*t.^3.*x.*y.^2.*(y.^2-1)-(2.*t.^3.*x.*(x.^2+y.^2).*(y.^2-1))./5+t.*x.*(x-y).*(-t.^2+x.*y).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)+t.*y.*(x-1).*(x+1).*(y-1).*(y+1),...
              t.*y.*(2.*x-y).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)-(7.*t.^2.*y.*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1))./5-x.*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)+(2.*t.*(x-1).*(x+1).*(y-1))./5+(2.*t.*(x-1).*(x+1).*(y+1))./5+2.*t.^3.*x.^2.*y.*(x.^2-1)-2.*t.^3.*x.^2.*y.*(y.^2-1)-(2.*t.^3.*y.*(x.^2+y.^2).*(x.^2-1))./5+t.*y.*(x-y).*(-t.^2+x.*y).*(t.*x.^2.*y.^2-t.*x.^2-t.*y.^2+t+1)-t.*x.*(x-1).*(x+1).*(y-1).*(y+1)];
Parameters.ForceEnergy=...                       % Force in energy equation
  @(x,y,z,t) (2.*t.^4.*x.^5.*y)./5-(2.*t.^4.*x.*y.^5)./5-3.*t.^3.*x.^4.*y.^3+3.*t.^3.*x.^4.*y+3.*t.^3.*x.^3.*y.^4-3.*t.^3.*x.^3.*y.^2+3.*t.^3.*x.^2.*y.^3-3.*t.^3.*x.^2.*y-3.*t.^3.*x.*y.^4+3.*t.^3.*x.*y.^2-(14.*t.^2.*x.^3.*y)./5-3.*t.^2.*x.^2.*y+(14.*t.^2.*x.*y.^3)./5+3.*t.^2.*x.*y.^2+x.^2.*y.^2-x.^2-y.^2+1;
Parameters.CurrentDensity=...                    % Current density
  @(x,y,z,t) [t.*(t-y-t.*y+2.*x.*y+t.*y.^3-x.^2+t.*x.^2.*y-t.*x.^2.*y.^3),...
              t.*(t+x+t.*x+2.*x.*y-t.*x.^3-y.^2-t.*x.*y.^2+t.*x.^3.*y.^2)];
Parameters.IncidentField=@(x,y,z,t) [0*x,0*x];   % Incident field for absorbing boundary conditions
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-10;                          % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/2;                           % Time step size
Time.BDFOrder=3;                                 % BDF order
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
  ['Results.DensityErrorL2             <1e-12 && ',...
   'Results.MomentumErrorL2            <1e-12 && ',...
   'Results.EnergyErrorL2              <1e-12 && ',...
   'Results.MagneticFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldErrorL2       <1e-12 && ',...
   'Results.MagneticFieldAuxErrorL2    <1e-12 && ',...
   'Results.ElectricFieldAuxErrorL2    <1e-12 && ',...
   'Results.ElectricFieldErrorHcurl    <1e-11 && ',...
   'Results.ElectricFieldPostErrorHcurl<1e-10'];
% --------------------------------------------------------------------------------------------------

%% Main

main