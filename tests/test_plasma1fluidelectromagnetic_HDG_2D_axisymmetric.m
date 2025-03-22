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
Parameters.Degree=6;                             % Degree
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.SpecificHeatRatio=7/5;                % Ratio of specific heats
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.ParticleCharge=1;                     % Particle charge
Parameters.ParticleMass=1;                       % Particle mass
Parameters.DampingFunction=@(x,y,z) [0*x,0*x];   % Damping function for PML
Parameters.Density=...                           % Density
  @(x,y,z,t) 1+x.^2.*y.^2*t;
Parameters.Momentum=...                          % Momentum
  @(x,y,z,t) [+(1+x.^2.*y.^2*t).*x*t,...
              -(1+x.^2.*y.^2*t).*y*t];
Parameters.Energy=...                            % Energy
  @(x,y,z,t) 1+x.^2.*y.^2*t;
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,t) (2*x.^2+3*y.^2).*x.^2.*y*t^2;
Parameters.ElectricField=...                     % Electric field
  @(x,y,z,t) [+2*x.^4.*y.^2*t,...
              -2*x.^3.*y.^3*t];
Parameters.MagneticFieldAux=@(x,y,z,t) 0*x;      % Auxiliary magnetic field for PML
Parameters.ElectricFieldAux=@(x,y,z,t) [0*x,0*x];% Auxiliary electric field for PML
Parameters.ForceDensity=...                      % Force in density equation
    @(x,y,z,t) x.^2.*y.^2-t*(t*x.^2.*y.^2+1);
Parameters.ForceMomentum=...                     % Force in momentum equation
  @(x,y,z,t) [x.*(t*x.^2.*y.^2+1)+t*x.^3.*y.^2-(2*t^2*x.*(t*x.^2.*y.^2+1))/5-2*t*x.^4.*y.^2.*(t*x.^2.*y.^2+1)-(2*t*x.*y.^2.*(t^2*x.^2+t^2*y.^2-2))/5+t^3*x.^2.*y.^2.*(2*x.^2+3*y.^2).*(t*x.^2.*y.^2+1),...
              (8*t^2*y)/5-y+(4*t*x.^2.*y)/5-2*t*x.^2.*y.^3-(2*t^3*x.^4.*y)/5+(6*t^3*x.^2.*y.^3)/5+2*t*x.^3.*y.^3.*(t*x.^2.*y.^2+1)+t^3*x.^3.*y.*(2*x.^2+3*y.^2).*(t*x.^2.*y.^2+1)];
Parameters.ForceEnergy=...                       % Force in energy equation
    @(x,y,z,t) (3*t^3*y.^2)/5-(t^3*x.^2)/5-(7*t)/5+x.^2.*y.^2-(7*t^2*x.^2.*y.^2)/5+(3*t^4*x.^2.*y.^4)/5-(t^4*x.^4.*y.^2)/5-2*t^2*x.^3.*y.^4.*(t*x.^2.*y.^2+1)-2*t^2*x.^5.*y.^2.*(t*x.^2.*y.^2+1);
Parameters.CurrentDensity=...                    % Current density
  @(x,y,z,t) [4*t^2*x.^4-2*x.^4.*y.^2-t*x.*(t*x.^2.*y.^2+1)+12*t^2*x.^2.*y.^2,...
              2*x.^3.*y.^3+t*y.*(t*x.^2.*y.^2+1)-6*t^2*x.*y.^3-8*t^2*x.^3.*y];
Parameters.IncidentField=@(x,y,z,t) [0*x,0*x];   % Incident field for absorbing boundary conditions
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_axisymmetric_euler_maxwell_1'}; % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-12;                          % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1/2^1;                         % Time step size
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
   'Results.EnergyErrorL2              <1e-11 && ',...
   'Results.MagneticFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldErrorL2       <1e-11 && ',...
   'Results.MagneticFieldAuxErrorL2    <1e-12 && ',...
   'Results.ElectricFieldAuxErrorL2    <1e-12 && ',...
   'Results.ElectricFieldErrorHcurl    <1e-10 && ',...
   'Results.ElectricFieldPostErrorHcurl<1e-10'];
% --------------------------------------------------------------------------------------------------

%% Main

main