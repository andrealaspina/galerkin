clearvars -except Files Passed CPUTime Test;

%% Input (thermal problem)

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Thermal';                    % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_CG';             % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=1e2;                   % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=11.5;             % Thermal conductivity
Parameters.ConvectionCoefficient=...             % Convection coefficient
  @(x,y,z,b) 15*((b== 2)||(b== 6)||(b== 7)||(b== 8)||(b== 9))+...
             20*((b==13))+...
             30*((b==12)||(b==14)||(b==15))+...
             40*((b== 1)||(b==10))+...
             50*((b==11))+...
           1000*((b== 3)||(b== 4)||(b== 5));
Parameters.AmbientTemperature=...                % Ambient temperature
  @(x,y,z,b)150*((b==12)||(b==14)||(b==15))+...
            300*((b== 3)||(b== 4)||(b== 5))+...
            400*((b== 2)||(b== 6)||(b== 7)||(b== 8)||(b== 9))+...
            800*((b== 1))+...
           1000*((b==10)||(b==11)||(b==13));
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) 0*x;
Parameters.TemperatureGradient=...               % Temperature gradient
  @(x,y,z,t) [0*x,0*x,0*x];
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=@(x,y,z,t) 0*x;            % Heat source               
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=importGeometry(createpde,'Blade.stl');  % Geometry definition
Mesh.MinElementSize=0.01;                        % Minimum element size
Mesh.MaxElementSize=0.01;                        % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
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
Boundaries.Dirichlet=[];                         % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
Boundaries.Robin=1:15;                           % Robin portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Temperature interpolant

TemperatureInterpolant=scatteredInterpolant(Mesh.Nodes(1,:)',Mesh.Nodes(2,:)',Mesh.Nodes(3,:)',...
  Results.Temperature,'linear','linear');

clearvars -except Files Passed CPUTime Test TemperatureInterpolant;

%% Input (structural problem)

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Elasticity_CG';          % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='LinearElasticity';             % Model
Parameters.Degree=2;                             % Degree
Parameters.NitschePenalty=1e12;                  % Nitsche's penalty parameter
Parameters.Density=0;                            % Density
Parameters.YoungsModulus=227e9;                  % Young's modulus
Parameters.PoissonsRatio=0.27;                   % Poisson's ratio
Parameters.CoefficientThermalExpansion=12.7e-6;  % Coefficient of thermal expansion
Parameters.ReferenceTemperature=300;             % Reference temperature
Parameters.Displacement=@(x,y,z,t) [0*x,0*x,0*x];% Displacement
Parameters.Traction=...                          % Traction
  @(x,y,z,t,b,nx,ny,nz) [-5.0e5*nx.*(b==11)-4.5e5*nx.*(b==10),...
                         -5.0e5*ny.*(b==11)-4.5e5*ny.*(b==10),...
                         -5.0e5*nz.*(b==11)-4.5e5*nz.*(b==10)];
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
Parameters.Force=@(x,y,z,t) [0*x,0*x,0*x];       % Force
Parameters.Temperature=...                       % Temperature
  @(x,y,z) TemperatureInterpolant(x,y,z);
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=importGeometry(createpde,'Blade.stl');  % Geometry definition
Mesh.MinElementSize=0.02;                        % Minimum element size
Mesh.MaxElementSize=0.02;                        % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
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
Boundaries.Dirichlet=3;                          % Dirichlet portion
Boundaries.Neumann=[10,11];                      % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Test=...                                 % Test
  ['abs(max(Results.Displacement(:,1))-1.192608856712264e-04)<1e-12 && ',...
   'abs(max(Results.Displacement(:,2))-1.439732238927635e-03)<1e-12 && ',...
   'abs(max(Results.Displacement(:,3))-3.936488411477511e-04)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main