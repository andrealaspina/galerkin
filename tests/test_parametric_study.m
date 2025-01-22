clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='ParametricStudy';               % Simulation type
Simulation.Problem='Thermal';                    % Problem
Simulation.ParameterStudy='StabTemperature';     % Parameter for parametric study
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Thermal_HDG';            % Formulation
Parameters.Problem='Thermal';                    % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=1;                             % Degree
Parameters.StabTemperature=10.^[0,1,2];          % Stabilization for temperature
Parameters.Density=0;                            % Density
Parameters.SpecificHeatCapacity=0;               % Specific heat capacity
Parameters.ThermalConductivity=1;                % Thermal conductivity
Parameters.ScaledTemperatureGradient=...         % Scaled temperature gradient
  @(x,y,z,t) [pi/2*x./sqrt(x.^2+y.^2).*sin(pi/2*sqrt(y.^2+x.^2)),...
              pi/2*y./sqrt(x.^2+y.^2).*sin(pi/2*sqrt(y.^2+x.^2))];
Parameters.Temperature=...                       % Temperature
  @(x,y,z,t) cos(pi/2*sqrt(x.^2+y.^2));
Parameters.ThermalFlux=@(x,y,t) 0*x;             % Thermal flux
Parameters.HeatSource=...                        % Heat source
  @(x,y,z,t) pi/2*(1./sqrt(x.^2+y.^2).*sin(pi/2*sqrt(x.^2+y.^2))+...
                                  pi/2*cos(pi/2*sqrt(x.^2+y.^2)));
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square_1'};                      % Mesh file
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
Boundaries.Dirichlet=[1,2,3,4];                  % Dirichlet portion
Boundaries.Neumann=[];                           % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeQuantityEnd=...                   % Compute quantity of interest (end of simulation)
  ['[~,NodesX]=find(Mesh.Nodes(1,:)==0);',...
   '[~,NodesY]=find(Mesh.Nodes(2,:)==0);',...
   'Nodes=intersect(NodesX,NodesY);',...
   'Results.CenterTemperature(iS2)=mean(Results.Temperature(Nodes,:));',...
   'Results.ConditionNumber(iS2)=cond(full(System.Lhs));'];
Options.ComputeError=...                         % Compute error
  {'ScaledTemperatureGradient';
   'Temperature';
   'TemperaturePost'};
Options.Test=...                                 % Test
  ['abs(Results.CenterTemperature(end)-1.284292372076179e+00)<1e-12 && ',...
   'abs(Results.ConditionNumber(end)  -1.579099980331070e+02)<1e-11'];
% --------------------------------------------------------------------------------------------------

%% Main

main