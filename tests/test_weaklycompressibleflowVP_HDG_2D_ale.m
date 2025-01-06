clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='FluidALE';                   % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters(1).Formulation='WeaklyCompressibleFlowVP_HDG';% Formulation
Parameters(1).Problem='Fluid';                   % Problem
Parameters(1).PostProcessingHDG='no';            % Perform HDG postprocessing
Parameters(1).ConvectiveFlow='yes';              % Convective flow
Parameters(1).ArbitraryLagrangianEulerian='yes'; % Arbitrary Lagrangian-Eulerian description
Parameters(1).Degree=1;                          % Degree
Parameters(1).StabVelocity=1;                    % Stabilization for velocity
Parameters(1).StabPressure=10;                   % Stabilization for pressure
Parameters(1).ReferenceDensity=1;                % Reference density
Parameters(1).ReferencePressure=0;               % Reference pressure
Parameters(1).CompressibilityCoeff=0.1;          % Compressibility coefficient
Parameters(1).DynamicViscosity=0.1;              % Dynamic viscosity
r0=Parameters(1).ReferenceDensity;      p0=Parameters(1).ReferencePressure;
eps=Parameters(1).CompressibilityCoeff; mu=Parameters(1).DynamicViscosity;
vx=@(x,y,t) sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0));
vy=@(x,y,t) cos(pi.*x).*cos(pi.*y).*sin(pi.*t)-eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0));
p= @(x,y,t) pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y);
Parameters(1).ScaledStrainRate=...               % Scaled strain rate
  @(x,y,z,t) [(mu.^(1./2).*pi.*(6.^(1./2).*eps.*pi.*cos(pi.*x).^2.*sin(pi.*t).^2-6.^(1./2).*eps.*pi.*sin(pi.*t).^2+6.^(1./2).*eps.*pi.*cos(pi.*y).^2.*sin(pi.*t).^2-6.*2.^(1./2).*r0.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)+6.^(1./2).*eps.*pi.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)))./(6.*r0),...
              (mu.^(1./2).*pi.*(6.^(1./2).*eps.*pi.*cos(pi.*x).^2.*sin(pi.*t).^2-6.^(1./2).*eps.*pi.*sin(pi.*t).^2+6.^(1./2).*eps.*pi.*cos(pi.*y).^2.*sin(pi.*t).^2+6.*2.^(1./2).*r0.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)+6.^(1./2).*eps.*pi.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)))./(6.*r0),...
              -(eps.*mu.^(1./2).*pi.^2.*(x.*pi.*sin(pi.*t).^2.*sin(2.*pi.*y)-2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x)+y.*pi.*sin(pi.*t).^2.*sin(2.*pi.*x)))./(2.*r0)];
Parameters(1).Velocity=...                       % Velocity
  @(x,y,z,t) [vx(x,y,t),vy(x,y,t)];
Parameters(1).Pressure=@(x,y,z,t) p(x,y,t);      % Pressure
Parameters(1).Traction=...                       % Traction
  @(x,y,z,t) [0*x, 0*x];
Parameters(1).ResidualContinuity=...             % Residual of continuity equation
  @(x,y,z,t) pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y)-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t));
Parameters(1).Force=...                          % Force
  @(x,y,z,t) [(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))-mu.*(eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0)+(pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))+eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^4.*x.*cos(2.*pi.*y).*sin(pi.*t).^2)./r0)+(2.*eps.*((pi.^3.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0)+(pi.^3.*cos(pi.*x).*sin(pi.*t).^2.*sin(pi.*x))./r0))./3-2.*pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-(eps.*((pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(4.*r0)-(pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*x.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0))-pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y)-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).^2+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)),...
              mu.*(eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^4.*y.*cos(2.*pi.*x).*sin(pi.*t).^2)./r0)-eps.*((pi.^3.*sin(pi.*t).^2.*sin(2.*pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))+(2.*eps.*((pi.^3.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0)-(pi.^3.*cos(pi.*y).*sin(pi.*t).^2.*sin(pi.*y))./r0))./3+2.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-(eps.*((pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t))./(2.*r0)+(pi.^2.*cos(pi.*t).*sin(pi.*t).*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(4.*r0))-pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.^2.*cos(pi.*t).*cos(pi.*y).*sin(pi.*x))./(2.*r0)-(pi.^3.*y.*sin(pi.*t).^2.*sin(2.*pi.*x))./(2.*r0))+pi.*cos(pi.*y).*sin(pi.*t).*sin(pi.*x)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)))+2.*(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))+pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+(eps.*((pi.*sin(pi.*t).^2.*(2.*pi.*cos(2.*pi.*x)+2.*pi.*cos(2.*pi.*y)))./(8.*r0)+(pi.^2.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y))./(2.*r0))-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y)).*(r0-eps.*(p0-pi.*cos(pi.*x).*sin(pi.*t).*sin(pi.*y))).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))+pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t)+pi.^2.*eps.*cos(pi.*x).*cos(pi.*y).*sin(pi.*t).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t)).^2-pi.^2.*eps.*cos(pi.*t).*cos(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))-pi.^2.*eps.*sin(pi.*t).*sin(pi.*x).*sin(pi.*y).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*x)+2.*pi.*x.*cos(2.*pi.*y)))./(8.*r0)+(pi.*cos(pi.*t).*sin(pi.*x).*sin(pi.*y))./(2.*r0))-sin(pi.*t).*sin(pi.*x).*sin(pi.*y)).*(eps.*((pi.*sin(pi.*t).^2.*(sin(2.*pi.*y)+2.*pi.*y.*cos(2.*pi.*x)))./(8.*r0)-(pi.*cos(pi.*t).*cos(pi.*x).*cos(pi.*y))./(2.*r0))-cos(pi.*x).*cos(pi.*y).*sin(pi.*t))];
Parameters(1).Density=@(Pressure)...             % Equation of state
  Parameters(1).ReferenceDensity+...
 (Pressure-Parameters(1).ReferencePressure)*Parameters(1).CompressibilityCoeff;
Parameters(1).DDensityDPressure=@(Pressure)...   % dDensity/dPressure
  Parameters(1).CompressibilityCoeff;
Parameters(2).Formulation='Elasticity_CG';       % Formulation
Parameters(2).Problem='Mesh';                    % Problem
Parameters(2).Model='LinearElasticity';          % Model
Parameters(2).Assumption='PlaneStrain';          % Assumption
Parameters(2).Degree=1;                          % Degree
Parameters(2).NitschePenalty=1000;               % Nitsche's penalty parameter
Parameters(2).Density=1;                         % Density
Parameters(2).YoungsModulus=@(x,y,z) 5/4;        % Young's modulus
Parameters(2).PoissonsRatio=@(x,y,z) 1/4;        % Poisson's ratio
Parameters(2).Displacement=...                   % Displacement
  @(x,y,z,t) [1/4*sin(2*pi*x).*(1-cos(2*pi*y)).*(1-cos(2*pi*t))*1/8,...
              1/4*sin(2*pi*y).*(1-cos(2*pi*x)).*(1-cos(2*pi*t))*1/8];
Parameters(2).Traction=@(x,y,z,t) [0*x,0*x];     % Traction
Parameters(2).Force=...                          % Force
  @(x,y,z,t) [-pi^2/16*sin(2*pi*x).*(6*cos(2*pi*y)+cos(2*pi*t)-4*cos(2*pi*y)*cos(2*pi*t)-3),...
              -pi^2/16*sin(2*pi*y).*(6*cos(2*pi*x)+cos(2*pi*t)-4*cos(2*pi*x)*cos(2*pi*t)-3)];
clear r0 p0 eps mu vx vy p
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_square01_ale_2'};                % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-3;                           % Tolerance
System.MaxIterations=100;                        % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=0.5;                              % Final time
Time.TimeStepSize=2^(-3);                        % Time step size
Time.BDFOrder=4;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries(1).Dirichlet_v_x=[1,2,3,4];           % Dirichlet portion for x-velocity
Boundaries(1).Dirichlet_v_y=[1,2,3,4];           % Dirichlet portion for y-velocity
Boundaries(1).Dirichlet_p=[1,2,3,4];             % Dirichlet portion for pressure
Boundaries(1).Neumann_t_x=[];                    % Neumann portion for x-traction
Boundaries(1).Neumann_t_y=[];                    % Neumann portion for y-traction
Boundaries(2).Dirichlet=[1,2,3,4];               % Dirichlet portion
Boundaries(2).Neumann=[];                        % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeQuantityIteration=...             % Compute quantity of interest at each iteration
  ['[~,Nodes12]=ismember(Mesh(1).NodesInitial'',Mesh(2).NodesInitial'',''Rows'');',...
   'Mesh(1).Nodes=Mesh(1).NodesInitial+[Block(2,2).SolutionGlobal(Nodes12,:)'';',...
   '                                    zeros(1,size(Mesh(1).NodesInitial,2))];',...
   'Elements(1).Nodes=mat2cell(Mesh(1).Nodes(:,Mesh(1).Elements(:)),',...
   '3,ones(Sizes(1).NumElements,1)*Sizes(1).NumElementNodes)'';'];
Options.ComputeError=...                         % Compute error
  {'ScaledStrainRate';
   'Velocity';
   'Pressure';
   'Displacement'};
Options.Test=...                                 % Test
  ['abs(Results(1).ScaledStrainRateErrorL2-4.247474847925973e-01)<1e-12 && ',...
   'abs(Results(1).VelocityErrorL2        -1.395975788997909e-01)<1e-12 && ',...
   'abs(Results(1).PressureErrorL2        -2.066939939590204e-01)<1e-12 && ',...
   'abs(Results(2).DisplacementErrorL2    -1.965976842183992e-02)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main