clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Electromagnetic';            % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Electromagnetic_HDG';    % Formulation
Parameters.Problem='Electromagnetic';            % Problem
Parameters.PostProcessingHDG='yes';              % Perform HDG postprocessing
Parameters.Degree=2;                             % Degree
Parameters.StabElectricField=1;                  % Stabilization for electric field
Parameters.MagneticPermeability=1;               % Magnetic permeability
Parameters.ElectricPermittivity=1;               % Electric permittivity
Parameters.ElectricConductivity=0;               % Electric conductivity
Parameters.MagneticField=...                     % Magnetic field
  @(x,y,z,t) [1/2*(+z/2-y/2)*t^2,...
              1/2*(-x/1-z/2)*t^2,...
              1/2*(+x/1+y/2)*t^2];
Parameters.ElectricField=...                     % Electric field
  @(x,y,z,t) [(+x.*x/2+x.*y/1+x.*z/1)*t,...
              (-x.*y/2-y.*y/4-y.*z/2)*t,...
              (-x.*z/2-y.*z/2-z.*z/4)*t];
Parameters.CurrentDensity=...                    % Current density
  @(x,y,z,t) [-x.*z/1-x.*y/1-x.^2/2+t^2/2,...
              +x.*y/2+y.*z/2+y.^2/4-t^2/4,...
              +x.*z/2+y.*z/2+z.^2/4-t^2/4];
tol=1e-9;
nx=@(x) (abs(mean(x)-(+1))<tol)-(abs(mean(x)-(0))<tol);
ny=@(y) (abs(mean(y)-(+1))<tol)-(abs(mean(y)-(0))<tol);
nz=@(z) (abs(mean(z)-(+1))<tol)-(abs(mean(z)-(0))<tol);
Parameters.IncidentField=...                     % Incident field for absorbing boundary conditions
 @(x,y,z,t) [-ny(y).*(ny(y).*t.*(x.*y+x.*z+x.^2./2)+nx(x).*t.*((x.*y)./2+(y.*z)./2+y.^2./4))-nz(z).*(nz(z).*t.*(x.*y+x.*z+x.^2./2)+nx(x).*t.*((x.*z)./2+(y.*z)./2+z.^2./4))-ny(y).*t.^2.*(x./2+y./4)-nz(z).*t.^2.*(x./2+z./4),...
             +nx(x).*(ny(y).*t.*(x.*y+x.*z+x.^2./2)+nx(x).*t.*((x.*y)./2+(y.*z)./2+y.^2./4))+nz(z).*(nz(z).*t.*((x.*y)./2+(y.*z)./2+y.^2./4)-ny(y).*t.*((x.*z)./2+(y.*z)./2+z.^2./4))+nx(x).*t.^2.*(x./2+y./4)+nz(z).*t.^2.*(y./4-z./4),...
             +nx(x).*(nz(z).*t.*(x.*y+x.*z+x.^2./2)+nx(x).*t.*((x.*z)./2+(y.*z)./2+z.^2./4))-ny(y).*(nz(z).*t.*((x.*y)./2+(y.*z)./2+y.^2./4)-ny(y).*t.*((x.*z)./2+(y.*z)./2+z.^2./4))+nx(x).*t.^2.*(x./2+z./4)-ny(y).*t.^2.*(y./4-z./4)];
clear tol nx ny nz
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Mesh.File={'Mesh_cube_1'};                       % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='yes';                        % Time dependent problem
Time.InitialTime=0;                              % Initial time
Time.FinalTime=1;                                % Final time
Time.TimeStepSize=1;                             % Time step size
Time.BDFOrder=3;                                 % BDF order
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=[1,2,3];                    % Dirichlet portion
Boundaries.Absorbing=[4,5,6];                    % Absorbing portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeError=...                         % Compute error
  {'MagneticField','L2';
   'ElectricField','L2';
   'ElectricFieldPost','L2';
   'ElectricFieldPost','Hcurl'};
Options.Test=...                                 % Test
  ['Results.MagneticFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldErrorL2       <1e-12 && ',...
   'Results.ElectricFieldPostErrorL2   <1e-12 && ',...
   'Results.ElectricFieldPostErrorHcurl<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main