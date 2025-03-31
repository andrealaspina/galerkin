clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Structural';                 % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='ElasticityModal_CG';     % Formulation
Parameters.Problem='Structural';                 % Problem
Parameters.Model='LinearElasticity';             % Model
Parameters.Degree=2;                             % Degree
Parameters.Density=1;                            % Density
Parameters.YoungsModulus=1;                      % Young's modulus
Parameters.PoissonsRatio=0;                      % Poisson's ratio
Parameters.Displacement=@(x,y,z,n) [0*x,0*x,0*x];% Displacement
Parameters.analyticalSolution=...                % Analytical solution
  @(Parameters,Mesh,Time) analyticalSolution(Parameters,Mesh,Time);
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
Geometry=multicuboid(30,1,1).translate([15,0,0]);% Geometry definition
Mesh.MinElementSize=1;                           % Minimum element size
Mesh.MaxElementSize=1;                           % Maximum element size
Mesh.MeshGradation=1.5;                          % Element size growth rate
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='yes';                   % Symmetrize matrix
System.Nonlinear='no';                           % Nonlinear
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
Time.Mode=[1,2,4];                               % Mode of modal analysis
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='eigs';                              % Type
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Fixed=5;                              % Fixed portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.ComputeQuantityStart=...                 % Compute quantity of interest (start of simulation)
  '[Parameters]=Parameters.analyticalSolution(Parameters,Mesh,Time);';
Options.ComputeQuantityError=...                 % Compute quantity of interest (error computation)
  ['for iM=1:numel(Time.Mode)',...
   '  if max(Results.Displacement(:,2,iM))>=max(Results.Displacement(:,3,iM));',...
   '    Results.MainDisplacement(:,:,iM)=Results.Displacement(:,2,iM);',...
   '  else;',...
   '    Results.MainDisplacement(:,:,iM)=Results.Displacement(:,3,iM);',...
   '  end;',...
   'end;'];
Options.ComputeError=...                         % Compute error
  {'MainDisplacement','L2';
   'NaturalFrequency','Number'};
Options.Test=...                                 % Test
  ['abs(Results.MainDisplacementErrorL2    -1.007416330309386e-02)<1e-12 && ',...
   'abs(Results.NaturalFrequencyErrorNumber-5.995651270992140e-03)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

main

%% Natural frequency

function [Parameters]=analyticalSolution(Parameters,Mesh,Time)

% Extract parameters
tol=1e-6;
d=1e-2;
E=Parameters.YoungsModulus;
r=Parameters.Density;
L=max(Mesh.Nodes(1,:))-min(Mesh.Nodes(1,:));
B=max(Mesh.Nodes(2,:))-min(Mesh.Nodes(2,:));
H=max(Mesh.Nodes(3,:))-min(Mesh.Nodes(3,:));

% Find roots of characteristic equation
beta=zeros(1,max(Time.Mode));
beta(1)=fzero(@(beta) cos(beta)*cosh(beta)+1,1);
for iM=2:max(Time.Mode)
  sol=beta(iM-1);
  guess=beta(iM-1);
  while abs(sol-beta(iM-1))<tol
    guess=guess+d;
    sol=fzero(@(beta) cos(beta)*cosh(beta)+1,guess);
  end
  beta(iM)=sol;
end

% Frequency-dependent coefficients
beta=@(n) beta(floor((n+1)/2));
C=@(n) 1/(cos(beta(n))-cosh(beta(n))+...
         (sin(beta(n))-sinh(beta(n)))./...
         (cos(beta(n))+cosh(beta(n))).*...
         (sin(beta(n))-sinh(beta(n))));

% Analytical solution
Parameters.NaturalFrequency=@(n) beta(n).^2*sqrt((E*B*H^3/12)/(r*B*H*L^4));
Parameters.MainDisplacement=@(x,y,z,n)  C(n)*(cos(beta(n)*x/L)-cosh(beta(n)*x/L)+...
                                             (sin(beta(n)    )-sinh(beta(n)    ))./...
                                             (cos(beta(n)    )+cosh(beta(n)    )).*...
                                             (sin(beta(n)*x/L)-sinh(beta(n)*x/L)));

end