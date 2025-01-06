clearvars -except Files Passed CPUTime Test;

%% Input

% Simulation ---------------------------------------------------------------------------------------
Simulation.Type='SingleSimulation';              % Simulation type
Simulation.Problem='Darcy';                      % Problem
% --------------------------------------------------------------------------------------------------

% Parameters ---------------------------------------------------------------------------------------
Parameters.Formulation='Darcy_CG';               % Formulation
Parameters.Problem='Darcy';                      % Problem
Parameters.UnitCellProblem='Navier-Stokes';      % Unit cell problem
Parameters.UnitCellFile='tests/UnitCell';        % Unit cell file
Parameters.Degree=1;                             % Degree
Parameters.NitschePenalty=100;                   % Nitsche's penalty parameter
Parameters.RelativeFiniteDifference=1e-6;        % Relative finite difference
Parameters.Pressure=@(x,y,z,t) 0*x;              % Pressure
Parameters.Flux=@(x,y,z,t) 10*(x<1e-6);          % Flux
Parameters.Source=@(x,y,z,t) 0*x;                % Source
Parameters.Permeability=...                      % Permeability
  @(Dpx,Dpy,Parameters,UnitCell) computePermeability(Dpx,Dpy,Parameters,UnitCell);
% --------------------------------------------------------------------------------------------------

% Geometry and mesh --------------------------------------------------------------------------------
MeshFile={'Mesh_test_2d_1elem_k1'};              % Mesh file
% --------------------------------------------------------------------------------------------------

% System -------------------------------------------------------------------------------------------
System.SymmetrizeMatrix='no';                    % Symmetrize matrix
System.Nonlinear='yes';                          % Nonlinear
System.Tolerance=1e-3;                           % Tolerance
System.MaxIterations=10;                         % Maximum number of iterations
% --------------------------------------------------------------------------------------------------

% Time parameters ----------------------------------------------------------------------------------
Time.TimeDependent='no';                         % Time dependent problem
% --------------------------------------------------------------------------------------------------

% Solver -------------------------------------------------------------------------------------------
Solver.Type='backslash';                         % Type
Solver.Equilibrate='yes';                        % Equilibrate matrix
% --------------------------------------------------------------------------------------------------

% Boundary splitting -------------------------------------------------------------------------------
Boundaries.Dirichlet=1;                          % Dirichlet portion
Boundaries.Neumann=[2,3];                        % Neumann portion
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.PlotGeometry='no';                       % Plot geometry
Options.PlotMesh='no';                           % Plot mesh
Options.Export2Paraview='no';                    % Export to Paraview
Options.ComputeError={'Pressure'};               % Compute error
Options.Test=...                                 % Test
  'abs(Results.PressureErrorL2-2.108084213203400e+06)<1e-07';
% --------------------------------------------------------------------------------------------------

%% Unit cell

PoolAux=gcp; if isempty(PoolAux); PoolAux=[]; PoolAux.NumWorkers=1; end
parfor iP=1:PoolAux.NumWorkers
  % Connect Comsol
  connectComsol(2);
  import com.comsol.model.*
  import com.comsol.model.util.*

  % Create model
  UnitCell = ModelUtil.create('Model');
  UnitCell.baseSystem('none');

  % Initialize model
  UnitCell.component.create('comp1', true);
  UnitCell.component('comp1').geom.create('geom1', 2);
  UnitCell.component('comp1').mesh.create('mesh1');
  UnitCell.component('comp1').physics.create('spf', 'LaminarFlow', 'geom1');
  UnitCell.study.create('std1');
  UnitCell.study('std1').create('stat', 'Stationary');
  UnitCell.study('std1').feature('stat').activate('spf', true);

  % Create geometry
  UnitCell.component('comp1').geom('geom1').create('sq1', 'Square');
  UnitCell.component('comp1').geom('geom1').feature('sq1').set('size', 0.01);
  UnitCell.component('comp1').geom('geom1').create('c1', 'Circle');
  UnitCell.component('comp1').geom('geom1').feature('c1').set('pos', {'0.01/2' '0.01/2'});
  UnitCell.component('comp1').geom('geom1').feature('c1').set('r', '0.01/5');
  UnitCell.component('comp1').geom('geom1').create('dif1', 'Difference');
  UnitCell.component('comp1').geom('geom1').feature('dif1').selection('input').set({'sq1'});
  UnitCell.component('comp1').geom('geom1').feature('dif1').selection('input2').set({'c1'});
  UnitCell.component('comp1').geom('geom1').run('fin');

  % Set parameters
  UnitCell.component('comp1').physics('spf').prop('PhysicalModelProperty').set('pref', 0);
  UnitCell.component('comp1').physics('spf').prop('PhysicalModelProperty').set('Tref', 0);
  UnitCell.component('comp1').physics('spf').feature('fp1').set('rho_mat', 'userdef');
  UnitCell.component('comp1').physics('spf').feature('fp1').set('mu_mat', 'userdef');
  UnitCell.component('comp1').physics('spf').feature('fp1').set('rho', 1000);
  UnitCell.component('comp1').physics('spf').feature('fp1').set('mu', 3);

  % Set boundary conditions
  UnitCell.component('comp1').physics('spf').create('pfc1', 'PeriodicFlowCondition', 1);
  UnitCell.component('comp1').physics('spf').feature('pfc1').selection.set([2 3]);
  UnitCell.component('comp1').physics('spf').create('pfc2', 'PeriodicFlowCondition', 1);
  UnitCell.component('comp1').physics('spf').feature('pfc2').selection.set([1 4]);
  UnitCell.component('comp1').physics('spf').create('prpc1', 'PressurePointConstraint', 0);
  UnitCell.component('comp1').physics('spf').feature('prpc1').selection.set([1 2 7 8]);
  UnitCell.component('comp1').physics('spf').create('vf1', 'VolumeForce', 2);
  UnitCell.component('comp1').physics('spf').feature('vf1').set('F', [0 0 0]);
  UnitCell.component('comp1').physics('spf').feature('vf1').selection.set(1);

  % Set mesh
  UnitCell.component('comp1').mesh('mesh1').feature('size').set('custom', true);
  UnitCell.component('comp1').mesh('mesh1').feature('size').set('hmax', '0.01/8');
  UnitCell.component('comp1').mesh('mesh1').feature('size').set('hmin', '0.01/8');
  UnitCell.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1.5);
  UnitCell.component('comp1').mesh('mesh1').create('dis1', 'Distribution');
  UnitCell.component('comp1').mesh('mesh1').feature('dis1').selection.set([1 2 3 4]);
  UnitCell.component('comp1').mesh('mesh1').feature('dis1').set('numelem', 8);
  UnitCell.component('comp1').mesh('mesh1').create('dis2', 'Distribution');
  UnitCell.component('comp1').mesh('mesh1').feature('dis2').selection.set([5 6 7 8]);
  UnitCell.component('comp1').mesh('mesh1').feature('dis2').set('numelem', 4);
  UnitCell.component('comp1').mesh('mesh1').create('edg1', 'Edge');
  UnitCell.component('comp1').mesh('mesh1').feature('edg1').selection.all;
  UnitCell.component('comp1').mesh('mesh1').create('fq1', 'FreeQuad');
  UnitCell.component('comp1').mesh('mesh1').create('ref1', 'Refine');
  UnitCell.component('comp1').mesh('mesh1').feature('ref1').set('numrefine', 1);
  UnitCell.component('comp1').mesh('mesh1').run;

  % Initialize solution
  UnitCell.sol.create('sol1');
  UnitCell.sol('sol1').study('std1');

  % Initialize study
  UnitCell.study('std1').feature('stat').set('notlistsolnum', 1);
  UnitCell.study('std1').feature('stat').set('notsolnum', '1');
  UnitCell.study('std1').feature('stat').set('listsolnum', 1);
  UnitCell.study('std1').feature('stat').set('solnum', '1');

  % Set solver
  UnitCell.sol('sol1').create('st1', 'StudyStep');
  UnitCell.sol('sol1').feature('st1').set('study', 'std1');
  UnitCell.sol('sol1').feature('st1').set('studystep', 'stat');
  UnitCell.sol('sol1').create('v1', 'Variables');
  UnitCell.sol('sol1').feature('v1').set('control', 'stat');
  UnitCell.sol('sol1').create('s1', 'Stationary');
  UnitCell.sol('sol1').feature('s1').set('stol', '1e-12');
  UnitCell.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
  UnitCell.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('dtech', 'auto');
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-4);
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
  UnitCell.sol('sol1').feature('s1').create('d1', 'Direct');
  UnitCell.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
  UnitCell.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
  UnitCell.sol('sol1').feature('s1').feature('d1').label('Direct, fluid flow variables (spf)');
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('dtech', 'auto');
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-4);
  UnitCell.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
  UnitCell.sol('sol1').feature('s1').feature.remove('fcDef');
  UnitCell.sol('sol1').attach('std1');

  % Set computation of average velocity
  UnitCell.result.numerical.create('av1', 'AvSurface');
  UnitCell.result.numerical('av1').set('intvolume', true);
  UnitCell.result.numerical('av1').selection.set(1);
  UnitCell.result.numerical('av1').setIndex('expr', 'u*(0.01^2-pi*(0.01/5)^2)/0.01^2', 0);
  UnitCell.result.numerical.create('av2', 'AvSurface');
  UnitCell.result.numerical('av2').set('intvolume', true);
  UnitCell.result.numerical('av2').selection.set(1);
  UnitCell.result.numerical('av2').setIndex('expr', 'v*(0.01^2-pi*(0.01/5)^2)/0.01^2', 0);

  % Compute linear permeability
  UnitCell.component('comp1').physics('spf').prop('PhysicalModelProperty').set('StokesFlowProp', true);
  UnitCell.component('comp1').physics('spf').feature('vf1').set('F',[1,1,0]);
  UnitCell.sol('sol1').runAll;
  PermeabilityLinearXXAux(iP)=UnitCell.result.numerical('av1').getReal;
  PermeabilityLinearYYAux(iP)=UnitCell.result.numerical('av2').getReal;

  % Set Stokes/Navier-Stokes unit cell problem
  if     strcmp(Parameters.UnitCellProblem,'Stokes')
    UnitCell.component('comp1').physics('spf').prop('PhysicalModelProperty').set('StokesFlowProp', true);
  elseif strcmp(Parameters.UnitCellProblem,'Navier-Stokes')
    UnitCell.component('comp1').physics('spf').prop('PhysicalModelProperty').set('StokesFlowProp', false);
  end

  % Plot unit cell mesh
  figure('Color','w');
  mphmesh(UnitCell);
  pause(eps);

  % Save unit cell
  mphsave(UnitCell,sprintf('%s%d',Parameters.UnitCellFile));
end

% Transfer unit cell to workers
Parameters.UnitCell=parallel.pool.Constant(@() mphopen(Parameters.UnitCellFile));

% Store linear permeability
Parameters.PermeabilityLinear(1)=PermeabilityLinearXXAux(1);
Parameters.PermeabilityLinear(2)=PermeabilityLinearYYAux(1);

%% Main

main

% Disconnect Comsol
PoolAux=gcp; if isempty(PoolAux); PoolAux=[]; PoolAux.NumWorkers=1; end
parfor iP=1:PoolAux.NumWorkers
  ModelUtil.disconnect
end

%% Compute permeability

function [Kappaxx,Kappayy]=computePermeability(Dpx,Dpy,Parameters,UnitCell)

% Compute permeability
UnitCell.component('comp1').physics('spf').feature('vf1').set('F',[Dpx,Dpy,0]);
UnitCell.sol('sol1').run;
Ux=UnitCell.result.numerical('av1').getReal;
Uy=UnitCell.result.numerical('av2').getReal;
Kappaxx=Ux/Dpx;
Kappayy=Uy/Dpy;
if     norm(Dpx)<norm(Dpy)*1e-6
  Kappaxx=Parameters.PermeabilityLinear(1);
elseif norm(Dpy)<norm(Dpx)*1e-6
  Kappayy=Parameters.PermeabilityLinear(2);
end

end