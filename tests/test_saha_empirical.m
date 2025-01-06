clearvars -except Files Passed CPUTime Test;

%% Input

% Parameters ---------------------------------------------------------------------------------------
Parameters.ElectronMass=9.1093837015e-31;        % Electron mass
Parameters.PlanckConstant=6.62607015e-34;        % Planck’s constant
Parameters.BoltzmannConstant=1.380649e-23;       % Boltzmann’s constant
Parameters.AvogadroNumber=6.02214076e23;         % Avogadro's number
Parameters.ElectronCharge=1.602176634e-19;       % Electron charge
Parameters.ElectricPermittivity=8.8541878128e-12;% Electric permittivity
Parameters.MolarFraction=1;                      % Molar fraction
Parameters.MolarMass=26.981539e-3;               % Molar mass
Parameters.IonizationEnergy=...                  % Ionization energy
  {[577.5,1816.7,2744.8,11577,14842,18379,23326,27465,31853,38473]*1e3/Parameters.AvogadroNumber};
Parameters.PartitionFunctionsExcitedStates='no'; % Include excited states in the partition functions
Parameters.StatisticalWeightsGroundState=...     % Statistical weights of the ground state
  {[6,1,2,1,6,9,4]};
Parameters.MaxIonization=6;                      % Maximum ionization
Parameters.SahaToleranceLinear=eps;              % Linear tolerance for Saha equation
Parameters.SahaNonIdeal='yes';                   % Non-ideality for Saha equation
Parameters.IonizationEnergyLowering=...          % Lowering of ionization energy
  @(nE,r) (6.96e-7)*(nE/1e6).^(1/3)*r^(2/3)*Parameters.ElectronCharge;
Parameters.PressureCorrection='None';            % Pressure correction
% --------------------------------------------------------------------------------------------------

% State variables ----------------------------------------------------------------------------------
IndependentVariables='DensityTemperature';       % Independent variables
Density=100;                                     % Density
Temperature=30000;                               % Temperature
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.Test=...                                 % Test
  ['abs(Results.AverageChargeState-2.299186769538464e-00)<1e-12 && ',...
   'abs(Results.Composition(1)    -6.425324059989132e-03)<1e-12 && ',...
   'abs(Results.Composition(2)    -5.374163629454972e-02)<1e-12 && ',...
   'abs(Results.Composition(3)    -5.740539856924786e-01)<1e-12 && ',...
   'abs(Results.Composition(4)    -3.657790539529738e-01)<1e-12 && ',...
   'abs(Results.Composition(5)    -9.028653893591346e-15)<1e-12 && ',...
   'abs(Results.Composition(6)    -7.521937211855238e-34)<1e-12 && ',...
   'abs(Results.Composition(7)    -7.457067750926234e-59)<1e-12 && ',...
   'abs(Results.Iterations        -0.000000000000000e-00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

Results=solveSahaEquation(Parameters,IndependentVariables,Density,Temperature);

fprintf('\n--------------------------------------------------------------------------------')

% Get file name
FileName=mfilename;

% Print file name
fprintf('\nInput file: %s',FileName);

fprintf('\n--------------------------------------------------------------------------------');
fprintf('\n');

% Perform test
Simulation.TestPassed=eval(Options.Test);
if Simulation.TestPassed; fprintf('\nTEST PASSED\n');
else;                     fprintf('\nTEST FAILED\n'); end