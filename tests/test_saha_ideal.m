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
Parameters.MolarMass=39.948e-3;                  % Molar mass
Parameters.IonizationEnergy=...                  % Ionization energy
  {[1520.6,2665.8,3931,5771,7238,8781,11995,13842,40760,46186]*1e3/Parameters.AvogadroNumber};
Parameters.PartitionFunctionsExcitedStates='no'; % Include excited states in the partition functions
Parameters.StatisticalWeightsGroundState=...     % Statistical weights of the ground state
  {[1,6,9,4,9,6,1]};
Parameters.MaxIonization=6;                      % Maximum ionization
Parameters.SahaToleranceLinear=eps;              % Linear tolerance for Saha equation
Parameters.SahaNonIdeal='no';                    % Non-ideality for Saha equation
% --------------------------------------------------------------------------------------------------

% State variables ----------------------------------------------------------------------------------
IndependentVariables='PressureTemperature';      % Independent variables
Pressure=101325;                                 % Pressure
Temperature=50000;                               % Temperature
% --------------------------------------------------------------------------------------------------

% Output options -----------------------------------------------------------------------------------
Options.Test=...                                 % Test
  ['abs(Results.AverageChargeState-3.461810360383030e-00)<1e-12 && ',...
   'abs(Results.Composition(1)    -3.426171348955255e-10)<1e-12 && ',...
   'abs(Results.Composition(2)    -2.513744405248731e-05)<1e-12 && ',...
   'abs(Results.Composition(3)    -2.933690729260222e-02)<1e-12 && ',...
   'abs(Results.Composition(4)    -4.836311208047152e-01)<1e-12 && ',...
   'abs(Results.Composition(5)    -4.828163468237952e-01)<1e-12 && ',...
   'abs(Results.Composition(6)    -4.190265108859138e-03)<1e-12 && ',...
   'abs(Results.Composition(7)    -2.221833586136860e-07)<1e-12 && ',...
   'abs(Results.Iterations        -0.000000000000000e-00)<1e-12'];
% --------------------------------------------------------------------------------------------------

%% Main

Results=solveSahaEquation(Parameters,IndependentVariables,Pressure,Temperature);

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