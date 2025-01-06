function [Ionization]=solveSahaEquation(...
         Parameters,IndependentVariable,DensityOrPressure,PressureOrTemperature)
         % Solve Saha equation

% Physical constants
me=Parameters.ElectronMass;
h=Parameters.PlanckConstant;
kb=Parameters.BoltzmannConstant;
NA=Parameters.AvogadroNumber;
e=Parameters.ElectronCharge;
eps0=Parameters.ElectricPermittivity;

% Number of species
J=length(Parameters.MolarFraction);

% Material-specific quantities
c=Parameters.MolarFraction;
M=Parameters.MolarMass;
I=Parameters.IonizationEnergy;
Z=Parameters.MaxIonization;
if strcmp(Parameters.PartitionFunctionsExcitedStates,'yes')
  g=Parameters.StatisticalWeight;
  E=Parameters.EnergyLevel;
else
  g0=Parameters.StatisticalWeightsGroundState;
  U=cell(J,1);
  for j=1:J
    for r=1:Z(j)+1
      U{j}{r}=@(Zav) g0{j}(r);
    end
  end
end
if strcmp(Parameters.SahaNonIdeal,'yes') && ...
   isa(Parameters.IonizationEnergyLowering,'function_handle')
  DI_fun=@(nE,r) Parameters.IonizationEnergyLowering(nE,r);
else
  DI_fun=@(nE,r) 0;
end

% Get independent variables
switch IndependentVariable
  case 'DensityTemperature'
    rho=DensityOrPressure;
    T=PressureOrTemperature;
  case 'PressureTemperature'
    P=DensityOrPressure;
    T=PressureOrTemperature;
  case 'DensityPressure'
    rho=DensityOrPressure;
    P=PressureOrTemperature;
end

% Algorithm parameters
tolL=Parameters.SahaToleranceLinear;
opts=optimset('TolX',tolL);
if strcmp(Parameters.SahaNonIdeal,'yes') && ...
   (not(isa(Parameters.IonizationEnergyLowering,'function_handle')) || ...
    not(strcmp(Parameters.PressureCorrection,'None')))
  tolN=Parameters.SahaToleranceNonLinear;
  itmax=Parameters.SahaMaxIterations;
  omega=Parameters.SahaRelaxation;
end

% Linear step --------------------------------------------------------------------------------------
% Compute total number density of heavy particles as a function of the unknown average charge state
switch IndependentVariable
  case 'DensityTemperature'
    nH=@(Zav) NA/sum(c.*M)*rho;
  case 'PressureTemperature'
    nH=@(Zav) P/((1+Zav)*kb*T);
  case 'DensityPressure'
    nH=@(Zav) NA/sum(c.*M)*rho;
end

% Define total number density of electrons as a function of the unknown average charge state
nE=@(Zav) Zav*nH(Zav);

% Compute temperature as a function of the unknown average charge state
switch IndependentVariable
  case 'DensityTemperature'
    T=@(Zav) T;
  case 'PressureTemperature'
    T=@(Zav) T;
  case 'DensityPressure'
    T=@(Zav) P/((1+Zav)*nH(Zav)*kb);
end

% Compute empirical lowering of ionization energy as a function of the unknown average charge state
DI=cell(1,max(Z));
for r=1:max(Z)
  DI{r}=@(Zav) DI_fun(nE(Zav),r);
end

% Compute species partion functions
if strcmp(Parameters.PartitionFunctionsExcitedStates,'yes')
  U=cell(J,1);
  for j=1:J
    for r=1:Z(j)
      nstar=find(E{j}{r}<=I{j}(r),1,'last');
      U{j}{r}=@(Zav) sum(g{j}{r}(1:nstar).*exp(-E{j}{r}(1:nstar)/(kb*T(Zav))));
    end
    U{j}{end+1}=@(Zav) 1; % not sure about this
  end
end

% Define species Saha function as a function of the unknown average charge state
f=cell(J,max(Z)+1);
f0max=1;
for j=1:J
  for r=1:Z(j)
    f{j,r+1}=@(Zav) 2*U{j}{r+1}(Zav)/U{j}{r}(Zav)*(2*pi*me*kb*T(Zav)/h^2)^(3/2)*...
                    exp(-(I{j}(r)-DI{r}(Zav))/(kb*T(Zav)));
    f0=f{j,r+1}(0);
    if f0>f0max
      f0max=f0;
    end
  end
end

% Compute scaling factor to avoid overflow when multiplying the species Saha function
s=max(1,f0max^(max(Z)/2));

% Compute average charge state by solving a trascendental equation
Ze=cell(J,1);
fun=@(Zav) 0;
for j=1:J
  aux1=@(Zav) 0;
  aux2=@(Zav) 0;
  for i=1:Z(j)
    aux0=@(Zav) 1;
    for m=1:i
      aux0=@(Zav) aux0(Zav)*f{j,m+1}(Zav)/s^(1/i);
    end
    aux1=@(Zav) aux1(Zav)+i*aux0(Zav)/(Zav*nH(Zav))^i;
    aux2=@(Zav) aux2(Zav)+1*aux0(Zav)/(Zav*nH(Zav))^i;
  end
  Ze{j}=@(Zav) c(j)*aux1(Zav)/(1/s+aux2(Zav));
  fun=@(Zav) fun(Zav)+Ze{j}(Zav);
end
fun=@(Zav) Zav-fun(Zav);
if     fun(tolL)>0      % No ionization
  Zav=tolL;
elseif fun(sum(c.*Z))<0 % Full ionization
  Zav=sum(c.*Z);
else                    % Partial ionization
  Zav=fzero(@(Zav) fun(Zav),[tolL,sum(c.*Z)],opts);
end

% Compute species contribution to average charge state
Ze=feval(@(Zav) cellfun(@(Ze) Ze(Zav),Ze),Zav); %#ok

% Compute total number density of heavy particles
nH=nH(Zav);

% Compute total number density of electrons
nE=nE(Zav);

% Compute temperature
T=T(Zav);

% Compute species composition of various ionization states
alpha=zeros(J,max(Z)+1);
for j=1:J
  aux=0;
  for i=1:Z(j)
    aux0=1;
    for m=1:i
      aux0=aux0*f{j,m+1}(Zav);
    end
    aux=aux+i*aux0/(Zav*nH)^i;
  end
  alpha(j,1)=Ze(j)/aux;
  for r=1:Z(j)
    alpha(j,r+1)=alpha(j,r)/(Zav*nH)*f{j,r+1}(Zav);
  end
end

% Compute species number density at various ionization states
n=alpha*nH;

% Compute species number density of heavy particles
nh=c*nH;

% Compute species number density of electrons
ne=Ze*nH;
% --------------------------------------------------------------------------------------------------

% Non-linear step ----------------------------------------------------------------------------------
% Initialization
it=0;
err=inf;
lambdaD=inf;

if strcmp(Parameters.SahaNonIdeal,'yes') && ...
   (not(isa(Parameters.IonizationEnergyLowering,'function_handle')) || ...
    not(strcmp(Parameters.PressureCorrection,'None')))
  % Iterative procedure
  while err>tolN && it<itmax
    
    % Update iteration
    it=it+1;
    
    % Save old Debye length
    lambdaD_old=lambdaD;
    
    % Compute de Broglie wavelength
    LambdaB=h/sqrt(2*pi*me*kb*T);
    
    % Compute lowering of ionization energy
    if strcmp(Parameters.IonizationEnergyLowering,'None')
      DI=zeros(1,max(Z));
    elseif strcmp(Parameters.IonizationEnergyLowering,'Griem')
      DI=(1:max(Z))*e^2/(4*pi*eps0*lambdaD);
    elseif strcmp(Parameters.IonizationEnergyLowering,'Ebeling')
      DI=(1:max(Z))*e^2/(4*pi*eps0*(lambdaD+LambdaB/8));
    end
    
    % Compute pressure correction
    if strcmp(Parameters.PressureCorrection,'None')
      DP=0;
    elseif strcmp(Parameters.PressureCorrection,'Griem')
      DP=kb*T/(24*pi*lambdaD^3);
    end
    
    % Compute total pressure from equation of state
    if strcmp(IndependentVariable,'DensityTemperature')
      P=(1+Zav)*nH*kb*T-DP;
    end
    
    % Compute temperature from equation of state
    if strcmp(IndependentVariable,'DensityPressure')
      T=(P+DP)/((1+Zav)*nH*kb);
    end
    
    % Compute species partion functions
    if strcmp(Parameters.PartitionFunctionsExcitedStates,'yes')
      U=cell(J,1);
      for j=1:J
        for r=1:Z(j)
          nstar=find(E{j}{r}<=I{j}(r)-DI(r),1,'last');
          U{j}{r}=@(Zav) sum(g{j}{r}(1:nstar).*exp(-E{j}{r}(1:nstar)/(kb*T)));
        end
        U{j}{end+1}=@(Zav) 1; % not sure about this
      end
    end
    
    % Define species modified Saha function
    eta=zeros(J,max(Z)+1);
    for j=1:J
      for r=1:Z(j)
        eta(j,r+1)=2*U{j}{r+1}(Zav)/U{j}{r}(Zav)*(2*pi*me/h^2)^(3/2)*(kb*T)^(5/2)*...
                   exp(-(I{j}(r)-DI(r))/(kb*T));
      end
    end
    etamax=max(max(eta));
    
    % Compute scaling factor to avoid overflow when multiplying the species modified Saha function
    s=max(1,etamax^(max(Z)/2));
    
    % Compute ratio of ideal total electron pressure to total pressure by solving a trascendental equation
    pe=cell(J,1);
    fun=@(Pe) 0;
    for j=1:J
      aux1=@(pE) 0;
      aux2=@(pE) 0;
      for i=1:Z(j)
        aux0=prod(eta(j,2:i+1)/s^(1/i));
        aux1=@(pE) aux1(pE)+i*aux0/(pE*P)^i;
        aux2=@(pE) aux2(pE)+1*aux0/(pE*P)^i;
      end
      pe{j}=@(pE) c(j)*(1-DP/P-pE)*aux1(pE)/(1/s+aux2(pE));
      fun=@(pE) fun(pE)+pe{j}(pE);
    end
    fun=@(pE) pE-fun(pE);
    if     fun(tolL)>0 % No ionization
      pE=tolL;
    elseif fun(1)<0    % Full ionization
      pE=1;
    else               % Partial ionization
      pE=fzero(@(pE) fun(pE),[tolL,1],opts);
    end
    
    % Compute species contribution to ratio of ideal total electron pressure to total pressure
    pe=feval(@(pE) cellfun(@(pe) pe(pE),pe),pE); %#ok
    
    % Compute total number density of heavy particles from equation of state
    nH=(P+DP)/((1+Zav)*kb*T);
    
    % Compute total number density of electrons
    nE=pE*P/(kb*T);
    
    % Compute species ratio of ideal pressure at various ionization states to total pressure
    p=zeros(J,max(Z)+1);
    for j=1:J
      aux=0;
      for i=1:Z(j)
        aux=aux+i*prod(eta(j,2:i+1))/(pE*P)^i;
      end
      p(j,1)=pe(j)/aux;
      for r=1:Z(j)
        p(j,r+1)=p(j,r)/(pE*P)*eta(j,r+1);
      end
    end
    
    % Compute species number density at various ionization states
    n=p*P/(kb*T);
    
    % Compute average charge state
    Zav=pE/(1-DP/P-pE); % or Zav=nE/sum(n,1:2); or Zav=pE/sum(p,1:2);
    
    % Compute species contribution to average charge state
    Ze=pe/(1-DP/P-pE); % or Ze=ne/sum(n,1:2); or Ze=pe/sum(p,1:2);
    
    % Compute species composition of various ionization states
    alpha=p/(1-DP/P-pE); % or alpha=n/sum(n,1:2); or alpha=p/sum(p,1:2);
    
    % Compute species number density of heavy particles
    nh=c*nH;
    
    % Compute species number density of electrons
    ne=pe*P/(kb*T);
    
    % Compute Debye length
    if strcmp(Parameters.DebyeLengthDefinition,'Incomplete')
      lambdaD=sqrt(eps0*kb*T/(e^2*nE));
    elseif strcmp(Parameters.DebyeLengthDefinition,'Complete')
      lambdaD=sqrt(eps0*kb*T/(e^2*(nE+sum((0:max(Z)).^2*n'))));
    end
    
    % Relaxation of Debye length
    if it>1
      lambdaD=omega*lambdaD+(1-omega)*lambdaD_old;
    end
    
    % Compute error
    err=abs(lambdaD/lambdaD_old-1);
  end
end
% --------------------------------------------------------------------------------------------------

% Store ionization info
Ionization.AverageChargeState=Zav;
Ionization.Composition=alpha;
Ionization.Iterations=it;
Ionization.TotalNumberDensityHeavyParticles=nH;
Ionization.TotalNumberDensityElectron=nE;
Ionization.SpeciesAverageChargeState=Ze;
Ionization.SpeciesNumberDensityIonization=n;
Ionization.SpeciesNumberDensityHeavyParticles=nh;
Ionization.SpeciesNumberDensityElectron=ne;

end