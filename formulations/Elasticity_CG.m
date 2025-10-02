classdef Elasticity_CG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim;

    % Discretization type
    DiscretizationType='CG';

    % Time derivative order
    TimeDerOrder=2;
    
    % Domain
    Domain='Time';
    
  end
  
  methods
    
    %% Initialize unknowns
    function [Block]=initializeUnknowns(~,iD,Block,Parameters,Time,Sizes)
      Block(iD,iD).SolutionGlobal=zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,1);
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionOld=zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,...
          Parameters(iD).TimeDerOrder);
      end
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,~)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionGlobal=Parameters(iD).Displacement(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime);
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionGlobal;
        Block(iD,iD).SolutionOld(:,:,2)=Parameters(iD).Displacement(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',-Time.TimeStepSize); % FIX THIS!!!
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=Parameters(iD).Displacement(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld);
          end
        end
      end
    end

    %% Evaluate solution at fixed DOFs
    function [Block]=evaluateSolutionFixedDofs(~,iD,Block,Parameters,Mesh,Time,Sizes)
      Block(iD,iD).SolutionGlobal(Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp),:)=...
        Parameters(iD).Displacement(...
        Mesh(iD).Nodes(1,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(2,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(3,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',Time.Time);
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,Simulation,Parameters,~,Faces,Time,...
        RefElement,Sizes)
      NodesElem=Elements(iD1).Nodes;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      
      % Extract coupling data
      SolutionGlobalElemCoupled=double.empty(Sizes.NumElements,0);
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        if strcmp(Simulation.Problem,'Structural')
          iD2=setdiff(1:2,iD1);
        elseif strcmp(Simulation.Problem,'FluidStructureInteraction')
          if strcmp(Parameters(iD1).Problem,'Mesh')
            iD2=find(contains({Parameters.Problem},'Structural'));
          elseif strcmp(Parameters(iD1).Problem,'Structural')
            iD2=find(contains({Parameters.Problem},'Fluid'));
          end
        end
        SolutionGlobalElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        iEF=sub2ind([Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces],...
          Faces(iD1,iD2).Interface(:,1),Faces(iD1,iD2).Interface(:,2));
        SolutionGlobalElemCoupled(iEF)=Elements(iD2).SolutionGlobal(Faces(iD1,iD2).Interface(:,3));
      end
      
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem]=...
          buildBlockElement(iD1,NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionOldElem{iElem},...
          SolutionGlobalElemCoupled(iElem,:),...
          Simulation,Parameters,Time,RefElement,Sizes);
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
      end
      Block(iD1,iD1).LhsGlobal=fsparse(Block(iD1,iD1).LhsRowIndices,...
                                       Block(iD1,iD1).LhsColIndices,LhsCoef(:));
      Block(iD1,iD1).RhsGlobal=fsparse(Block(iD1,iD1).RhsRowIndices,1,RhsCoef(:));
    end
    
    %% Do coupling
    function [Block]=doCoupling(~,iD1,iD2,Block,~,Simulation,Parameters,Mesh,Faces,~,RefElement,...
        Sizes)
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface)) && ...
         (strcmp(Simulation.Problem,'Structural') || ...
          (strcmp(Simulation.Problem,'FluidStructureInteraction') && ...
           ((strcmp(Parameters(iD1).Problem,'Mesh') && ...
             strcmp(Parameters(iD2).Problem,'Structural')) || ...
            (strcmp(Parameters(iD1).Problem,'Structural') && ...
             strcmp(Parameters(iD2).Problem,'Fluid')))))
        LhsCoupCoef=zeros(Sizes(iD1).NumElementLhsCoupCoef(iD2),Sizes(iD1).NumFacesInterface(iD2));
        for iFaceInterface=1:Sizes(iD1).NumFacesInterface(iD2)
          [LhsCoupElem]=...
            doCouplingElement(iFaceInterface,iD1,iD2,Block,...
            Simulation,Parameters,Mesh,Faces,RefElement,Sizes); 
          LhsCoupCoef(:,iFaceInterface)=reshape(LhsCoupElem',[],1);
        end
        Block(iD1,iD2).LhsGlobal=fsparse(Block(iD1,iD2).LhsRowIndices,...
                                         Block(iD1,iD2).LhsColIndices,[LhsCoupCoef(:);0]);
      else
        Block(iD1,iD2).LhsGlobal=fsparse(Block(iD1,iD2).LhsRowIndices,...
                                         Block(iD1,iD2).LhsColIndices,0);
      end
    end
    
    %% Store results
    function [Results]=storeResults(~,iD,iST,Results,Block,~,~,~,Time,~,~)
      if iST==1
        Results(iD).Time=[];
        Results(iD).Displacement=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).Displacement(:,:,iST)=Block(iD,iD).SolutionGlobal;
    end

    %% Evaluate stress
    function [Results]=evaluateStress(~,iD,Results,Elements,Parameters,Mesh,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      ResultsElem=mat2cell(...
        Results(iD).Displacement(Mesh(iD).Elements(:),:),...
        ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes,...
        Sizes(iD).NumSpaceDim);
      Stress=zeros(Sizes(iD).NumElementNodes*Sizes(iD).NumSpaceDim^2,Sizes(iD).NumElements);
      parfor iElem=1:Sizes(iD).NumElements
        [StressElem]=...
          evaluateStressElement(iD,NodesElem{iElem},...
          ResultsElem{iElem},...
          Parameters,RefElement,Sizes);
        Stress(:,iElem)=reshape(StressElem',[],1);
      end
      Results(iD).StressDisc=reshape(Stress,Sizes(iD).NumSpaceDim^2,[])';
      for iC=1:Sizes(iD).NumSpaceDim^2
        Results(iD).Stress(:,iC)=accumarray(Mesh(iD).Elements(:),Results(iD).StressDisc(:,iC),[],...
          @mean);
      end
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal]=buildBlockElement(...
  iD1,Nodes,Faces,SolutionGlobal,SolutionOld,SolutionGlobalCoupled,...
  Simulation,Parameters,Time,RefElement,Sizes)

% Get general parameters
isStructural=strcmp(Simulation.Problem,'Structural');
isFSI=strcmp(Simulation.Problem,'FluidStructureInteraction');
if isFSI
  isMesh=strcmp(Parameters(iD1).Problem,'Mesh');
  isStructure=strcmp(Parameters(iD1).Problem,'Structural');
end
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes=Sizes(iD1).NumElementNodes;
NumElementFaces=Sizes(iD1).NumElementFaces;
rho=Parameters(iD1).Density;
uD=Parameters(iD1).Displacement;
tN=Parameters(iD1).Traction;
f=Parameters(iD1).Force;
Xe=Nodes';
Xem=sum(Xe(1:NumElementFaces,:),1)/NumElementFaces;
gamma=Parameters(iD1).NitschePenalty;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
  beta=Time.BDF2ndDerEff;
end

% Get solution
ue=reshape(SolutionGlobal,[],1);
if isTimeDependent
  uolde=reshape(SolutionOld,[],BDFo+1);
end

% Initialize lhs
Kuu=zeros(nsd*NumElementNodes,nsd*NumElementNodes);

% Initialize rhs
fu=zeros(nsd*NumElementNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,~,pinvNe]=mapShapeFunctions(1,RefElement(iD1,iD1),RefElement(iD1,iD1),Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
uxe=ue(ne1);
uye=ue(ne2);
if nsd==3
  uze=ue(ne3);
end
if isTimeDependent
  uoldxe=uolde(ne1,:);
  uoldye=uolde(ne2,:);
  if nsd==3
    uoldze=uolde(ne3,:);
  end
end
if nsd==2
  Fxxe=pinvNe*(Nex*uxe+1); Fxye=pinvNe*(Ney*uxe);
  Fyxe=pinvNe*(Nex*uye);   Fyye=pinvNe*(Ney*uye+1);
elseif nsd==3
  Fxxe=pinvNe*(Nex*uxe+1); Fxye=pinvNe*(Ney*uxe);   Fxze=pinvNe*(Nez*uxe);
  Fyxe=pinvNe*(Nex*uye);   Fyye=pinvNe*(Ney*uye+1); Fyze=pinvNe*(Nez*uye);
  Fzxe=pinvNe*(Nex*uze);   Fzye=pinvNe*(Ney*uze);   Fzze=pinvNe*(Nez*uze+1);  
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
uxeg=Ne*uxe;
uyeg=Ne*uye;
if nsd==3
  uzeg=Ne*uze;
end
if isTimeDependent
  uoldxeg=Ne*uoldxe;
  uoldyeg=Ne*uoldye;
  if nsd==3
    uoldzeg=Ne*uoldze;
  end
end
if nsd==2
  Fxxeg=Ne*Fxxe; Fxyeg=Ne*Fxye;
  Fyxeg=Ne*Fyxe; Fyyeg=Ne*Fyye;
elseif nsd==3
  Fxxeg=Ne*Fxxe; Fxyeg=Ne*Fxye; Fxzeg=Ne*Fxze;
  Fyxeg=Ne*Fyxe; Fyyeg=Ne*Fyye; Fyzeg=Ne*Fyze;
  Fzxeg=Ne*Fzxe; Fzyeg=Ne*Fzye; Fzzeg=Ne*Fzze;
end
feg=f(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
fxeg=feg(:,1);
fyeg=feg(:,2);
if nsd==3
  fzeg=feg(:,3);
end

% Compute basic matrices
NweT=(weg.*Ne)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
Me=NweT*Ne;
Me=(Me+Me')/2;

% Compute stress and its linearization
if nsd==2
  [seg,dsdFeg]=...
    computeStress('yes','yes','no',Parameters(iD1),Xem,...
    [Fxxeg,Fxyeg,Fyxeg,Fyyeg],Sizes(iD1));
  sxxeg=seg(:,1); sxyeg=seg(:,2);
  syxeg=seg(:,3); syyeg=seg(:,4);
  dsxxdFxxeg=dsdFeg(:,1);  dsxxdFxyeg=dsdFeg(:,2);  dsxxdFyxeg=dsdFeg(:,3);  dsxxdFyyeg=dsdFeg(:,4);
                           dsxydFxyeg=dsdFeg(:,6);  dsxydFyxeg=dsdFeg(:,7);  dsxydFyyeg=dsdFeg(:,8);
                                                    dsyxdFyxeg=dsdFeg(:,11); dsyxdFyyeg=dsdFeg(:,12);
                                                                             dsyydFyyeg=dsdFeg(:,16);

elseif nsd==3
  [seg,dsdFeg]=...
    computeStress('yes','yes','no',Parameters(iD1),Xem,...
    [Fxxeg,Fxyeg,Fxzeg,Fyxeg,Fyyeg,Fyzeg,Fzxeg,Fzyeg,Fzzeg],Sizes(iD1));
  sxxeg=seg(:,1); sxyeg=seg(:,2); sxzeg=seg(:,3);
  syxeg=seg(:,4); syyeg=seg(:,5); syzeg=seg(:,6);
  szxeg=seg(:,7); szyeg=seg(:,8); szzeg=seg(:,9);
  dsxxdFxxeg=dsdFeg(:,1);  dsxxdFxyeg=dsdFeg(:,2);  dsxxdFxzeg=dsdFeg(:,3);  dsxxdFyxeg=dsdFeg(:,4);  dsxxdFyyeg=dsdFeg(:,5);  dsxxdFyzeg=dsdFeg(:,6);  dsxxdFzxeg=dsdFeg(:,7);  dsxxdFzyeg=dsdFeg(:,8);  dsxxdFzzeg=dsdFeg(:,9);
                           dsxydFxyeg=dsdFeg(:,11); dsxydFxzeg=dsdFeg(:,12); dsxydFyxeg=dsdFeg(:,13); dsxydFyyeg=dsdFeg(:,14); dsxydFyzeg=dsdFeg(:,15); dsxydFzxeg=dsdFeg(:,16); dsxydFzyeg=dsdFeg(:,17); dsxydFzzeg=dsdFeg(:,18);
                                                    dsxzdFxzeg=dsdFeg(:,21); dsxzdFyxeg=dsdFeg(:,22); dsxzdFyyeg=dsdFeg(:,23); dsxzdFyzeg=dsdFeg(:,24); dsxzdFzxeg=dsdFeg(:,25); dsxzdFzyeg=dsdFeg(:,26); dsxzdFzzeg=dsdFeg(:,27);
                                                                             dsyxdFyxeg=dsdFeg(:,31); dsyxdFyyeg=dsdFeg(:,32); dsyxdFyzeg=dsdFeg(:,33); dsyxdFzxeg=dsdFeg(:,34); dsyxdFzyeg=dsdFeg(:,35); dsyxdFzzeg=dsdFeg(:,36);
                                                                                                      dsyydFyyeg=dsdFeg(:,41); dsyydFyzeg=dsdFeg(:,42); dsyydFzxeg=dsdFeg(:,43); dsyydFzyeg=dsdFeg(:,44); dsyydFzzeg=dsdFeg(:,45);
                                                                                                                               dsyzdFyzeg=dsdFeg(:,51); dsyzdFzxeg=dsdFeg(:,52); dsyzdFzyeg=dsdFeg(:,53); dsyzdFzzeg=dsdFeg(:,54);
                                                                                                                                                        dszxdFzxeg=dsdFeg(:,61); dszxdFzyeg=dsdFeg(:,62); dszxdFzzeg=dsdFeg(:,63);
                                                                                                                                                                                 dszydFzyeg=dsdFeg(:,71); dszydFzzeg=dsdFeg(:,72);
                                                                                                                                                                                                          dszzdFzzeg=dsdFeg(:,81);
end

% Compute lhs
if nsd==2
  Kuu(ne1,ne1)=+NwexT*(dsxxdFxxeg.*Nex+dsxxdFxyeg.*Ney)...
               +NweyT*(dsxxdFxyeg.*Nex+dsxydFxyeg.*Ney);
  Kuu(ne1,ne2)=+NwexT*(dsxxdFyxeg.*Nex+dsxxdFyyeg.*Ney)...
               +NweyT*(dsxydFyxeg.*Nex+dsxydFyyeg.*Ney);
  Kuu(ne2,ne1)=+NwexT*(dsxxdFyxeg.*Nex+dsxydFyxeg.*Ney)...
               +NweyT*(dsxxdFyyeg.*Nex+dsxydFyyeg.*Ney);
  Kuu(ne2,ne2)=+NwexT*(dsyxdFyxeg.*Nex+dsyxdFyyeg.*Ney)...
               +NweyT*(dsyxdFyyeg.*Nex+dsyydFyyeg.*Ney);
elseif nsd==3
  Kuu(ne1,ne1)=+NwexT*(dsxxdFxxeg.*Nex+dsxxdFxyeg.*Ney+dsxxdFxzeg.*Nez)...
               +NweyT*(dsxxdFxyeg.*Nex+dsxydFxyeg.*Ney+dsxydFxzeg.*Nez)...
               +NwezT*(dsxxdFxzeg.*Nex+dsxydFxzeg.*Ney+dsxzdFxzeg.*Nez);
  Kuu(ne1,ne2)=+NwexT*(dsxxdFyxeg.*Nex+dsxxdFyyeg.*Ney+dsxxdFyzeg.*Nez)...
               +NweyT*(dsxydFyxeg.*Nex+dsxydFyyeg.*Ney+dsxydFyzeg.*Nez)...
               +NwezT*(dsxzdFyxeg.*Nex+dsxzdFyyeg.*Ney+dsxzdFyzeg.*Nez);
  Kuu(ne1,ne3)=+NwexT*(dsxxdFzxeg.*Nex+dsxxdFzyeg.*Ney+dsxxdFzzeg.*Nez)...
               +NweyT*(dsxydFzxeg.*Nex+dsxydFzyeg.*Ney+dsxydFzzeg.*Nez)...
               +NwezT*(dsxzdFzxeg.*Nex+dsxzdFzyeg.*Ney+dsxzdFzzeg.*Nez);
  Kuu(ne2,ne1)=+NwexT*(dsxxdFyxeg.*Nex+dsxydFyxeg.*Ney+dsxzdFyxeg.*Nez)...
               +NweyT*(dsxxdFyyeg.*Nex+dsxydFyyeg.*Ney+dsxzdFyyeg.*Nez)...
               +NwezT*(dsxxdFyzeg.*Nex+dsxydFyzeg.*Ney+dsxzdFyzeg.*Nez);
  Kuu(ne2,ne2)=+NwexT*(dsyxdFyxeg.*Nex+dsyxdFyyeg.*Ney+dsyxdFyzeg.*Nez)...
               +NweyT*(dsyxdFyyeg.*Nex+dsyydFyyeg.*Ney+dsyydFyzeg.*Nez)...
               +NwezT*(dsyxdFyzeg.*Nex+dsyydFyzeg.*Ney+dsyzdFyzeg.*Nez);
  Kuu(ne2,ne3)=+NwexT*(dsyxdFzxeg.*Nex+dsyxdFzyeg.*Ney+dsyxdFzzeg.*Nez)...
               +NweyT*(dsyydFzxeg.*Nex+dsyydFzyeg.*Ney+dsyydFzzeg.*Nez)...
               +NwezT*(dsyzdFzxeg.*Nex+dsyzdFzyeg.*Ney+dsyzdFzzeg.*Nez);
  Kuu(ne3,ne1)=+NwexT*(dsxxdFzxeg.*Nex+dsxydFzxeg.*Ney+dsxzdFzxeg.*Nez)...
               +NweyT*(dsxxdFzyeg.*Nex+dsxydFzyeg.*Ney+dsxzdFzyeg.*Nez)...
               +NwezT*(dsxxdFzzeg.*Nex+dsxydFzzeg.*Ney+dsxzdFzzeg.*Nez);
  Kuu(ne3,ne2)=+NwexT*(dsyxdFzxeg.*Nex+dsyydFzxeg.*Ney+dsyzdFzxeg.*Nez)...
               +NweyT*(dsyxdFzyeg.*Nex+dsyydFzyeg.*Ney+dsyzdFzyeg.*Nez)...
               +NwezT*(dsyxdFzzeg.*Nex+dsyydFzzeg.*Ney+dsyzdFzzeg.*Nez);
  Kuu(ne3,ne3)=+NwexT*(dszxdFzxeg.*Nex+dszxdFzyeg.*Ney+dszxdFzzeg.*Nez)...
               +NweyT*(dszxdFzyeg.*Nex+dszydFzyeg.*Ney+dszydFzzeg.*Nez)...
               +NwezT*(dszxdFzzeg.*Nex+dszydFzzeg.*Ney+dszzdFzzeg.*Nez);
end

if isTimeDependent
  Kuu(ne1,ne1)=Kuu(ne1,ne1)+beta(1)/dt^2*rho*Me;
  Kuu(ne2,ne2)=Kuu(ne2,ne2)+beta(1)/dt^2*rho*Me;
  if nsd==3
    Kuu(ne3,ne3)=Kuu(ne3,ne3)+beta(1)/dt^2*rho*Me;
  end
end

% Compute rhs
if nsd==2
  fu(ne1,1)=-NwexT*(sxxeg)...
            -NweyT*(sxyeg)...
            +NweT*(fxeg);
  fu(ne2,1)=-NwexT*(syxeg)...
            -NweyT*(syyeg)...
            +NweT*(fyeg);
elseif nsd==3
  fu(ne1,1)=-NwexT*(sxxeg)...
            -NweyT*(sxyeg)...
            -NwezT*(sxzeg)...
            +NweT*(fxeg);
  fu(ne2,1)=-NwexT*(syxeg)...
            -NweyT*(syyeg)...
            -NwezT*(syzeg)...
            +NweT*(fyeg);
  fu(ne3,1)=-NwexT*(szxeg)...
            -NweyT*(szyeg)...
            -NwezT*(szzeg)...
            +NweT*(fzeg);
end

if isTimeDependent
  fu(ne1,1)=fu(ne1,1)-NweT*(rho/dt^2*uxeg*beta(1)...
                           +rho/dt^2*uoldxeg*beta(2:BDFo+2,1));
  fu(ne2,1)=fu(ne2,1)-NweT*(rho/dt^2*uyeg*beta(1)...
                           +rho/dt^2*uoldyeg*beta(2:BDFo+2,1));
  if nsd==3
    fu(ne3,1)=fu(ne3,1)-NweT*(rho/dt^2*uzeg*beta(1)...
                             +rho/dt^2*uoldzeg*beta(2:BDFo+2,1));
  end
end

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=Faces.Exterior(iFace);
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement(iD1,iD1).FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg,~,pinvNf]=mapShapeFunctions(0,RefElement(iD1,iD1),RefElement(iD1,iD1),Xf,nsd);
    
    % Compute characteristic element size
    h=sum(wfg);
    
    % Check boundary
    Boundary=Faces.Boundary(iFace);
    isDirichlet=Faces.Dirichlet(iFace);
    isNeumann=Faces.Neumann(iFace);
    if matchField(Faces,'Interface')
      isInterface=Faces.Interface(1,iFace);
    else
      isInterface=false;
    end
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    
    % Compute derivatives of shape functions
    if isDirichlet || isInterface
      Nxe=pinvNe*Nex;
      Nye=pinvNe*Ney;
      if nsd==3
        Nze=pinvNe*Nez;
      end
      Nxf=Nf*Nxe(nf1,:);
      Nyf=Nf*Nye(nf1,:);
      if nsd==3
        Nzf=Nf*Nze(nf1,:);
      end
    end
    
    % Compute variables at nodes
    uxf=uxe(nf1);
    uyf=uye(nf1);
    if nsd==3
      uzf=uze(nf1);
    end
    if isFSI && isStructure && isInterface && isTimeDependent
      uoldxf=uoldxe(nf1,:);
      uoldyf=uoldye(nf1,:);
      if nsd==3
        uoldzf=uoldze(nf1,:);
      end
    end
    if nsd==2
      Fxxf=Fxxe(nf1); Fxyf=Fxye(nf1);
      Fyxf=Fyxe(nf1); Fyyf=Fyye(nf1);
    elseif nsd==3
      Fxxf=Fxxe(nf1); Fxyf=Fxye(nf1); Fxzf=Fxze(nf1);
      Fyxf=Fyxe(nf1); Fyyf=Fyye(nf1); Fyzf=Fyze(nf1);
      Fzxf=Fzxe(nf1); Fzyf=Fzye(nf1); Fzzf=Fzze(nf1);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    uxfg=Nf*uxf;
    uyfg=Nf*uyf;
    if nsd==3
      uzfg=Nf*uzf;
    end
    if isFSI && isStructure && isInterface && isTimeDependent
      uoldxfg=Nf*uoldxf;
      uoldyfg=Nf*uoldyf;
      if nsd==3
        uoldzfg=Nf*uoldzf;
      end
      vxfg=1/dt*uxfg*alpha(1)+1/dt*uoldxfg(:,1:BDFo)*alpha(2:BDFo+1,1);
      vyfg=1/dt*uyfg*alpha(1)+1/dt*uoldyfg(:,1:BDFo)*alpha(2:BDFo+1,1);
      if nsd==3
        vzfg=1/dt*uzfg*alpha(1)+1/dt*uoldzfg(:,1:BDFo)*alpha(2:BDFo+1,1);
      end
    end
    if isDirichlet || isInterface
      if nsd==2
        Fxxfg=Nf*Fxxf; Fxyfg=Nf*Fxyf;
        Fyxfg=Nf*Fyxf; Fyyfg=Nf*Fyyf;
      elseif nsd==3
        Fxxfg=Nf*Fxxf; Fxyfg=Nf*Fxyf; Fxzfg=Nf*Fxzf;
        Fyxfg=Nf*Fyxf; Fyyfg=Nf*Fyyf; Fyzfg=Nf*Fyzf;
        Fzxfg=Nf*Fzxf; Fzyfg=Nf*Fzyf; Fzzfg=Nf*Fzzf;
      end
    end
    if isDirichlet
      uDfg=uD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      uDxfg=uDfg(:,1);
      uDyfg=uDfg(:,2);
      if nsd==3
        uDzfg=uDfg(:,3);
      end
    elseif isNeumann
      tNfg=tN(Xfg(:,1),Xfg(:,2),Xfg(:,3),t,Boundary,nx,ny,nz);
      tNxfg=tNfg(:,1);
      tNyfg=tNfg(:,2);
      if nsd==3
        tNzfg=tNfg(:,3);
      end
    end
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    Mf=NwfT*Nf;
    Mf=(Mf+Mf')/2;
    if isDirichlet || isInterface
      NwxfT=(wfg.*Nxf)';
      NwyfT=(wfg.*Nyf)';
      if nsd==3
        NwzfT=(wfg.*Nzf)';
      end
    end
    
    % Compute stress and its linearization
    if isDirichlet || isInterface
      if nsd==2
        [sfg,dsdFfg]=...
          computeStress('yes','yes','no',Parameters(iD1),Xem,...
          [Fxxfg,Fxyfg,Fyxfg,Fyyfg],Sizes(iD1));
        sxxfg=sfg(:,1); sxyfg=sfg(:,2);
        syxfg=sfg(:,3); syyfg=sfg(:,4);
        dsxxdFxxfg=dsdFfg(:,1);  dsxxdFxyfg=dsdFfg(:,2);  dsxxdFyxfg=dsdFfg(:,3);  dsxxdFyyfg=dsdFfg(:,4);
                                 dsxydFxyfg=dsdFfg(:,6);  dsxydFyxfg=dsdFfg(:,7);  dsxydFyyfg=dsdFfg(:,8);
                                                          dsyxdFyxfg=dsdFfg(:,11); dsyxdFyyfg=dsdFfg(:,12);
                                                                                   dsyydFyyfg=dsdFfg(:,16);
      elseif nsd==3
        [sfg,dsdFfg]=...
          computeStress('yes','yes','no',Parameters(iD1),Xem,...
          [Fxxfg,Fxyfg,Fxzfg,Fyxfg,Fyyfg,Fyzfg,Fzxfg,Fzyfg,Fzzfg],Sizes(iD1));
        sxxfg=sfg(:,1); sxyfg=sfg(:,2); sxzfg=sfg(:,3);
        syxfg=sfg(:,4); syyfg=sfg(:,5); syzfg=sfg(:,6);
        szxfg=sfg(:,7); szyfg=sfg(:,8); szzfg=sfg(:,9);
        dsxxdFxxfg=dsdFfg(:,1);  dsxxdFxyfg=dsdFfg(:,2);  dsxxdFxzfg=dsdFfg(:,3);  dsxxdFyxfg=dsdFfg(:,4);  dsxxdFyyfg=dsdFfg(:,5);  dsxxdFyzfg=dsdFfg(:,6);  dsxxdFzxfg=dsdFfg(:,7);  dsxxdFzyfg=dsdFfg(:,8);  dsxxdFzzfg=dsdFfg(:,9);
                                 dsxydFxyfg=dsdFfg(:,11); dsxydFxzfg=dsdFfg(:,12); dsxydFyxfg=dsdFfg(:,13); dsxydFyyfg=dsdFfg(:,14); dsxydFyzfg=dsdFfg(:,15); dsxydFzxfg=dsdFfg(:,16); dsxydFzyfg=dsdFfg(:,17); dsxydFzzfg=dsdFfg(:,18);
                                                          dsxzdFxzfg=dsdFfg(:,21); dsxzdFyxfg=dsdFfg(:,22); dsxzdFyyfg=dsdFfg(:,23); dsxzdFyzfg=dsdFfg(:,24); dsxzdFzxfg=dsdFfg(:,25); dsxzdFzyfg=dsdFfg(:,26); dsxzdFzzfg=dsdFfg(:,27);
                                                                                   dsyxdFyxfg=dsdFfg(:,31); dsyxdFyyfg=dsdFfg(:,32); dsyxdFyzfg=dsdFfg(:,33); dsyxdFzxfg=dsdFfg(:,34); dsyxdFzyfg=dsdFfg(:,35); dsyxdFzzfg=dsdFfg(:,36);
                                                                                                            dsyydFyyfg=dsdFfg(:,41); dsyydFyzfg=dsdFfg(:,42); dsyydFzxfg=dsdFfg(:,43); dsyydFzyfg=dsdFfg(:,44); dsyydFzzfg=dsdFfg(:,45);
                                                                                                                                     dsyzdFyzfg=dsdFfg(:,51); dsyzdFzxfg=dsdFfg(:,52); dsyzdFzyfg=dsdFfg(:,53); dsyzdFzzfg=dsdFfg(:,54);
                                                                                                                                                              dszxdFzxfg=dsdFfg(:,61); dszxdFzyfg=dsdFfg(:,62); dszxdFzzfg=dsdFfg(:,63);
                                                                                                                                                                                       dszydFzyfg=dsdFfg(:,71); dszydFzzfg=dsdFfg(:,72);
                                                                                                                                                                                                                dszzdFzzfg=dsdFfg(:,81);
      end
    end
    
    % Get quantities for coupling ------------------------------------------------------------------
    if isInterface
      % Get general parameters
      if isStructural
        iD2=setdiff(1:2,iD1);
      elseif isFSI
        if isMesh
          iD2=find(contains({Parameters.Problem},'Structural'));
        elseif isStructure
          iD2=find(contains({Parameters.Problem},'Fluid'));
          isFluidDM=strcmp(Parameters(iD2).Formulation,'WeaklyCompressibleFlowDM_HDG');
          isFluidVP=strcmp(Parameters(iD2).Formulation,'WeaklyCompressibleFlowVP_HDG');
          isFluidFCFV=strcmp(Parameters(iD2).Formulation,'IncompressibleFlow_FCFV');
        end
      end
      iFace2=Faces.Interface(2,iFace);
      NumElementNodes2=Sizes(iD2).NumElementNodes;
      NumFaceNodes2=Sizes(iD2).NumFaceNodes;
      
      % Compute weights at Gauss points
      FaceNodes2=RefElement(iD2,iD2).FaceNodesElem;
      [N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions(0,RefElement(iD1,iD2),...
                                                      RefElement(iD1,iD2),Xf,nsd);
      N21f=RefElement(iD2,iD1).ShapeFunctionsFace;
      
      % Compute derivatives of shape functions
      N12xf=N12f*Nxe(nf1,:);
      N12yf=N12f*Nye(nf1,:);
      if nsd==3
        N12zf=N12f*Nze(nf1,:);
      end
      
      if isStructural
        % Get solution
        U2e=SolutionGlobalCoupled{iFace};
        
        % Indices
        n2ef1=(iFace2-1)*nsd*NumFaceNodes2+(1:NumFaceNodes2);
        n2ef2=n2ef1+NumFaceNodes2;
        n2ef3=n2ef2+NumFaceNodes2;
        
        % Flip face
        Node2Match1stNode1=Faces.Interface(3,iFace);
        order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
        n2ef1=n2ef1(order);
        n2ef2=n2ef2(order);
        n2ef3=n2ef3(order);
        
        % Compute variables at nodes
        u2xf=U2e(n2ef1);
        u2yf=U2e(n2ef2);
        if nsd==3
          u2zf=U2e(n2ef3);
        end
      elseif isFSI && isMesh
        % Get solution
        u2e=reshape(SolutionGlobalCoupled{iFace},[],1);
        
        % Indices
        n2e1=1:NumElementNodes2;
        n2e2=n2e1+NumElementNodes2;
        n2e3=n2e2+NumElementNodes2;
        
        % Compute variables at nodes
        u2xe=u2e(n2e1);
        u2ye=u2e(n2e2);
        if nsd==3
          u2ze=u2e(n2e3);
        end
        
        % Indices
        n2f1=FaceNodes2(iFace2,:);
        
        % Flip face
        Node2Match1stNode1=Faces.Interface(3,iFace);
        order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
        n2f1=n2f1(order);
        
        % Compute variables at nodes
        u2xf=u2xe(n2f1);
        u2yf=u2ye(n2f1);
        if nsd==3
          u2zf=u2ze(n2f1);
        end
      elseif isFSI && isStructure
        % Get solution
        U2e=SolutionGlobalCoupled{iFace};
        
        % Indices
        n2ef1=(iFace2-1)*(nsd+1)*NumFaceNodes2+(1:NumFaceNodes2);
        n2ef2=n2ef1+NumFaceNodes2;
        n2ef3=n2ef2+NumFaceNodes2;
        n2ef4=n2ef3+NumFaceNodes2;
        
        % Flip face
        Node2Match1stNode1=Faces.Interface(3,iFace);
        order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
        n2ef1=n2ef1(order);
        n2ef2=n2ef2(order);
        n2ef3=n2ef3(order);
        n2ef4=n2ef4(order);
        
        % Compute variables at nodes
        if isFluidVP || isFluidFCFV
          v2xf=U2e(n2ef1);
          v2yf=U2e(n2ef2);
          if nsd==3
            v2zf=U2e(n2ef3);
          end
        elseif isFluidDM
          R2f=U2e(n2ef1);
          W2xf=U2e(n2ef2);
          W2yf=U2e(n2ef3);
          if nsd==3
            W2zf=U2e(n2ef4);
          end
          v2xf=W2xf./R2f;
          v2yf=W2yf./R2f;
          if nsd==3
            v2zf=W2zf./R2f;
          end
        end
      end
      
      % Compute variables at Gauss points
      if isStructural || (isFSI && isMesh)
        u2xfg=N21f*u2xf;
        u2yfg=N21f*u2yf;
        if nsd==3
          u2zfg=N21f*u2zf;
        end
      end
      if isFSI && isStructure
        v2xfg=N21f*v2xf;
        v2yfg=N21f*v2yf;
        if nsd==3
          v2zfg=N21f*v2zf;
        end
      end
      if nsd==2
        dsxxdFxx12fg=N12f*(pinvNf*dsxxdFxxfg);
        dsxxdFxy12fg=N12f*(pinvNf*dsxxdFxyfg);
        dsxxdFyx12fg=N12f*(pinvNf*dsxxdFyxfg);
        dsxxdFyy12fg=N12f*(pinvNf*dsxxdFyyfg);
        dsxydFxy12fg=N12f*(pinvNf*dsxydFxyfg);
        dsxydFyx12fg=N12f*(pinvNf*dsxydFyxfg);
        dsxydFyy12fg=N12f*(pinvNf*dsxydFyyfg);
        dsyxdFyx12fg=N12f*(pinvNf*dsyxdFyxfg);
        dsyxdFyy12fg=N12f*(pinvNf*dsyxdFyyfg);
        dsyydFyy12fg=N12f*(pinvNf*dsyydFyyfg);
      elseif nsd==3
        dsxxdFxx12fg=N12f*(pinvNf*dsxxdFxxfg);
        dsxxdFxy12fg=N12f*(pinvNf*dsxxdFxyfg);
        dsxxdFxz12fg=N12f*(pinvNf*dsxxdFxzfg);
        dsxxdFyx12fg=N12f*(pinvNf*dsxxdFyxfg);
        dsxxdFyy12fg=N12f*(pinvNf*dsxxdFyyfg);
        dsxxdFyz12fg=N12f*(pinvNf*dsxxdFyzfg);
        dsxxdFzx12fg=N12f*(pinvNf*dsxxdFzxfg);
        dsxxdFzy12fg=N12f*(pinvNf*dsxxdFzyfg);
        dsxxdFzz12fg=N12f*(pinvNf*dsxxdFzzfg);
        dsxydFxy12fg=N12f*(pinvNf*dsxydFxyfg);
        dsxydFxz12fg=N12f*(pinvNf*dsxydFxzfg);
        dsxydFyx12fg=N12f*(pinvNf*dsxydFyxfg);
        dsxydFyy12fg=N12f*(pinvNf*dsxydFyyfg);
        dsxydFyz12fg=N12f*(pinvNf*dsxydFyzfg);
        dsxydFzx12fg=N12f*(pinvNf*dsxydFzxfg);
        dsxydFzy12fg=N12f*(pinvNf*dsxydFzyfg);
        dsxydFzz12fg=N12f*(pinvNf*dsxydFzzfg);
        dsxzdFxz12fg=N12f*(pinvNf*dsxzdFxzfg);
        dsxzdFyx12fg=N12f*(pinvNf*dsxzdFyxfg);
        dsxzdFyy12fg=N12f*(pinvNf*dsxzdFyyfg);
        dsxzdFyz12fg=N12f*(pinvNf*dsxzdFyzfg);
        dsxzdFzx12fg=N12f*(pinvNf*dsxzdFzxfg);
        dsxzdFzy12fg=N12f*(pinvNf*dsxzdFzyfg);
        dsxzdFzz12fg=N12f*(pinvNf*dsxzdFzzfg);
        dsyxdFyx12fg=N12f*(pinvNf*dsyxdFyxfg);
        dsyxdFyy12fg=N12f*(pinvNf*dsyxdFyyfg);
        dsyxdFyz12fg=N12f*(pinvNf*dsyxdFyzfg);
        dsyxdFzx12fg=N12f*(pinvNf*dsyxdFzxfg);
        dsyxdFzy12fg=N12f*(pinvNf*dsyxdFzyfg);
        dsyxdFzz12fg=N12f*(pinvNf*dsyxdFzzfg);
        dsyydFyy12fg=N12f*(pinvNf*dsyydFyyfg);
        dsyydFyz12fg=N12f*(pinvNf*dsyydFyzfg);
        dsyydFzx12fg=N12f*(pinvNf*dsyydFzxfg);
        dsyydFzy12fg=N12f*(pinvNf*dsyydFzyfg);
        dsyydFzz12fg=N12f*(pinvNf*dsyydFzzfg);
        dsyzdFyz12fg=N12f*(pinvNf*dsyzdFyzfg);
        dsyzdFzx12fg=N12f*(pinvNf*dsyzdFzxfg);
        dsyzdFzy12fg=N12f*(pinvNf*dsyzdFzyfg);
        dsyzdFzz12fg=N12f*(pinvNf*dsyzdFzzfg);
        dszxdFzx12fg=N12f*(pinvNf*dszxdFzxfg);
        dszxdFzy12fg=N12f*(pinvNf*dszxdFzyfg);
        dszxdFzz12fg=N12f*(pinvNf*dszxdFzzfg);
        dszydFzy12fg=N12f*(pinvNf*dszydFzyfg);
        dszydFzz12fg=N12f*(pinvNf*dszydFzzfg);
        dszzdFzz12fg=N12f*(pinvNf*dszzdFzzfg);
      end
      
      % Compute basic matrices
      Nw12fT=(w12fg.*N12f)';
      Nw12xfT=(w12fg.*N12xf)';
      Nw12yfT=(w12fg.*N12yf)';
      if nsd==3
        Nw12zfT=(w12fg.*N12zf)';
      end
    end
    % ----------------------------------------------------------------------------------------------
    
    % Compute lhs
    if isDirichlet || isInterface
      if nsd==2
        Kuu(nf1,ne1)=Kuu(nf1,ne1)+NwfT*((-dsxxdFxxfg.*nx-dsxxdFxyfg.*ny).*Nxf...
                                       +(-dsxxdFxyfg.*nx-dsxydFxyfg.*ny).*Nyf);
        Kuu(nf1,ne2)=Kuu(nf1,ne2)+NwfT*((-dsxxdFyxfg.*nx-dsxydFyxfg.*ny).*Nxf...
                                       +(-dsxxdFyyfg.*nx-dsxydFyyfg.*ny).*Nyf);
        Kuu(nf2,ne1)=Kuu(nf2,ne1)+NwfT*((-dsxxdFyxfg.*nx-dsxxdFyyfg.*ny).*Nxf...
                                       +(-dsxydFyxfg.*nx-dsxydFyyfg.*ny).*Nyf);
        Kuu(nf2,ne2)=Kuu(nf2,ne2)+NwfT*((-dsyxdFyxfg.*nx-dsyxdFyyfg.*ny).*Nxf...
                                       +(-dsyxdFyyfg.*nx-dsyydFyyfg.*ny).*Nyf);
      elseif nsd==3
        Kuu(nf1,ne1)=Kuu(nf1,ne1)+NwfT*((-dsxxdFxxfg.*nx-dsxxdFxyfg.*ny-dsxxdFxzfg.*nz).*Nxf...
                                       +(-dsxxdFxyfg.*nx-dsxydFxyfg.*ny-dsxydFxzfg.*nz).*Nyf...
                                       +(-dsxxdFxzfg.*nx-dsxydFxzfg.*ny-dsxzdFxzfg.*nz).*Nzf);
        Kuu(nf1,ne2)=Kuu(nf1,ne2)+NwfT*((-dsxxdFyxfg.*nx-dsxydFyxfg.*ny-dsxzdFyxfg.*nz).*Nxf...
                                       +(-dsxxdFyyfg.*nx-dsxydFyyfg.*ny-dsxzdFyyfg.*nz).*Nyf...
                                       +(-dsxxdFyzfg.*nx-dsxydFyzfg.*ny-dsxzdFyzfg.*nz).*Nzf);
        Kuu(nf1,ne3)=Kuu(nf1,ne3)+NwfT*((-dsxxdFzxfg.*nx-dsxydFzxfg.*ny-dsxzdFzxfg.*nz).*Nxf...
                                       +(-dsxxdFzyfg.*nx-dsxydFzyfg.*ny-dsxzdFzyfg.*nz).*Nyf...
                                       +(-dsxxdFzzfg.*nx-dsxydFzzfg.*ny-dsxzdFzzfg.*nz).*Nzf);
        Kuu(nf2,ne1)=Kuu(nf2,ne1)+NwfT*((-dsxxdFyxfg.*nx-dsxxdFyyfg.*ny-dsxxdFyzfg.*nz).*Nxf...
                                       +(-dsxydFyxfg.*nx-dsxydFyyfg.*ny-dsxydFyzfg.*nz).*Nyf...
                                       +(-dsxzdFyxfg.*nx-dsxzdFyyfg.*ny-dsxzdFyzfg.*nz).*Nzf);
        Kuu(nf2,ne2)=Kuu(nf2,ne2)+NwfT*((-dsyxdFyxfg.*nx-dsyxdFyyfg.*ny-dsyxdFyzfg.*nz).*Nxf...
                                       +(-dsyxdFyyfg.*nx-dsyydFyyfg.*ny-dsyydFyzfg.*nz).*Nyf...
                                       +(-dsyxdFyzfg.*nx-dsyydFyzfg.*ny-dsyzdFyzfg.*nz).*Nzf);
        Kuu(nf2,ne3)=Kuu(nf2,ne3)+NwfT*((-dsyxdFzxfg.*nx-dsyydFzxfg.*ny-dsyzdFzxfg.*nz).*Nxf...
                                       +(-dsyxdFzyfg.*nx-dsyydFzyfg.*ny-dsyzdFzyfg.*nz).*Nyf...
                                       +(-dsyxdFzzfg.*nx-dsyydFzzfg.*ny-dsyzdFzzfg.*nz).*Nzf);
        Kuu(nf3,ne1)=Kuu(nf3,ne1)+NwfT*((-dsxxdFzxfg.*nx-dsxxdFzyfg.*ny-dsxxdFzzfg.*nz).*Nxf...
                                       +(-dsxydFzxfg.*nx-dsxydFzyfg.*ny-dsxydFzzfg.*nz).*Nyf...
                                       +(-dsxzdFzxfg.*nx-dsxzdFzyfg.*ny-dsxzdFzzfg.*nz).*Nzf);
        Kuu(nf3,ne2)=Kuu(nf3,ne2)+NwfT*((-dsyxdFzxfg.*nx-dsyxdFzyfg.*ny-dsyxdFzzfg.*nz).*Nxf...
                                       +(-dsyydFzxfg.*nx-dsyydFzyfg.*ny-dsyydFzzfg.*nz).*Nyf...
                                       +(-dsyzdFzxfg.*nx-dsyzdFzyfg.*ny-dsyzdFzzfg.*nz).*Nzf);
        Kuu(nf3,ne3)=Kuu(nf3,ne3)+NwfT*((-dszxdFzxfg.*nx-dszxdFzyfg.*ny-dszxdFzzfg.*nz).*Nxf...
                                       +(-dszxdFzyfg.*nx-dszydFzyfg.*ny-dszydFzzfg.*nz).*Nyf...
                                       +(-dszxdFzzfg.*nx-dszydFzzfg.*ny-dszzdFzzfg.*nz).*Nzf);
      end
    end
    
    if isDirichlet || (isStructural && isInterface) || (isFSI && isMesh && isInterface)
      Kuu(nf1,nf1)=Kuu(nf1,nf1)+gamma/h*Mf;
      Kuu(nf2,nf2)=Kuu(nf2,nf2)+gamma/h*Mf;
      if nsd==3
        Kuu(nf3,nf3)=Kuu(nf3,nf3)+gamma/h*Mf;
      end
    end
    
    if isFSI && isStructure && isInterface && isTimeDependent
      Kuu(nf1,nf1)=Kuu(nf1,nf1)+gamma/h*alpha(1)/dt*Mf;
      Kuu(nf2,nf2)=Kuu(nf2,nf2)+gamma/h*alpha(1)/dt*Mf;
      if nsd==3
        Kuu(nf3,nf3)=Kuu(nf3,nf3)+gamma/h*alpha(1)/dt*Mf;
      end
    end
    
    if isDirichlet || (isStructural && isInterface) || (isFSI && isMesh && isInterface)
      if nsd==2
        Kuu(ne1,nf1)=Kuu(ne1,nf1)+NwxfT*((-dsxxdFxxfg.*nx-dsxxdFxyfg.*ny).*Nf)...
                                 +NwyfT*((-dsxxdFxyfg.*nx-dsxydFxyfg.*ny).*Nf);
        Kuu(ne1,nf2)=Kuu(ne1,nf2)+NwxfT*((-dsxxdFyxfg.*nx-dsxxdFyyfg.*ny).*Nf)...
                                 +NwyfT*((-dsxydFyxfg.*nx-dsxydFyyfg.*ny).*Nf);
        Kuu(ne2,nf1)=Kuu(ne2,nf1)+NwxfT*((-dsxxdFyxfg.*nx-dsxydFyxfg.*ny).*Nf)...
                                 +NwyfT*((-dsxxdFyyfg.*nx-dsxydFyyfg.*ny).*Nf);
        Kuu(ne2,nf2)=Kuu(ne2,nf2)+NwxfT*((-dsyxdFyxfg.*nx-dsyxdFyyfg.*ny).*Nf)...
                                 +NwyfT*((-dsyxdFyyfg.*nx-dsyydFyyfg.*ny).*Nf);
      elseif nsd==3
        Kuu(ne1,nf1)=Kuu(ne1,nf1)+NwxfT*((-dsxxdFxxfg.*nx-dsxxdFxyfg.*ny-dsxxdFxzfg.*nz).*Nf)...
                                 +NwyfT*((-dsxxdFxyfg.*nx-dsxydFxyfg.*ny-dsxydFxzfg.*nz).*Nf)...
                                 +NwzfT*((-dsxxdFxzfg.*nx-dsxydFxzfg.*ny-dsxzdFxzfg.*nz).*Nf);
        Kuu(ne1,nf2)=Kuu(ne1,nf2)+NwxfT*((-dsxxdFyxfg.*nx-dsxxdFyyfg.*ny-dsxxdFyzfg.*nz).*Nf)...
                                 +NwyfT*((-dsxydFyxfg.*nx-dsxydFyyfg.*ny-dsxydFyzfg.*nz).*Nf)...
                                 +NwzfT*((-dsxzdFyxfg.*nx-dsxzdFyyfg.*ny-dsxzdFyzfg.*nz).*Nf);
        Kuu(ne1,nf3)=Kuu(ne1,nf3)+NwxfT*((-dsxxdFzxfg.*nx-dsxxdFzyfg.*ny-dsxxdFzzfg.*nz).*Nf)...
                                 +NwyfT*((-dsxydFzxfg.*nx-dsxydFzyfg.*ny-dsxydFzzfg.*nz).*Nf)...
                                 +NwzfT*((-dsxzdFzxfg.*nx-dsxzdFzyfg.*ny-dsxzdFzzfg.*nz).*Nf);
        Kuu(ne2,nf1)=Kuu(ne2,nf1)+NwxfT*((-dsxxdFyxfg.*nx-dsxydFyxfg.*ny-dsxzdFyxfg.*nz).*Nf)...
                                 +NwyfT*((-dsxxdFyyfg.*nx-dsxydFyyfg.*ny-dsxzdFyyfg.*nz).*Nf)...
                                 +NwzfT*((-dsxxdFyzfg.*nx-dsxydFyzfg.*ny-dsxzdFyzfg.*nz).*Nf);
        Kuu(ne2,nf2)=Kuu(ne2,nf2)+NwxfT*((-dsyxdFyxfg.*nx-dsyxdFyyfg.*ny-dsyxdFyzfg.*nz).*Nf)...
                                 +NwyfT*((-dsyxdFyyfg.*nx-dsyydFyyfg.*ny-dsyydFyzfg.*nz).*Nf)...
                                 +NwzfT*((-dsyxdFyzfg.*nx-dsyydFyzfg.*ny-dsyzdFyzfg.*nz).*Nf);
        Kuu(ne2,nf3)=Kuu(ne2,nf3)+NwxfT*((-dsyxdFzxfg.*nx-dsyxdFzyfg.*ny-dsyxdFzzfg.*nz).*Nf)...
                                 +NwyfT*((-dsyydFzxfg.*nx-dsyydFzyfg.*ny-dsyydFzzfg.*nz).*Nf)...
                                 +NwzfT*((-dsyzdFzxfg.*nx-dsyzdFzyfg.*ny-dsyzdFzzfg.*nz).*Nf);
        Kuu(ne3,nf1)=Kuu(ne3,nf1)+NwxfT*((-dsxxdFzxfg.*nx-dsxydFzxfg.*ny-dsxzdFzxfg.*nz).*Nf)...
                                 +NwyfT*((-dsxxdFzyfg.*nx-dsxydFzyfg.*ny-dsxzdFzyfg.*nz).*Nf)...
                                 +NwzfT*((-dsxxdFzzfg.*nx-dsxydFzzfg.*ny-dsxzdFzzfg.*nz).*Nf);
        Kuu(ne3,nf2)=Kuu(ne3,nf2)+NwxfT*((-dsyxdFzxfg.*nx-dsyydFzxfg.*ny-dsyzdFzxfg.*nz).*Nf)...
                                 +NwyfT*((-dsyxdFzyfg.*nx-dsyydFzyfg.*ny-dsyzdFzyfg.*nz).*Nf)...
                                 +NwzfT*((-dsyxdFzzfg.*nx-dsyydFzzfg.*ny-dsyzdFzzfg.*nz).*Nf);
        Kuu(ne3,nf3)=Kuu(ne3,nf3)+NwxfT*((-dszxdFzxfg.*nx-dszxdFzyfg.*ny-dszxdFzzfg.*nz).*Nf)...
                                 +NwyfT*((-dszxdFzyfg.*nx-dszydFzyfg.*ny-dszydFzzfg.*nz).*Nf)...
                                 +NwzfT*((-dszxdFzzfg.*nx-dszydFzzfg.*ny-dszzdFzzfg.*nz).*Nf);
      end
    end
    
    if isFSI && isStructure && isInterface && isTimeDependent
      if nsd==2
        Kuu(ne1,nf1)=Kuu(ne1,nf1)...
                     +NwxfT*(((-dsxxdFxxfg.*nx-dsxxdFxyfg.*ny)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxxdFxyfg.*nx-dsxydFxyfg.*ny)*alpha(1)/dt).*Nf);
        Kuu(ne1,nf2)=Kuu(ne1,nf2)...
                     +NwxfT*(((-dsxxdFyxfg.*nx-dsxxdFyyfg.*ny)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxydFyxfg.*nx-dsxydFyyfg.*ny)*alpha(1)/dt).*Nf);
        Kuu(ne2,nf1)=Kuu(ne2,nf1)...
                     +NwxfT*(((-dsxxdFyxfg.*nx-dsxydFyxfg.*ny)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxxdFyyfg.*nx-dsxydFyyfg.*ny)*alpha(1)/dt).*Nf);
        Kuu(ne2,nf2)=Kuu(ne2,nf2)...
                     +NwxfT*(((-dsyxdFyxfg.*nx-dsyxdFyyfg.*ny)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsyxdFyyfg.*nx-dsyydFyyfg.*ny)*alpha(1)/dt).*Nf);
      elseif nsd==3
        Kuu(ne1,nf1)=Kuu(ne1,nf1)...
                     +NwxfT*(((-dsxxdFxxfg.*nx-dsxxdFxyfg.*ny-dsxxdFxzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxxdFxyfg.*nx-dsxydFxyfg.*ny-dsxydFxzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsxxdFxzfg.*nx-dsxydFxzfg.*ny-dsxzdFxzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne1,nf2)=Kuu(ne1,nf2)...
                     +NwxfT*(((-dsxxdFyxfg.*nx-dsxxdFyyfg.*ny-dsxxdFyzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxydFyxfg.*nx-dsxydFyyfg.*ny-dsxydFyzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsxzdFyxfg.*nx-dsxzdFyyfg.*ny-dsxzdFyzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne1,nf3)=Kuu(ne1,nf3)...
                     +NwxfT*(((-dsxxdFzxfg.*nx-dsxxdFzyfg.*ny-dsxxdFzzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxydFzxfg.*nx-dsxydFzyfg.*ny-dsxydFzzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsxzdFzxfg.*nx-dsxzdFzyfg.*ny-dsxzdFzzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne2,nf1)=Kuu(ne2,nf1)...
                     +NwxfT*(((-dsxxdFyxfg.*nx-dsxydFyxfg.*ny-dsxzdFyxfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxxdFyyfg.*nx-dsxydFyyfg.*ny-dsxzdFyyfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsxxdFyzfg.*nx-dsxydFyzfg.*ny-dsxzdFyzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne2,nf2)=Kuu(ne2,nf2)...
                     +NwxfT*(((-dsyxdFyxfg.*nx-dsyxdFyyfg.*ny-dsyxdFyzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsyxdFyyfg.*nx-dsyydFyyfg.*ny-dsyydFyzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsyxdFyzfg.*nx-dsyydFyzfg.*ny-dsyzdFyzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne2,nf3)=Kuu(ne2,nf3)...
                     +NwxfT*(((-dsyxdFzxfg.*nx-dsyxdFzyfg.*ny-dsyxdFzzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsyydFzxfg.*nx-dsyydFzyfg.*ny-dsyydFzzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsyzdFzxfg.*nx-dsyzdFzyfg.*ny-dsyzdFzzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne3,nf1)=Kuu(ne3,nf1)...
                     +NwxfT*(((-dsxxdFzxfg.*nx-dsxydFzxfg.*ny-dsxzdFzxfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsxxdFzyfg.*nx-dsxydFzyfg.*ny-dsxzdFzyfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsxxdFzzfg.*nx-dsxydFzzfg.*ny-dsxzdFzzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne3,nf2)=Kuu(ne3,nf2)...
                     +NwxfT*(((-dsyxdFzxfg.*nx-dsyydFzxfg.*ny-dsyzdFzxfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dsyxdFzyfg.*nx-dsyydFzyfg.*ny-dsyzdFzyfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dsyxdFzzfg.*nx-dsyydFzzfg.*ny-dsyzdFzzfg.*nz)*alpha(1)/dt).*Nf);
        Kuu(ne3,nf3)=Kuu(ne3,nf3)...
                     +NwxfT*(((-dszxdFzxfg.*nx-dszxdFzyfg.*ny-dszxdFzzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwyfT*(((-dszxdFzyfg.*nx-dszydFzyfg.*ny-dszydFzzfg.*nz)*alpha(1)/dt).*Nf)...
                     +NwzfT*(((-dszxdFzzfg.*nx-dszydFzzfg.*ny-dszzdFzzfg.*nz)*alpha(1)/dt).*Nf);
      end
    end
    
    % Compute rhs
    if isDirichlet
      if nsd==2
        fu(nf1,1)=fu(nf1,1)+NwfT*(+sxxfg.*nx+sxyfg.*ny...
                                  -gamma/h*(uxfg-uDxfg));
        fu(nf2,1)=fu(nf2,1)+NwfT*(+syxfg.*nx+syyfg.*ny...
                                  -gamma/h*(uyfg-uDyfg));
      elseif nsd==3
        fu(nf1,1)=fu(nf1,1)+NwfT*(+sxxfg.*nx+sxyfg.*ny+sxzfg.*nz...
                                  -gamma/h*(uxfg-uDxfg));
        fu(nf2,1)=fu(nf2,1)+NwfT*(+syxfg.*nx+syyfg.*ny+syzfg.*nz...
                                  -gamma/h*(uyfg-uDyfg));
        fu(nf3,1)=fu(nf3,1)+NwfT*(+szxfg.*nx+szyfg.*ny+szzfg.*nz...
                                  -gamma/h*(uzfg-uDzfg));
      end
      
      if nsd==2
        fu(ne1,1)=fu(ne1,1)...
                  +NwxfT*(+(+dsxxdFxxfg.*nx+dsxxdFxyfg.*ny).*(uxfg-uDxfg)...
                          +(+dsxxdFyxfg.*nx+dsxxdFyyfg.*ny).*(uyfg-uDyfg))...
                  +NwyfT*(+(+dsxxdFxyfg.*nx+dsxydFxyfg.*ny).*(uxfg-uDxfg)...
                          +(+dsxydFyxfg.*nx+dsxydFyyfg.*ny).*(uyfg-uDyfg));
        fu(ne2,1)=fu(ne2,1)...
                  +NwxfT*(+(+dsxxdFyxfg.*nx+dsxydFyxfg.*ny).*(uxfg-uDxfg)...
                          +(+dsyxdFyxfg.*nx+dsyxdFyyfg.*ny).*(uyfg-uDyfg))...
                  +NwyfT*(+(+dsxxdFyyfg.*nx+dsxydFyyfg.*ny).*(uxfg-uDxfg)...
                          +(+dsyxdFyyfg.*nx+dsyydFyyfg.*ny).*(uyfg-uDyfg));
      elseif nsd==3
        fu(ne1,1)=fu(ne1,1)...
                  +NwxfT*(+(+dsxxdFxxfg.*nx+dsxxdFxyfg.*ny+dsxxdFxzfg.*nz).*(uxfg-uDxfg)...
                          +(+dsxxdFyxfg.*nx+dsxxdFyyfg.*ny+dsxxdFyzfg.*nz).*(uyfg-uDyfg)...
                          +(+dsxxdFzxfg.*nx+dsxxdFzyfg.*ny+dsxxdFzzfg.*nz).*(uzfg-uDzfg))...
                  +NwyfT*(+(+dsxxdFxyfg.*nx+dsxydFxyfg.*ny+dsxydFxzfg.*nz).*(uxfg-uDxfg)...
                          +(+dsxydFyxfg.*nx+dsxydFyyfg.*ny+dsxydFyzfg.*nz).*(uyfg-uDyfg)...
                          +(+dsxydFzxfg.*nx+dsxydFzyfg.*ny+dsxydFzzfg.*nz).*(uzfg-uDzfg))...
                  +NwzfT*(+(+dsxxdFxzfg.*nx+dsxydFxzfg.*ny+dsxzdFxzfg.*nz).*(uxfg-uDxfg)...
                          +(+dsxzdFyxfg.*nx+dsxzdFyyfg.*ny+dsxzdFyzfg.*nz).*(uyfg-uDyfg)...
                          +(+dsxzdFzxfg.*nx+dsxzdFzyfg.*ny+dsxzdFzzfg.*nz).*(uzfg-uDzfg));
        fu(ne2,1)=fu(ne2,1)...
                  +NwxfT*(+(+dsxxdFyxfg.*nx+dsxydFyxfg.*ny+dsxzdFyxfg.*nz).*(uxfg-uDxfg)...
                          +(+dsyxdFyxfg.*nx+dsyxdFyyfg.*ny+dsyxdFyzfg.*nz).*(uyfg-uDyfg)...
                          +(+dsyxdFzxfg.*nx+dsyxdFzyfg.*ny+dsyxdFzzfg.*nz).*(uzfg-uDzfg))...
                  +NwyfT*(+(+dsxxdFyyfg.*nx+dsxydFyyfg.*ny+dsxzdFyyfg.*nz).*(uxfg-uDxfg)...
                          +(+dsyxdFyyfg.*nx+dsyydFyyfg.*ny+dsyydFyzfg.*nz).*(uyfg-uDyfg)...
                          +(+dsyydFzxfg.*nx+dsyydFzyfg.*ny+dsyydFzzfg.*nz).*(uzfg-uDzfg))...
                  +NwzfT*(+(+dsxxdFyzfg.*nx+dsxydFyzfg.*ny+dsxzdFyzfg.*nz).*(uxfg-uDxfg)...
                          +(+dsyxdFyzfg.*nx+dsyydFyzfg.*ny+dsyzdFyzfg.*nz).*(uyfg-uDyfg)...
                          +(+dsyzdFzxfg.*nx+dsyzdFzyfg.*ny+dsyzdFzzfg.*nz).*(uzfg-uDzfg));
        fu(ne3,1)=fu(ne3,1)...
                  +NwxfT*(+(+dsxxdFzxfg.*nx+dsxydFzxfg.*ny+dsxzdFzxfg.*nz).*(uxfg-uDxfg)...
                          +(+dsyxdFzxfg.*nx+dsyydFzxfg.*ny+dsyzdFzxfg.*nz).*(uyfg-uDyfg)...
                          +(+dszxdFzxfg.*nx+dszxdFzyfg.*ny+dszxdFzzfg.*nz).*(uzfg-uDzfg))...
                  +NwyfT*(+(+dsxxdFzyfg.*nx+dsxydFzyfg.*ny+dsxzdFzyfg.*nz).*(uxfg-uDxfg)...
                          +(+dsyxdFzyfg.*nx+dsyydFzyfg.*ny+dsyzdFzyfg.*nz).*(uyfg-uDyfg)...
                          +(+dszxdFzyfg.*nx+dszydFzyfg.*ny+dszydFzzfg.*nz).*(uzfg-uDzfg))...
                  +NwzfT*(+(+dsxxdFzzfg.*nx+dsxydFzzfg.*ny+dsxzdFzzfg.*nz).*(uxfg-uDxfg)...
                          +(+dsyxdFzzfg.*nx+dsyydFzzfg.*ny+dsyzdFzzfg.*nz).*(uyfg-uDyfg)...
                          +(+dszxdFzzfg.*nx+dszydFzzfg.*ny+dszzdFzzfg.*nz).*(uzfg-uDzfg));
      end
    end
    
    if (isStructural && isInterface) || (isFSI && isMesh && isInterface)
      if nsd==2
        fu(nf1,1)=fu(nf1,1)+NwfT*(+sxxfg.*nx+sxyfg.*ny...
                                  -gamma/h*uxfg);
        fu(nf2,1)=fu(nf2,1)+NwfT*(+syxfg.*nx+syyfg.*ny...
                                  -gamma/h*uyfg);
      elseif nsd==3
        fu(nf1,1)=fu(nf1,1)+NwfT*(+sxxfg.*nx+sxyfg.*ny+sxzfg.*nz...
                                  -gamma/h*uxfg);
        fu(nf2,1)=fu(nf2,1)+NwfT*(+syxfg.*nx+syyfg.*ny+syzfg.*nz...
                                  -gamma/h*uyfg);
        fu(nf3,1)=fu(nf3,1)+NwfT*(+szxfg.*nx+szyfg.*ny+szzfg.*nz...
                                  -gamma/h*uzfg);
      end
      
      fu(nf1,1)=fu(nf1,1)+Nw12fT*(gamma/h*u2xfg);
      fu(nf2,1)=fu(nf2,1)+Nw12fT*(gamma/h*u2yfg);
      if nsd==3
        fu(nf3,1)=fu(nf3,1)+Nw12fT*(gamma/h*u2zfg);
      end
      
      if nsd==2
        fu(ne1,1)=fu(ne1,1)+NwxfT*(+(+dsxxdFxxfg.*nx+dsxxdFxyfg.*ny).*uxfg...
                                   +(+dsxxdFyxfg.*nx+dsxxdFyyfg.*ny).*uyfg)...
                           +NwyfT*(+(+dsxxdFxyfg.*nx+dsxydFxyfg.*ny).*uxfg...
                                   +(+dsxydFyxfg.*nx+dsxydFyyfg.*ny).*uyfg);
        fu(ne2,1)=fu(ne2,1)+NwxfT*(+(+dsxxdFyxfg.*nx+dsxydFyxfg.*ny).*uxfg...
                                   +(+dsyxdFyxfg.*nx+dsyxdFyyfg.*ny).*uyfg)...
                           +NwyfT*(+(+dsxxdFyyfg.*nx+dsxydFyyfg.*ny).*uxfg...
                                   +(+dsyxdFyyfg.*nx+dsyydFyyfg.*ny).*uyfg);
      elseif nsd==3
        fu(ne1,1)=fu(ne1,1)+NwxfT*(+(+dsxxdFxxfg.*nx+dsxxdFxyfg.*ny+dsxxdFxzfg.*nz).*uxfg...
                                   +(+dsxxdFyxfg.*nx+dsxxdFyyfg.*ny+dsxxdFyzfg.*nz).*uyfg...
                                   +(+dsxxdFzxfg.*nx+dsxxdFzyfg.*ny+dsxxdFzzfg.*nz).*uzfg)...
                           +NwyfT*(+(+dsxxdFxyfg.*nx+dsxydFxyfg.*ny+dsxydFxzfg.*nz).*uxfg...
                                   +(+dsxydFyxfg.*nx+dsxydFyyfg.*ny+dsxydFyzfg.*nz).*uyfg...
                                   +(+dsxydFzxfg.*nx+dsxydFzyfg.*ny+dsxydFzzfg.*nz).*uzfg)...
                           +NwzfT*(+(+dsxxdFxzfg.*nx+dsxydFxzfg.*ny+dsxzdFxzfg.*nz).*uxfg...
                                   +(+dsxzdFyxfg.*nx+dsxzdFyyfg.*ny+dsxzdFyzfg.*nz).*uyfg...
                                   +(+dsxzdFzxfg.*nx+dsxzdFzyfg.*ny+dsxzdFzzfg.*nz).*uzfg);
        fu(ne2,1)=fu(ne2,1)+NwxfT*(+(+dsxxdFyxfg.*nx+dsxydFyxfg.*ny+dsxzdFyxfg.*nz).*uxfg...
                                   +(+dsyxdFyxfg.*nx+dsyxdFyyfg.*ny+dsyxdFyzfg.*nz).*uyfg...
                                   +(+dsyxdFzxfg.*nx+dsyxdFzyfg.*ny+dsyxdFzzfg.*nz).*uzfg)...
                           +NwyfT*(+(+dsxxdFyyfg.*nx+dsxydFyyfg.*ny+dsxzdFyyfg.*nz).*uxfg...
                                   +(+dsyxdFyyfg.*nx+dsyydFyyfg.*ny+dsyydFyzfg.*nz).*uyfg...
                                   +(+dsyydFzxfg.*nx+dsyydFzyfg.*ny+dsyydFzzfg.*nz).*uzfg)...
                           +NwzfT*(+(+dsxxdFyzfg.*nx+dsxydFyzfg.*ny+dsxzdFyzfg.*nz).*uxfg...
                                   +(+dsyxdFyzfg.*nx+dsyydFyzfg.*ny+dsyzdFyzfg.*nz).*uyfg...
                                   +(+dsyzdFzxfg.*nx+dsyzdFzyfg.*ny+dsyzdFzzfg.*nz).*uzfg);
        fu(ne3,1)=fu(ne3,1)+NwxfT*(+(+dsxxdFzxfg.*nx+dsxydFzxfg.*ny+dsxzdFzxfg.*nz).*uxfg...
                                   +(+dsyxdFzxfg.*nx+dsyydFzxfg.*ny+dsyzdFzxfg.*nz).*uyfg...
                                   +(+dszxdFzxfg.*nx+dszxdFzyfg.*ny+dszxdFzzfg.*nz).*uzfg)...
                           +NwyfT*(+(+dsxxdFzyfg.*nx+dsxydFzyfg.*ny+dsxzdFzyfg.*nz).*uxfg...
                                   +(+dsyxdFzyfg.*nx+dsyydFzyfg.*ny+dsyzdFzyfg.*nz).*uyfg...
                                   +(+dszxdFzyfg.*nx+dszydFzyfg.*ny+dszydFzzfg.*nz).*uzfg)...
                           +NwzfT*(+(+dsxxdFzzfg.*nx+dsxydFzzfg.*ny+dsxzdFzzfg.*nz).*uxfg...
                                   +(+dsyxdFzzfg.*nx+dsyydFzzfg.*ny+dsyzdFzzfg.*nz).*uyfg...
                                   +(+dszxdFzzfg.*nx+dszydFzzfg.*ny+dszzdFzzfg.*nz).*uzfg);
      end
      
      if nsd==2
        fu(ne1,1)=fu(ne1,1)...
                   -Nw12xfT*(+(+dsxxdFxx12fg.*n12x+dsxxdFxy12fg.*n12y).*u2xfg...
                             +(+dsxxdFyx12fg.*n12x+dsxxdFyy12fg.*n12y).*u2yfg)...
                   -Nw12yfT*(+(+dsxxdFxy12fg.*n12x+dsxydFxy12fg.*n12y).*u2xfg...
                             +(+dsxydFyx12fg.*n12x+dsxydFyy12fg.*n12y).*u2yfg);
        fu(ne2,1)=fu(ne2,1)...
                   -Nw12xfT*(+(+dsxxdFyx12fg.*n12x+dsxydFyx12fg.*n12y).*u2xfg...
                             +(+dsyxdFyx12fg.*n12x+dsyxdFyy12fg.*n12y).*u2yfg)...
                   -Nw12yfT*(+(+dsxxdFyy12fg.*n12x+dsxydFyy12fg.*n12y).*u2xfg...
                             +(+dsyxdFyy12fg.*n12x+dsyydFyy12fg.*n12y).*u2yfg);
      elseif nsd==3
        fu(ne1,1)=fu(ne1,1)...
                   -Nw12xfT*(+(+dsxxdFxx12fg.*n12x+dsxxdFxy12fg.*n12y+dsxxdFxz12fg.*n12z).*u2xfg...
                             +(+dsxxdFyx12fg.*n12x+dsxxdFyy12fg.*n12y+dsxxdFyz12fg.*n12z).*u2yfg...
                             +(+dsxxdFzx12fg.*n12x+dsxxdFzy12fg.*n12y+dsxxdFzz12fg.*n12z).*u2zfg)...
                   -Nw12yfT*(+(+dsxxdFxy12fg.*n12x+dsxydFxy12fg.*n12y+dsxydFxz12fg.*n12z).*u2xfg...
                             +(+dsxydFyx12fg.*n12x+dsxydFyy12fg.*n12y+dsxydFyz12fg.*n12z).*u2yfg...
                             +(+dsxydFzx12fg.*n12x+dsxydFzy12fg.*n12y+dsxydFzz12fg.*n12z).*u2zfg)...
                   -Nw12zfT*(+(+dsxxdFxz12fg.*n12x+dsxydFxz12fg.*n12y+dsxzdFxz12fg.*n12z).*u2xfg...
                             +(+dsxzdFyx12fg.*n12x+dsxzdFyy12fg.*n12y+dsxzdFyz12fg.*n12z).*u2yfg...
                             +(+dsxzdFzx12fg.*n12x+dsxzdFzy12fg.*n12y+dsxzdFzz12fg.*n12z).*u2zfg);
        fu(ne2,1)=fu(ne2,1)...
                   -Nw12xfT*(+(+dsxxdFyx12fg.*n12x+dsxydFyx12fg.*n12y+dsxzdFyx12fg.*n12z).*u2xfg...
                             +(+dsyxdFyx12fg.*n12x+dsyxdFyy12fg.*n12y+dsyxdFyz12fg.*n12z).*u2yfg...
                             +(+dsyxdFzx12fg.*n12x+dsyxdFzy12fg.*n12y+dsyxdFzz12fg.*n12z).*u2zfg)...
                   -Nw12yfT*(+(+dsxxdFyy12fg.*n12x+dsxydFyy12fg.*n12y+dsxzdFyy12fg.*n12z).*u2xfg...
                             +(+dsyxdFyy12fg.*n12x+dsyydFyy12fg.*n12y+dsyydFyz12fg.*n12z).*u2yfg...
                             +(+dsyydFzx12fg.*n12x+dsyydFzy12fg.*n12y+dsyydFzz12fg.*n12z).*u2zfg)...
                   -Nw12zfT*(+(+dsxxdFyz12fg.*n12x+dsxydFyz12fg.*n12y+dsxzdFyz12fg.*n12z).*u2xfg...
                             +(+dsyxdFyz12fg.*n12x+dsyydFyz12fg.*n12y+dsyzdFyz12fg.*n12z).*u2yfg...
                             +(+dsyzdFzx12fg.*n12x+dsyzdFzy12fg.*n12y+dsyzdFzz12fg.*n12z).*u2zfg);
        fu(ne3,1)=fu(ne3,1)...
                   -Nw12xfT*(+(+dsxxdFzx12fg.*n12x+dsxydFzx12fg.*n12y+dsxzdFzx12fg.*n12z).*u2xfg...
                             +(+dsyxdFzx12fg.*n12x+dsyydFzx12fg.*n12y+dsyzdFzx12fg.*n12z).*u2yfg...
                             +(+dszxdFzx12fg.*n12x+dszxdFzy12fg.*n12y+dszxdFzz12fg.*n12z).*u2zfg)...
                   -Nw12yfT*(+(+dsxxdFzy12fg.*n12x+dsxydFzy12fg.*n12y+dsxzdFzy12fg.*n12z).*u2xfg...
                             +(+dsyxdFzy12fg.*n12x+dsyydFzy12fg.*n12y+dsyzdFzy12fg.*n12z).*u2yfg...
                             +(+dszxdFzy12fg.*n12x+dszydFzy12fg.*n12y+dszydFzz12fg.*n12z).*u2zfg)...
                   -Nw12zfT*(+(+dsxxdFzz12fg.*n12x+dsxydFzz12fg.*n12y+dsxzdFzz12fg.*n12z).*u2xfg...
                             +(+dsyxdFzz12fg.*n12x+dsyydFzz12fg.*n12y+dsyzdFzz12fg.*n12z).*u2yfg...
                             +(+dszxdFzz12fg.*n12x+dszydFzz12fg.*n12y+dszzdFzz12fg.*n12z).*u2zfg);
      end
    end
    
    if isFSI && isStructure && isInterface
      if nsd==2
        fu(nf1,1)=fu(nf1,1)+NwfT*(+sxxfg.*nx+sxyfg.*ny);
        fu(nf2,1)=fu(nf2,1)+NwfT*(+syxfg.*nx+syyfg.*ny);
      elseif nsd==3
        fu(nf1,1)=fu(nf1,1)+NwfT*(+sxxfg.*nx+sxyfg.*ny+sxzfg.*nz);
        fu(nf2,1)=fu(nf2,1)+NwfT*(+syxfg.*nx+syyfg.*ny+syzfg.*nz);
        fu(nf3,1)=fu(nf3,1)+NwfT*(+szxfg.*nx+szyfg.*ny+szzfg.*nz);
      end
    end
    
    if isFSI && isStructure && isInterface && isTimeDependent
      fu(nf1,1)=fu(nf1,1)-NwfT*(gamma/h*vxfg);
      fu(nf2,1)=fu(nf2,1)-NwfT*(gamma/h*vyfg);
      if nsd==3
        fu(nf3,1)=fu(nf3,1)-NwfT*(gamma/h*vzfg);
      end
    end
    
    if isFSI && isStructure && isInterface
      fu(nf1,1)=fu(nf1,1)+Nw12fT*(gamma/h*v2xfg);
      fu(nf2,1)=fu(nf2,1)+Nw12fT*(gamma/h*v2yfg);
      if nsd==3
        fu(nf3,1)=fu(nf3,1)+Nw12fT*(gamma/h*v2zfg);
      end
    end
    
    if isFSI && isStructure && isInterface && isTimeDependent
      if nsd==2
        fu(ne1,1)=fu(ne1,1)+NwxfT*(+(+dsxxdFxxfg.*nx+dsxxdFxyfg.*ny).*vxfg...
                                   +(+dsxxdFyxfg.*nx+dsxxdFyyfg.*ny).*vyfg)...
                           +NwyfT*(+(+dsxxdFxyfg.*nx+dsxydFxyfg.*ny).*vxfg...
                                   +(+dsxydFyxfg.*nx+dsxydFyyfg.*ny).*vyfg);
        fu(ne2,1)=fu(ne2,1)+NwxfT*(+(+dsxxdFyxfg.*nx+dsxydFyxfg.*ny).*vxfg...
                                   +(+dsyxdFyxfg.*nx+dsyxdFyyfg.*ny).*vyfg)...
                           +NwyfT*(+(+dsxxdFyyfg.*nx+dsxydFyyfg.*ny).*vxfg...
                                   +(+dsyxdFyyfg.*nx+dsyydFyyfg.*ny).*vyfg);
      elseif nsd==3
        fu(ne1,1)=fu(ne1,1)+NwxfT*(+(+dsxxdFxxfg.*nx+dsxxdFxyfg.*ny+dsxxdFxzfg.*nz).*vxfg...
                                   +(+dsxxdFyxfg.*nx+dsxxdFyyfg.*ny+dsxxdFyzfg.*nz).*vyfg...
                                   +(+dsxxdFzxfg.*nx+dsxxdFzyfg.*ny+dsxxdFzzfg.*nz).*vzfg)...
                           +NwyfT*(+(+dsxxdFxyfg.*nx+dsxydFxyfg.*ny+dsxydFxzfg.*nz).*vxfg...
                                   +(+dsxydFyxfg.*nx+dsxydFyyfg.*ny+dsxydFyzfg.*nz).*vyfg...
                                   +(+dsxydFzxfg.*nx+dsxydFzyfg.*ny+dsxydFzzfg.*nz).*vzfg)...
                           +NwzfT*(+(+dsxxdFxzfg.*nx+dsxydFxzfg.*ny+dsxzdFxzfg.*nz).*vxfg...
                                   +(+dsxzdFyxfg.*nx+dsxzdFyyfg.*ny+dsxzdFyzfg.*nz).*vyfg...
                                   +(+dsxzdFzxfg.*nx+dsxzdFzyfg.*ny+dsxzdFzzfg.*nz).*vzfg);
        fu(ne2,1)=fu(ne2,1)+NwxfT*(+(+dsxxdFyxfg.*nx+dsxydFyxfg.*ny+dsxzdFyxfg.*nz).*vxfg...
                                   +(+dsyxdFyxfg.*nx+dsyxdFyyfg.*ny+dsyxdFyzfg.*nz).*vyfg...
                                   +(+dsyxdFzxfg.*nx+dsyxdFzyfg.*ny+dsyxdFzzfg.*nz).*vzfg)...
                           +NwyfT*(+(+dsxxdFyyfg.*nx+dsxydFyyfg.*ny+dsxzdFyyfg.*nz).*vxfg...
                                   +(+dsyxdFyyfg.*nx+dsyydFyyfg.*ny+dsyydFyzfg.*nz).*vyfg...
                                   +(+dsyydFzxfg.*nx+dsyydFzyfg.*ny+dsyydFzzfg.*nz).*vzfg)...
                           +NwzfT*(+(+dsxxdFyzfg.*nx+dsxydFyzfg.*ny+dsxzdFyzfg.*nz).*vxfg...
                                   +(+dsyxdFyzfg.*nx+dsyydFyzfg.*ny+dsyzdFyzfg.*nz).*vyfg...
                                   +(+dsyzdFzxfg.*nx+dsyzdFzyfg.*ny+dsyzdFzzfg.*nz).*vzfg);
        fu(ne3,1)=fu(ne3,1)+NwxfT*(+(+dsxxdFzxfg.*nx+dsxydFzxfg.*ny+dsxzdFzxfg.*nz).*vxfg...
                                   +(+dsyxdFzxfg.*nx+dsyydFzxfg.*ny+dsyzdFzxfg.*nz).*vyfg...
                                   +(+dszxdFzxfg.*nx+dszxdFzyfg.*ny+dszxdFzzfg.*nz).*vzfg)...
                           +NwyfT*(+(+dsxxdFzyfg.*nx+dsxydFzyfg.*ny+dsxzdFzyfg.*nz).*vxfg...
                                   +(+dsyxdFzyfg.*nx+dsyydFzyfg.*ny+dsyzdFzyfg.*nz).*vyfg...
                                   +(+dszxdFzyfg.*nx+dszydFzyfg.*ny+dszydFzzfg.*nz).*vzfg)...
                           +NwzfT*(+(+dsxxdFzzfg.*nx+dsxydFzzfg.*ny+dsxzdFzzfg.*nz).*vxfg...
                                   +(+dsyxdFzzfg.*nx+dsyydFzzfg.*ny+dsyzdFzzfg.*nz).*vyfg...
                                   +(+dszxdFzzfg.*nx+dszydFzzfg.*ny+dszzdFzzfg.*nz).*vzfg);
      end
    end
    
    if isFSI && isStructure && isInterface
      if nsd==2
        fu(ne1,1)=fu(ne1,1)...
                  -Nw12xfT*(+(+dsxxdFxx12fg.*n12x+dsxxdFxy12fg.*n12y).*v2xfg...
                            +(+dsxxdFyx12fg.*n12x+dsxxdFyy12fg.*n12y).*v2yfg)...
                  -Nw12yfT*(+(+dsxxdFxy12fg.*n12x+dsxydFxy12fg.*n12y).*v2xfg...
                            +(+dsxydFyx12fg.*n12x+dsxydFyy12fg.*n12y).*v2yfg);
        fu(ne2,1)=fu(ne2,1)...
                  -Nw12xfT*(+(+dsxxdFyx12fg.*n12x+dsxydFyx12fg.*n12y).*v2xfg...
                            +(+dsyxdFyx12fg.*n12x+dsyxdFyy12fg.*n12y).*v2yfg)...
                  -Nw12yfT*(+(+dsxxdFyy12fg.*n12x+dsxydFyy12fg.*n12y).*v2xfg...
                            +(+dsyxdFyy12fg.*n12x+dsyydFyy12fg.*n12y).*v2yfg);
      elseif nsd==3
        fu(ne1,1)=fu(ne1,1)...
                  -Nw12xfT*(+(+dsxxdFxx12fg.*n12x+dsxxdFxy12fg.*n12y+dsxxdFxz12fg.*n12z).*v2xfg...
                            +(+dsxxdFyx12fg.*n12x+dsxxdFyy12fg.*n12y+dsxxdFyz12fg.*n12z).*v2yfg...
                            +(+dsxxdFzx12fg.*n12x+dsxxdFzy12fg.*n12y+dsxxdFzz12fg.*n12z).*v2zfg)...
                  -Nw12yfT*(+(+dsxxdFxy12fg.*n12x+dsxydFxy12fg.*n12y+dsxydFxz12fg.*n12z).*v2xfg...
                            +(+dsxydFyx12fg.*n12x+dsxydFyy12fg.*n12y+dsxydFyz12fg.*n12z).*v2yfg...
                            +(+dsxydFzx12fg.*n12x+dsxydFzy12fg.*n12y+dsxydFzz12fg.*n12z).*v2zfg)...
                  -Nw12zfT*(+(+dsxxdFxz12fg.*n12x+dsxydFxz12fg.*n12y+dsxzdFxz12fg.*n12z).*v2xfg...
                            +(+dsxzdFyx12fg.*n12x+dsxzdFyy12fg.*n12y+dsxzdFyz12fg.*n12z).*v2yfg...
                            +(+dsxzdFzx12fg.*n12x+dsxzdFzy12fg.*n12y+dsxzdFzz12fg.*n12z).*v2zfg);
        fu(ne2,1)=fu(ne2,1)...
                  -Nw12xfT*(+(+dsxxdFyx12fg.*n12x+dsxydFyx12fg.*n12y+dsxzdFyx12fg.*n12z).*v2xfg...
                            +(+dsyxdFyx12fg.*n12x+dsyxdFyy12fg.*n12y+dsyxdFyz12fg.*n12z).*v2yfg...
                            +(+dsyxdFzx12fg.*n12x+dsyxdFzy12fg.*n12y+dsyxdFzz12fg.*n12z).*v2zfg)...
                  -Nw12yfT*(+(+dsxxdFyy12fg.*n12x+dsxydFyy12fg.*n12y+dsxzdFyy12fg.*n12z).*v2xfg...
                            +(+dsyxdFyy12fg.*n12x+dsyydFyy12fg.*n12y+dsyydFyz12fg.*n12z).*v2yfg...
                            +(+dsyydFzx12fg.*n12x+dsyydFzy12fg.*n12y+dsyydFzz12fg.*n12z).*v2zfg)...
                  -Nw12zfT*(+(+dsxxdFyz12fg.*n12x+dsxydFyz12fg.*n12y+dsxzdFyz12fg.*n12z).*v2xfg...
                            +(+dsyxdFyz12fg.*n12x+dsyydFyz12fg.*n12y+dsyzdFyz12fg.*n12z).*v2yfg...
                            +(+dsyzdFzx12fg.*n12x+dsyzdFzy12fg.*n12y+dsyzdFzz12fg.*n12z).*v2zfg);
        fu(ne3,1)=fu(ne3,1)...
                  -Nw12xfT*(+(+dsxxdFzx12fg.*n12x+dsxydFzx12fg.*n12y+dsxzdFzx12fg.*n12z).*v2xfg...
                            +(+dsyxdFzx12fg.*n12x+dsyydFzx12fg.*n12y+dsyzdFzx12fg.*n12z).*v2yfg...
                            +(+dszxdFzx12fg.*n12x+dszxdFzy12fg.*n12y+dszxdFzz12fg.*n12z).*v2zfg)...
                  -Nw12yfT*(+(+dsxxdFzy12fg.*n12x+dsxydFzy12fg.*n12y+dsxzdFzy12fg.*n12z).*v2xfg...
                            +(+dsyxdFzy12fg.*n12x+dsyydFzy12fg.*n12y+dsyzdFzy12fg.*n12z).*v2yfg...
                            +(+dszxdFzy12fg.*n12x+dszydFzy12fg.*n12y+dszydFzz12fg.*n12z).*v2zfg)...
                  -Nw12zfT*(+(+dsxxdFzz12fg.*n12x+dsxydFzz12fg.*n12y+dsxzdFzz12fg.*n12z).*v2xfg...
                            +(+dsyxdFzz12fg.*n12x+dsyydFzz12fg.*n12y+dsyzdFzz12fg.*n12z).*v2yfg...
                            +(+dszxdFzz12fg.*n12x+dszydFzz12fg.*n12y+dszzdFzz12fg.*n12z).*v2zfg);
      end
    end
    
    if isNeumann
      fu(nf1,1)=fu(nf1,1)+NwfT*(tNxfg);
      fu(nf2,1)=fu(nf2,1)+NwfT*(tNyfg);
      if nsd==3
        fu(nf3,1)=fu(nf3,1)+NwfT*(tNzfg);
      end
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Lhs for global problem
LhsGlobal=Kuu;
if not(isFSI && isStructure)
  LhsGlobal=(LhsGlobal+LhsGlobal.')/2;
end

% Rhs for global problem
RhsGlobal=fu;

end

%% Do coupling element
function [LhsCoup]=doCouplingElement(...
  iFaceInterface,iD1,iD2,Block,Simulation,Parameters,Mesh,Faces,RefElement,Sizes)       

% Get general parameters
isStructural=strcmp(Simulation.Problem,'Structural');
isFSI=strcmp(Simulation.Problem,'FluidStructureInteraction');
if isFSI
  isMesh=strcmp(Parameters(iD1).Problem,'Mesh');
  isStructure=strcmp(Parameters(iD1).Problem,'Structural');
  isFluidDM=strcmp(Parameters(iD2).Formulation,'WeaklyCompressibleFlowDM_HDG');
  isFluidVP=strcmp(Parameters(iD2).Formulation,'WeaklyCompressibleFlowVP_HDG');
  isFluidFCFV=strcmp(Parameters(iD2).Formulation,'IncompressibleFlow_FCFV');
end
iElem1=Faces(iD1,iD2).Interface(iFaceInterface,1);
iFace1=Faces(iD1,iD2).Interface(iFaceInterface,2);
iFace2=Faces(iD1,iD2).Interface(iFaceInterface,4);
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes1=Sizes(iD1).NumElementNodes;
NumElementNodes2=Sizes(iD2).NumElementNodes;
NumElementFaces1=Sizes(iD1).NumElementFaces;
NumElementFaces2=Sizes(iD2).NumElementFaces;
NumFaceNodes2=Sizes(iD2).NumFaceNodes;
C1e=Mesh(iD1).Elements(:,iElem1)';
X1e=Mesh(iD1).Nodes(:,C1e)';
X1em=sum(X1e(1:NumElementFaces1,:),1)/NumElementFaces1;
gamma=Parameters(iD1).NitschePenalty;
if isFSI && isMesh && matchField(Parameters(iD1),'RelaxationParameter')
  omega=Parameters(iD1).RelaxationParameter;
else
  omega=1;
end

% Get solution
u1e=reshape(Block(iD1,iD1).SolutionGlobal(C1e,:),[],1);
if isFSI && isStructure && isFluidDM
  U2e=Block(iD2,iD2).SolutionGlobal;
end

% Initialize lhs
if isStructural
  Ku1U2=zeros(nsd*NumElementNodes1,nsd*NumElementFaces2*NumFaceNodes2);
elseif isFSI && isMesh
  Ku1u2=zeros(nsd*NumElementNodes1,nsd*NumElementNodes2);
elseif isFSI && isStructure && (isFluidVP || isFluidFCFV)
  Ku1V2=zeros(nsd*NumElementNodes1,nsd*NumElementFaces2*NumFaceNodes2);
elseif isFSI && isStructure && isFluidDM
  Ku1R2=zeros(nsd*NumElementNodes1,NumElementFaces2*NumFaceNodes2);
  Ku1W2=zeros(nsd*NumElementNodes1,nsd*NumElementFaces2*NumFaceNodes2);
end

% Compute weights at Gauss points
[~,N1ex,N1ey,N1ez,~,~,pinvN1e]=mapShapeFunctions(1,RefElement(iD1,iD1),RefElement(iD1,iD1),X1e,nsd);

% Indices
n1e1=1:NumElementNodes1;
n1e2=n1e1+NumElementNodes1;
n1e3=n1e2+NumElementNodes1;

% Compute variables at nodes
u1xe=u1e(n1e1);
u1ye=u1e(n1e2);
if nsd==3
  u1ze=u1e(n1e3);
end
if nsd==2
  F1xxe=pinvN1e*(N1ex*u1xe+1); F1xye=pinvN1e*(N1ey*u1xe);
  F1yxe=pinvN1e*(N1ex*u1ye);   F1yye=pinvN1e*(N1ey*u1ye+1);
elseif nsd==3
  F1xxe=pinvN1e*(N1ex*u1xe+1); F1xye=pinvN1e*(N1ey*u1xe);   F1xze=pinvN1e*(N1ez*u1xe);
  F1yxe=pinvN1e*(N1ex*u1ye);   F1yye=pinvN1e*(N1ey*u1ye+1); F1yze=pinvN1e*(N1ez*u1ye);
  F1zxe=pinvN1e*(N1ex*u1ze);   F1zye=pinvN1e*(N1ey*u1ze);   F1zze=pinvN1e*(N1ez*u1ze+1);
end

% Compute weights at Gauss points
FaceNodes1=RefElement(iD1,iD1).FaceNodesElem;
FaceNodes2=RefElement(iD2,iD2).FaceNodesElem;
X1f=X1e(FaceNodes1(iFace1,:),:);
[~,~,~,~,w1fg]=mapShapeFunctions(0,RefElement(iD1,iD1),RefElement(iD1,iD1),X1f,nsd);
[N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions(0,RefElement(iD1,iD2),RefElement(iD1,iD2),X1f,nsd);
N21f=RefElement(iD2,iD1).ShapeFunctionsFace;

% Compute characteristic element size
h=sum(w1fg);

% Indices
n1f1=FaceNodes1(iFace1,:);
n1f2=n1f1+NumElementNodes1;
n1f3=n1f2+NumElementNodes1;
if isStructural
  n2ef1=(iFace2-1)*nsd*NumFaceNodes2+(1:NumFaceNodes2);
elseif isFSI && isMesh
  n2f1=FaceNodes2(iFace2,:);
  n2f2=n2f1+NumElementNodes2;
  n2f3=n2f2+NumElementNodes2;
elseif isFSI && isStructure
  n2ef1=(iFace2-1)*(nsd+1)*NumFaceNodes2+(1:NumFaceNodes2);
end
if isStructural || (isFSI && isStructure)
  n2ef2=n2ef1+NumFaceNodes2;
  n2ef3=n2ef2+NumFaceNodes2;
  n2ef4=n2ef3+NumFaceNodes2;
end
if isFSI && isStructure && (isFluidVP || isFluidFCFV)
  n2efV1=(iFace2-1)*nsd*NumFaceNodes2+(1:NumFaceNodes2);
  n2efV2=n2efV1+NumFaceNodes2;
  n2efV3=n2efV2+NumFaceNodes2;
elseif isFSI && isStructure && isFluidDM
  n2efR1=(iFace2-1)*NumFaceNodes2+(1:NumFaceNodes2);
  n2efW1=(iFace2-1)*nsd*NumFaceNodes2+(1:NumFaceNodes2);
  n2efW2=n2efW1+NumFaceNodes2;
  n2efW3=n2efW2+NumFaceNodes2;
end

% Flip face
Node2Match1stNode1=Faces(iD1,iD2).Interface(iFaceInterface,5);
order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
if isStructural || (isFSI && isStructure)
  n2ef1=n2ef1(order);
  n2ef2=n2ef2(order);
  n2ef3=n2ef3(order);
  n2ef4=n2ef4(order);
elseif isFSI && isMesh
  n2f1=n2f1(order);
  n2f2=n2f2(order);
  n2f3=n2f3(order);
end
if isFSI && isStructure && (isFluidVP || isFluidFCFV)
  n2efV1=n2efV1(order);
  n2efV2=n2efV2(order);
  n2efV3=n2efV3(order);
elseif isFSI && isStructure && isFluidDM
  n2efR1=n2efR1(order);
  n2efW1=n2efW1(order);
  n2efW2=n2efW2(order);
  n2efW3=n2efW3(order);
end

% Compute derivatives of shape functions
N1xe=pinvN1e*N1ex;
N1ye=pinvN1e*N1ey;
if nsd==3
  N1ze=pinvN1e*N1ez;
end
N12xf=N12f*N1xe(n1f1,:);
N12yf=N12f*N1ye(n1f1,:);
if nsd==3
  N12zf=N12f*N1ze(n1f1,:);
end

% Compute variables at nodes
if nsd==2
  F1xxf=F1xxe(n1f1); F1xyf=F1xye(n1f1);
  F1yxf=F1yxe(n1f1); F1yyf=F1yye(n1f1);
elseif nsd==3
  F1xxf=F1xxe(n1f1); F1xyf=F1xye(n1f1); F1xzf=F1xze(n1f1);
  F1yxf=F1yxe(n1f1); F1yyf=F1yye(n1f1); F1yzf=F1yze(n1f1);
  F1zxf=F1zxe(n1f1); F1zyf=F1zye(n1f1); F1zzf=F1zze(n1f1);
end
if isFSI && isStructure && isFluidDM
  R2f=U2e(n2ef1);
  W2xf=U2e(n2ef2);
  W2yf=U2e(n2ef3);
  if nsd==3
    W2zf=U2e(n2ef4);
  end
end

% Compute variables at Gauss points
if nsd==2
  F1xxfg=N12f*F1xxf; F1xyfg=N12f*F1xyf;
  F1yxfg=N12f*F1yxf; F1yyfg=N12f*F1yyf;
elseif nsd==3
  F1xxfg=N12f*F1xxf; F1xyfg=N12f*F1xyf; F1xzfg=N12f*F1xzf;
  F1yxfg=N12f*F1yxf; F1yyfg=N12f*F1yyf; F1yzfg=N12f*F1yzf;
  F1zxfg=N12f*F1zxf; F1zyfg=N12f*F1zyf; F1zzfg=N12f*F1zzf;
end
if isFSI && isStructure && isFluidDM
  R2fg=N21f*R2f;
  W2xfg=N21f*W2xf;
  W2yfg=N21f*W2yf;
  if nsd==3
    W2zfg=N21f*W2zf;
  end
end

% Compute basic matrices
Nw12fT=(w12fg.*N12f)';
M12f=Nw12fT*N21f;
Nw12xfT=(w12fg.*N12xf)';
Nw12yfT=(w12fg.*N12yf)';
if nsd==3
  Nw12zfT=(w12fg.*N12zf)';
end

% Compute stress linearization
if nsd==2
  [~,ds1dF1fg]=...
    computeStress('no','yes','no',Parameters(iD1),X1em,...
    [F1xxfg,F1xyfg,F1yxfg,F1yyfg],Sizes(iD1));
  ds1xxdF1xxfg=ds1dF1fg(:,1);  ds1xxdF1xyfg=ds1dF1fg(:,2);  ds1xxdF1yxfg=ds1dF1fg(:,3);  ds1xxdF1yyfg=ds1dF1fg(:,4);
                               ds1xydF1xyfg=ds1dF1fg(:,6);  ds1xydF1yxfg=ds1dF1fg(:,7);  ds1xydF1yyfg=ds1dF1fg(:,8);
                                                            ds1yxdF1yxfg=ds1dF1fg(:,11); ds1yxdF1yyfg=ds1dF1fg(:,12);
                                                                                         ds1yydF1yyfg=ds1dF1fg(:,16);
elseif nsd==3
  [~,ds1dF1fg]=...
    computeStress('no','yes','no',Parameters(iD1),X1em,...
    [F1xxfg,F1xyfg,F1xzfg,F1yxfg,F1yyfg,F1yzfg,F1zxfg,F1zyfg,F1zzfg],Sizes(iD1));
  ds1xxdF1xxfg=ds1dF1fg(:,1);  ds1xxdF1xyfg=ds1dF1fg(:,2);  ds1xxdF1xzfg=ds1dF1fg(:,3);  ds1xxdF1yxfg=ds1dF1fg(:,4);  ds1xxdF1yyfg=ds1dF1fg(:,5);  ds1xxdF1yzfg=ds1dF1fg(:,6);  ds1xxdF1zxfg=ds1dF1fg(:,7);  ds1xxdF1zyfg=ds1dF1fg(:,8);  ds1xxdF1zzfg=ds1dF1fg(:,9);
                               ds1xydF1xyfg=ds1dF1fg(:,11); ds1xydF1xzfg=ds1dF1fg(:,12); ds1xydF1yxfg=ds1dF1fg(:,13); ds1xydF1yyfg=ds1dF1fg(:,14); ds1xydF1yzfg=ds1dF1fg(:,15); ds1xydF1zxfg=ds1dF1fg(:,16); ds1xydF1zyfg=ds1dF1fg(:,17); ds1xydF1zzfg=ds1dF1fg(:,18);
                                                            ds1xzdF1xzfg=ds1dF1fg(:,21); ds1xzdF1yxfg=ds1dF1fg(:,22); ds1xzdF1yyfg=ds1dF1fg(:,23); ds1xzdF1yzfg=ds1dF1fg(:,24); ds1xzdF1zxfg=ds1dF1fg(:,25); ds1xzdF1zyfg=ds1dF1fg(:,26); ds1xzdF1zzfg=ds1dF1fg(:,27);
                                                                                         ds1yxdF1yxfg=ds1dF1fg(:,31); ds1yxdF1yyfg=ds1dF1fg(:,32); ds1yxdF1yzfg=ds1dF1fg(:,33); ds1yxdF1zxfg=ds1dF1fg(:,34); ds1yxdF1zyfg=ds1dF1fg(:,35); ds1yxdF1zzfg=ds1dF1fg(:,36);
                                                                                                                      ds1yydF1yyfg=ds1dF1fg(:,41); ds1yydF1yzfg=ds1dF1fg(:,42); ds1yydF1zxfg=ds1dF1fg(:,43); ds1yydF1zyfg=ds1dF1fg(:,44); ds1yydF1zzfg=ds1dF1fg(:,45);
                                                                                                                                                   ds1yzdF1yzfg=ds1dF1fg(:,51); ds1yzdF1zxfg=ds1dF1fg(:,52); ds1yzdF1zyfg=ds1dF1fg(:,53); ds1yzdF1zzfg=ds1dF1fg(:,54);
                                                                                                                                                                                ds1zxdF1zxfg=ds1dF1fg(:,61); ds1zxdF1zyfg=ds1dF1fg(:,62); ds1zxdF1zzfg=ds1dF1fg(:,63);
                                                                                                                                                                                                             ds1zydF1zyfg=ds1dF1fg(:,71); ds1zydF1zzfg=ds1dF1fg(:,72);
                                                                                                                                                                                                                                          ds1zzdF1zzfg=ds1dF1fg(:,81);
end

% Compute lhs
if isStructural
  Ku1U2(n1f1,n2ef1)=Ku1U2(n1f1,n2ef1)-gamma/h*M12f;
  Ku1U2(n1f2,n2ef2)=Ku1U2(n1f2,n2ef2)-gamma/h*M12f;
  if nsd==3
    Ku1U2(n1f3,n2ef3)=Ku1U2(n1f3,n2ef3)-gamma/h*M12f;
  end
  
  if nsd==2
    Ku1U2(n1e1,n2ef1)=Ku1U2(n1e1,n2ef1)...
                     +Nw12xfT*((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y).*N21f);
    Ku1U2(n1e1,n2ef2)=Ku1U2(n1e1,n2ef2)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y).*N21f);
    Ku1U2(n1e2,n2ef1)=Ku1U2(n1e2,n2ef1)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y).*N21f);
    Ku1U2(n1e2,n2ef2)=Ku1U2(n1e2,n2ef2)...
                     +Nw12xfT*((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y).*N21f);
  elseif nsd==3
    Ku1U2(n1e1,n2ef1)=Ku1U2(n1e1,n2ef1)...
                     +Nw12xfT*((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y+ds1xxdF1xzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y+ds1xydF1xzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xxdF1xzfg.*n12x+ds1xydF1xzfg.*n12y+ds1xzdF1xzfg.*n12z).*N21f);
    Ku1U2(n1e1,n2ef2)=Ku1U2(n1e1,n2ef2)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y+ds1xxdF1yzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y+ds1xydF1yzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xzdF1yxfg.*n12x+ds1xzdF1yyfg.*n12y+ds1xzdF1yzfg.*n12z).*N21f);
    Ku1U2(n1e1,n2ef3)=Ku1U2(n1e1,n2ef3)...
                     +Nw12xfT*((+ds1xxdF1zxfg.*n12x+ds1xxdF1zyfg.*n12y+ds1xxdF1zzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xydF1zxfg.*n12x+ds1xydF1zyfg.*n12y+ds1xydF1zzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xzdF1zxfg.*n12x+ds1xzdF1zyfg.*n12y+ds1xzdF1zzfg.*n12z).*N21f);
    Ku1U2(n1e2,n2ef1)=Ku1U2(n1e2,n2ef1)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y+ds1xzdF1yxfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y+ds1xzdF1yyfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xxdF1yzfg.*n12x+ds1xydF1yzfg.*n12y+ds1xzdF1yzfg.*n12z).*N21f);
    Ku1U2(n1e2,n2ef2)=Ku1U2(n1e2,n2ef2)...
                     +Nw12xfT*((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y+ds1yxdF1yzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y+ds1yydF1yzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1yxdF1yzfg.*n12x+ds1yydF1yzfg.*n12y+ds1yzdF1yzfg.*n12z).*N21f);
    Ku1U2(n1e2,n2ef3)=Ku1U2(n1e2,n2ef3)...
                     +Nw12xfT*((+ds1yxdF1zxfg.*n12x+ds1yxdF1zyfg.*n12y+ds1yxdF1zzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1yydF1zxfg.*n12x+ds1yydF1zyfg.*n12y+ds1yydF1zzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1yzdF1zxfg.*n12x+ds1yzdF1zyfg.*n12y+ds1yzdF1zzfg.*n12z).*N21f);
    Ku1U2(n1e3,n2ef1)=Ku1U2(n1e3,n2ef1)...
                     +Nw12xfT*((+ds1xxdF1zxfg.*n12x+ds1xydF1zxfg.*n12y+ds1xzdF1zxfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xxdF1zyfg.*n12x+ds1xydF1zyfg.*n12y+ds1xzdF1zyfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xxdF1zzfg.*n12x+ds1xydF1zzfg.*n12y+ds1xzdF1zzfg.*n12z).*N21f);
    Ku1U2(n1e3,n2ef2)=Ku1U2(n1e3,n2ef2)...
                     +Nw12xfT*((+ds1yxdF1zxfg.*n12x+ds1yydF1zxfg.*n12y+ds1yzdF1zxfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1yxdF1zyfg.*n12x+ds1yydF1zyfg.*n12y+ds1yzdF1zyfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1yxdF1zzfg.*n12x+ds1yydF1zzfg.*n12y+ds1yzdF1zzfg.*n12z).*N21f);
    Ku1U2(n1e3,n2ef3)=Ku1U2(n1e3,n2ef3)...
                     +Nw12xfT*((+ds1zxdF1zxfg.*n12x+ds1zxdF1zyfg.*n12y+ds1zxdF1zzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1zxdF1zyfg.*n12x+ds1zydF1zyfg.*n12y+ds1zydF1zzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1zxdF1zzfg.*n12x+ds1zydF1zzfg.*n12y+ds1zzdF1zzfg.*n12z).*N21f);
  end
end

if isFSI && isMesh
  Ku1u2(n1f1,n2f1)=Ku1u2(n1f1,n2f1)-gamma/h*omega*M12f;
  Ku1u2(n1f2,n2f2)=Ku1u2(n1f2,n2f2)-gamma/h*omega*M12f;
  if nsd==3
    Ku1u2(n1f3,n2f3)=Ku1u2(n1f3,n2f3)-gamma/h*omega*M12f;
  end
  
  if nsd==2
    Ku1u2(n1e1,n2f1)=Ku1u2(n1e1,n2f1)...
               +omega*Nw12xfT*((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y).*N21f)...
               +omega*Nw12yfT*((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y).*N21f);
    Ku1u2(n1e1,n2f2)=Ku1u2(n1e1,n2f2)...
               +omega*Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y).*N21f)...
               +omega*Nw12yfT*((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y).*N21f);
    Ku1u2(n1e2,n2f1)=Ku1u2(n1e2,n2f1)...
               +omega*Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y).*N21f)...
               +omega*Nw12yfT*((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y).*N21f);
    Ku1u2(n1e2,n2f2)=Ku1u2(n1e2,n2f2)...
               +omega*Nw12xfT*((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y).*N21f)...
               +omega*Nw12yfT*((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y).*N21f);
  elseif nsd==3
    Ku1u2(n1e1,n2f1)=Ku1u2(n1e1,n2f1)...
               +omega*Nw12xfT*((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y+ds1xxdF1xzfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y+ds1xydF1xzfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1xxdF1xzfg.*n12x+ds1xydF1xzfg.*n12y+ds1xzdF1xzfg.*n12z).*N21f);
    Ku1u2(n1e1,n2f2)=Ku1u2(n1e1,n2f2)...
               +omega*Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y+ds1xxdF1yzfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y+ds1xydF1yzfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1xzdF1yxfg.*n12x+ds1xzdF1yyfg.*n12y+ds1xzdF1yzfg.*n12z).*N21f);
    Ku1u2(n1e1,n2f3)=Ku1u2(n1e1,n2f3)...
               +omega*Nw12xfT*((+ds1xxdF1zxfg.*n12x+ds1xxdF1zyfg.*n12y+ds1xxdF1zzfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1xydF1zxfg.*n12x+ds1xydF1zyfg.*n12y+ds1xydF1zzfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1xzdF1zxfg.*n12x+ds1xzdF1zyfg.*n12y+ds1xzdF1zzfg.*n12z).*N21f);
    Ku1u2(n1e2,n2f1)=Ku1u2(n1e2,n2f1)...
               +omega*Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y+ds1xzdF1yxfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y+ds1xzdF1yyfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1xxdF1yzfg.*n12x+ds1xydF1yzfg.*n12y+ds1xzdF1yzfg.*n12z).*N21f);
    Ku1u2(n1e2,n2f2)=Ku1u2(n1e2,n2f2)...
               +omega*Nw12xfT*((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y+ds1yxdF1yzfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y+ds1yydF1yzfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1yxdF1yzfg.*n12x+ds1yydF1yzfg.*n12y+ds1yzdF1yzfg.*n12z).*N21f);
    Ku1u2(n1e2,n2f3)=Ku1u2(n1e2,n2f3)...
               +omega*Nw12xfT*((+ds1yxdF1zxfg.*n12x+ds1yxdF1zyfg.*n12y+ds1yxdF1zzfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1yydF1zxfg.*n12x+ds1yydF1zyfg.*n12y+ds1yydF1zzfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1yzdF1zxfg.*n12x+ds1yzdF1zyfg.*n12y+ds1yzdF1zzfg.*n12z).*N21f);
    Ku1u2(n1e3,n2f1)=Ku1u2(n1e3,n2f1)...
               +omega*Nw12xfT*((+ds1xxdF1zxfg.*n12x+ds1xydF1zxfg.*n12y+ds1xzdF1zxfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1xxdF1zyfg.*n12x+ds1xydF1zyfg.*n12y+ds1xzdF1zyfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1xxdF1zzfg.*n12x+ds1xydF1zzfg.*n12y+ds1xzdF1zzfg.*n12z).*N21f);
    Ku1u2(n1e3,n2f2)=Ku1u2(n1e3,n2f2)...
               +omega*Nw12xfT*((+ds1yxdF1zxfg.*n12x+ds1yydF1zxfg.*n12y+ds1yzdF1zxfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1yxdF1zyfg.*n12x+ds1yydF1zyfg.*n12y+ds1yzdF1zyfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1yxdF1zzfg.*n12x+ds1yydF1zzfg.*n12y+ds1yzdF1zzfg.*n12z).*N21f);
    Ku1u2(n1e3,n2f3)=Ku1u2(n1e3,n2f3)...
               +omega*Nw12xfT*((+ds1zxdF1zxfg.*n12x+ds1zxdF1zyfg.*n12y+ds1zxdF1zzfg.*n12z).*N21f)...
               +omega*Nw12yfT*((+ds1zxdF1zyfg.*n12x+ds1zydF1zyfg.*n12y+ds1zydF1zzfg.*n12z).*N21f)...
               +omega*Nw12zfT*((+ds1zxdF1zzfg.*n12x+ds1zydF1zzfg.*n12y+ds1zzdF1zzfg.*n12z).*N21f);
  end
end

if isFSI && isStructure && (isFluidVP || isFluidFCFV)
  Ku1V2(n1f1,n2efV1)=Ku1V2(n1f1,n2efV1)-gamma/h*M12f;
  Ku1V2(n1f2,n2efV2)=Ku1V2(n1f2,n2efV2)-gamma/h*M12f;
  if nsd==3
    Ku1V2(n1f3,n2efV3)=Ku1V2(n1f3,n2efV3)-gamma/h*M12f;
  end
  
  if nsd==2
    Ku1V2(n1e1,n2efV1)=Ku1V2(n1e1,n2efV1)...
                     +Nw12xfT*((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y).*N21f);
    Ku1V2(n1e1,n2efV2)=Ku1V2(n1e1,n2efV2)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y).*N21f);
    Ku1V2(n1e2,n2efV1)=Ku1V2(n1e2,n2efV1)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y).*N21f);
    Ku1V2(n1e2,n2efV2)=Ku1V2(n1e2,n2efV2)...
                     +Nw12xfT*((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y).*N21f)...
                     +Nw12yfT*((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y).*N21f);
  elseif nsd==3
    Ku1V2(n1e1,n2efV1)=Ku1V2(n1e1,n2efV1)...
                     +Nw12xfT*((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y+ds1xxdF1xzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y+ds1xydF1xzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xxdF1xzfg.*n12x+ds1xydF1xzfg.*n12y+ds1xzdF1xzfg.*n12z).*N21f);
    Ku1V2(n1e1,n2efV2)=Ku1V2(n1e1,n2efV2)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y+ds1xxdF1yzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y+ds1xydF1yzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xzdF1yxfg.*n12x+ds1xzdF1yyfg.*n12y+ds1xzdF1yzfg.*n12z).*N21f);
    Ku1V2(n1e1,n2efV3)=Ku1V2(n1e1,n2efV3)...
                     +Nw12xfT*((+ds1xxdF1zxfg.*n12x+ds1xxdF1zyfg.*n12y+ds1xxdF1zzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xydF1zxfg.*n12x+ds1xydF1zyfg.*n12y+ds1xydF1zzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xzdF1zxfg.*n12x+ds1xzdF1zyfg.*n12y+ds1xzdF1zzfg.*n12z).*N21f);
    Ku1V2(n1e2,n2efV1)=Ku1V2(n1e2,n2efV1)...
                     +Nw12xfT*((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y+ds1xzdF1yxfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y+ds1xzdF1yyfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xxdF1yzfg.*n12x+ds1xydF1yzfg.*n12y+ds1xzdF1yzfg.*n12z).*N21f);
    Ku1V2(n1e2,n2efV2)=Ku1V2(n1e2,n2efV2)...
                     +Nw12xfT*((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y+ds1yxdF1yzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y+ds1yydF1yzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1yxdF1yzfg.*n12x+ds1yydF1yzfg.*n12y+ds1yzdF1yzfg.*n12z).*N21f);
    Ku1V2(n1e2,n2efV3)=Ku1V2(n1e2,n2efV3)...
                     +Nw12xfT*((+ds1yxdF1zxfg.*n12x+ds1yxdF1zyfg.*n12y+ds1yxdF1zzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1yydF1zxfg.*n12x+ds1yydF1zyfg.*n12y+ds1yydF1zzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1yzdF1zxfg.*n12x+ds1yzdF1zyfg.*n12y+ds1yzdF1zzfg.*n12z).*N21f);
    Ku1V2(n1e3,n2efV1)=Ku1V2(n1e3,n2efV1)...
                     +Nw12xfT*((+ds1xxdF1zxfg.*n12x+ds1xydF1zxfg.*n12y+ds1xzdF1zxfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1xxdF1zyfg.*n12x+ds1xydF1zyfg.*n12y+ds1xzdF1zyfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1xxdF1zzfg.*n12x+ds1xydF1zzfg.*n12y+ds1xzdF1zzfg.*n12z).*N21f);
    Ku1V2(n1e3,n2efV2)=Ku1V2(n1e3,n2efV2)...
                     +Nw12xfT*((+ds1yxdF1zxfg.*n12x+ds1yydF1zxfg.*n12y+ds1yzdF1zxfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1yxdF1zyfg.*n12x+ds1yydF1zyfg.*n12y+ds1yzdF1zyfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1yxdF1zzfg.*n12x+ds1yydF1zzfg.*n12y+ds1yzdF1zzfg.*n12z).*N21f);
    Ku1V2(n1e3,n2efV3)=Ku1V2(n1e3,n2efV3)...
                     +Nw12xfT*((+ds1zxdF1zxfg.*n12x+ds1zxdF1zyfg.*n12y+ds1zxdF1zzfg.*n12z).*N21f)...
                     +Nw12yfT*((+ds1zxdF1zyfg.*n12x+ds1zydF1zyfg.*n12y+ds1zydF1zzfg.*n12z).*N21f)...
                     +Nw12zfT*((+ds1zxdF1zzfg.*n12x+ds1zydF1zzfg.*n12y+ds1zzdF1zzfg.*n12z).*N21f);
  end
end

if isFSI && isStructure && isFluidDM
  Ku1R2(n1f1,n2efR1)=Ku1R2(n1f1,n2efR1)+gamma/h*Nw12fT*((W2xfg./R2fg.^2).*N21f);
  Ku1R2(n1f2,n2efR1)=Ku1R2(n1f2,n2efR1)+gamma/h*Nw12fT*((W2yfg./R2fg.^2).*N21f);
  if nsd==3
    Ku1R2(n1f3,n2efR1)=Ku1R2(n1f3,n2efR1)+gamma/h*Nw12fT*((W2zfg./R2fg.^2).*N21f);
  end
  
  Ku1W2(n1f1,n2efW1)=Ku1W2(n1f1,n2efW1)-gamma/h*Nw12fT*((1./R2fg).*N21f);
  Ku1W2(n1f2,n2efW2)=Ku1W2(n1f2,n2efW2)-gamma/h*Nw12fT*((1./R2fg).*N21f);
  if nsd==3
    Ku1W2(n1f3,n2efW3)=Ku1W2(n1f3,n2efW3)-gamma/h*Nw12fT*((1./R2fg).*N21f);
  end
  
  if nsd==2
    Ku1R2(n1e1,n2efR1)=Ku1R2(n1e1,n2efR1)...
-Nw12xfT*(((+(+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y).*W2xfg...
            +(+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y).*W2yfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y).*W2xfg...
            +(+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y).*W2yfg)./R2fg.^2).*N21f);
    Ku1R2(n1e2,n2efR1)=Ku1R2(n1e2,n2efR1)...
-Nw12xfT*(((+(+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y).*W2xfg...
            +(+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y).*W2yfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y).*W2xfg...
            +(+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y).*W2yfg)./R2fg.^2).*N21f);
  elseif nsd==3
    Ku1R2(n1e1,n2efR1)=Ku1R2(n1e1,n2efR1)...
-Nw12xfT*(((+(+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y+ds1xxdF1xzfg.*n12z).*W2xfg...
            +(+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y+ds1xxdF1yzfg.*n12z).*W2yfg...
            +(+ds1xxdF1zxfg.*n12x+ds1xxdF1zyfg.*n12y+ds1xxdF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y+ds1xydF1xzfg.*n12z).*W2xfg...
            +(+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y+ds1xydF1yzfg.*n12z).*W2yfg...
            +(+ds1xydF1zxfg.*n12x+ds1xydF1zyfg.*n12y+ds1xydF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f)...
-Nw12zfT*(((+(+ds1xxdF1xzfg.*n12x+ds1xydF1xzfg.*n12y+ds1xzdF1xzfg.*n12z).*W2xfg...
            +(+ds1xzdF1yxfg.*n12x+ds1xzdF1yyfg.*n12y+ds1xzdF1yzfg.*n12z).*W2yfg...
            +(+ds1xzdF1zxfg.*n12x+ds1xzdF1zyfg.*n12y+ds1xzdF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f);
    Ku1R2(n1e2,n2efR1)=Ku1R2(n1e2,n2efR1)...
-Nw12xfT*(((+(+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y+ds1xzdF1yxfg.*n12z).*W2xfg...
            +(+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y+ds1yxdF1yzfg.*n12z).*W2yfg...
            +(+ds1yxdF1zxfg.*n12x+ds1yxdF1zyfg.*n12y+ds1yxdF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y+ds1xzdF1yyfg.*n12z).*W2xfg...
            +(+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y+ds1yydF1yzfg.*n12z).*W2yfg...
            +(+ds1yydF1zxfg.*n12x+ds1yydF1zyfg.*n12y+ds1yydF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1yzfg.*n12x+ds1xydF1yzfg.*n12y+ds1xzdF1yzfg.*n12z).*W2xfg...
            +(+ds1yxdF1yzfg.*n12x+ds1yydF1yzfg.*n12y+ds1yzdF1yzfg.*n12z).*W2yfg...
            +(+ds1yzdF1zxfg.*n12x+ds1yzdF1zyfg.*n12y+ds1yzdF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f);
    Ku1R2(n1e3,n2efR1)=Ku1R2(n1e3,n2efR1)...
-Nw12xfT*(((+(+ds1xxdF1zxfg.*n12x+ds1xydF1zxfg.*n12y+ds1xzdF1zxfg.*n12z).*W2xfg...
            +(+ds1yxdF1zxfg.*n12x+ds1yydF1zxfg.*n12y+ds1yzdF1zxfg.*n12z).*W2yfg...
            +(+ds1zxdF1zxfg.*n12x+ds1zxdF1zyfg.*n12y+ds1zxdF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1zyfg.*n12x+ds1xydF1zyfg.*n12y+ds1xzdF1zyfg.*n12z).*W2xfg...
            +(+ds1yxdF1zyfg.*n12x+ds1yydF1zyfg.*n12y+ds1yzdF1zyfg.*n12z).*W2yfg...
            +(+ds1zxdF1zyfg.*n12x+ds1zydF1zyfg.*n12y+ds1zydF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f)...
-Nw12yfT*(((+(+ds1xxdF1zzfg.*n12x+ds1xydF1zzfg.*n12y+ds1xzdF1zzfg.*n12z).*W2xfg...
            +(+ds1yxdF1zzfg.*n12x+ds1yydF1zzfg.*n12y+ds1yzdF1zzfg.*n12z).*W2yfg...
            +(+ds1zxdF1zzfg.*n12x+ds1zydF1zzfg.*n12y+ds1zzdF1zzfg.*n12z).*W2zfg)./R2fg.^2).*N21f);
  end
  
  if nsd==2
    Ku1W2(n1e1,n2efW1)=Ku1W2(n1e1,n2efW1)...
             +Nw12xfT*(((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y)./R2fg).*N21f);
    Ku1W2(n1e1,n2efW2)=Ku1W2(n1e1,n2efW2)...
             +Nw12xfT*(((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y)./R2fg).*N21f);
    Ku1W2(n1e2,n2efW1)=Ku1W2(n1e2,n2efW1)...
             +Nw12xfT*(((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y)./R2fg).*N21f);
    Ku1W2(n1e2,n2efW2)=Ku1W2(n1e2,n2efW2)...
             +Nw12xfT*(((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y)./R2fg).*N21f);
  elseif nsd==3
    Ku1W2(n1e1,n2efW1)=Ku1W2(n1e1,n2efW1)...
             +Nw12xfT*(((+ds1xxdF1xxfg.*n12x+ds1xxdF1xyfg.*n12y+ds1xxdF1xzfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xxdF1xyfg.*n12x+ds1xydF1xyfg.*n12y+ds1xydF1xzfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1xxdF1xzfg.*n12x+ds1xydF1xzfg.*n12y+ds1xzdF1xzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e1,n2efW2)=Ku1W2(n1e1,n2efW2)...
             +Nw12xfT*(((+ds1xxdF1yxfg.*n12x+ds1xxdF1yyfg.*n12y+ds1xxdF1yzfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xydF1yxfg.*n12x+ds1xydF1yyfg.*n12y+ds1xydF1yzfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1xzdF1yxfg.*n12x+ds1xzdF1yyfg.*n12y+ds1xzdF1yzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e1,n2efW3)=Ku1W2(n1e1,n2efW3)...
             +Nw12xfT*(((+ds1xxdF1zxfg.*n12x+ds1xxdF1zyfg.*n12y+ds1xxdF1zzfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xydF1zxfg.*n12x+ds1xydF1zyfg.*n12y+ds1xydF1zzfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1xzdF1zxfg.*n12x+ds1xzdF1zyfg.*n12y+ds1xzdF1zzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e2,n2efW1)=Ku1W2(n1e2,n2efW1)...
             +Nw12xfT*(((+ds1xxdF1yxfg.*n12x+ds1xydF1yxfg.*n12y+ds1xzdF1yxfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xxdF1yyfg.*n12x+ds1xydF1yyfg.*n12y+ds1xzdF1yyfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1xxdF1yzfg.*n12x+ds1xydF1yzfg.*n12y+ds1xzdF1yzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e2,n2efW2)=Ku1W2(n1e2,n2efW2)...
             +Nw12xfT*(((+ds1yxdF1yxfg.*n12x+ds1yxdF1yyfg.*n12y+ds1yxdF1yzfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1yxdF1yyfg.*n12x+ds1yydF1yyfg.*n12y+ds1yydF1yzfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1yxdF1yzfg.*n12x+ds1yydF1yzfg.*n12y+ds1yzdF1yzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e2,n2efW3)=Ku1W2(n1e2,n2efW3)...
             +Nw12xfT*(((+ds1yxdF1zxfg.*n12x+ds1yxdF1zyfg.*n12y+ds1yxdF1zzfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1yydF1zxfg.*n12x+ds1yydF1zyfg.*n12y+ds1yydF1zzfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1yzdF1zxfg.*n12x+ds1yzdF1zyfg.*n12y+ds1yzdF1zzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e3,n2efW1)=Ku1W2(n1e3,n2efW1)...
             +Nw12xfT*(((+ds1xxdF1zxfg.*n12x+ds1xydF1zxfg.*n12y+ds1xzdF1zxfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1xxdF1zyfg.*n12x+ds1xydF1zyfg.*n12y+ds1xzdF1zyfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1xxdF1zzfg.*n12x+ds1xydF1zzfg.*n12y+ds1xzdF1zzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e3,n2efW2)=Ku1W2(n1e3,n2efW2)...
             +Nw12xfT*(((+ds1yxdF1zxfg.*n12x+ds1yydF1zxfg.*n12y+ds1yzdF1zxfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1yxdF1zyfg.*n12x+ds1yydF1zyfg.*n12y+ds1yzdF1zyfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1yxdF1zzfg.*n12x+ds1yydF1zzfg.*n12y+ds1yzdF1zzfg.*n12z)./R2fg).*N21f);
    Ku1W2(n1e3,n2efW3)=Ku1W2(n1e3,n2efW3)...
             +Nw12xfT*(((+ds1zxdF1zxfg.*n12x+ds1zxdF1zyfg.*n12y+ds1zxdF1zzfg.*n12z)./R2fg).*N21f)...
             +Nw12yfT*(((+ds1zxdF1zyfg.*n12x+ds1zydF1zyfg.*n12y+ds1zydF1zzfg.*n12z)./R2fg).*N21f)...
             +Nw12zfT*(((+ds1zxdF1zzfg.*n12x+ds1zydF1zzfg.*n12y+ds1zzdF1zzfg.*n12z)./R2fg).*N21f);
  end
end

% Indices
if isFSI && isStructure && (isFluidVP || isFluidFCFV)
  iu1=1:nsd*NumElementNodes1;
  iV2=reshape((0:NumElementFaces2-1)*(nsd+1)*NumFaceNodes2+repmat((1:nsd*NumFaceNodes2)',...
    1,NumElementFaces2),1,[]);
elseif isFSI && isStructure && isFluidDM
  iu1=1:nsd*NumElementNodes1;
  iR2=reshape((0:NumElementFaces2-1)*(1+nsd)*NumFaceNodes2+repmat((1:NumFaceNodes2)',...
    1,NumElementFaces2),1,[]);
  iW2=reshape((0:NumElementFaces2-1)*(1+nsd)*NumFaceNodes2+repmat((1:nsd*NumFaceNodes2)',...
    1,NumElementFaces2),1,[])+NumFaceNodes2;
end

% Initialization of lhs and rhs
if isFSI && isStructure
  LhsCoup=zeros(nsd*NumElementNodes1,(nsd+1)*NumElementFaces2*NumFaceNodes2);
end

% Compute elemental contributions to lhs
if isStructural
  LhsCoup=Ku1U2;
elseif isFSI && isMesh
  LhsCoup=Ku1u2;
elseif isFSI && isStructure && (isFluidVP || isFluidFCFV)
  LhsCoup(iu1,iV2)=Ku1V2;
elseif isFSI && isStructure && isFluidDM
  LhsCoup(iu1,iR2)=Ku1R2;
  LhsCoup(iu1,iW2)=Ku1W2;
end

end

%% Evaluate stress element
function [Stress]=evaluateStressElement(...
  iD,Nodes,Results,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes(iD).NumSpaceDim;
NumElementNodes=Sizes(iD).NumElementNodes;
NumElementFaces=Sizes(iD).NumElementFaces;
Xe=Nodes';
Xem=sum(Xe(1:NumElementFaces,:),1)/NumElementFaces;

% Get solution
ue=reshape(Results,[],1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,~,~,pinvNe]=mapShapeFunctions(1,RefElement(iD,iD),RefElement(iD,iD),Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
uxe=ue(ne1);
uye=ue(ne2);
if nsd==3
  uze=ue(ne3);
end
if nsd==2
  Fxxe=pinvNe*(Nex*uxe+1); Fxye=pinvNe*(Ney*uxe);
  Fyxe=pinvNe*(Nex*uye);   Fyye=pinvNe*(Ney*uye+1);
elseif nsd==3
  Fxxe=pinvNe*(Nex*uxe+1); Fxye=pinvNe*(Ney*uxe);   Fxze=pinvNe*(Nez*uxe);
  Fyxe=pinvNe*(Nex*uye);   Fyye=pinvNe*(Ney*uye+1); Fyze=pinvNe*(Nez*uye);
  Fzxe=pinvNe*(Nex*uze);   Fzye=pinvNe*(Ney*uze);   Fzze=pinvNe*(Nez*uze+1);  
end

% Compute variables at Gauss points
if nsd==2
  Fxxeg=Ne*Fxxe; Fxyeg=Ne*Fxye;
  Fyxeg=Ne*Fyxe; Fyyeg=Ne*Fyye;
elseif nsd==3
  Fxxeg=Ne*Fxxe; Fxyeg=Ne*Fxye; Fxzeg=Ne*Fxze;
  Fyxeg=Ne*Fyxe; Fyyeg=Ne*Fyye; Fyzeg=Ne*Fyze;
  Fzxeg=Ne*Fzxe; Fzyeg=Ne*Fzye; Fzzeg=Ne*Fzze;
end

% Compute stress
if nsd==2
  [seg]=...
    computeStress('yes','no','no',Parameters(iD),Xem,...
    [Fxxeg,Fxyeg,Fyxeg,Fyyeg],Sizes(iD));
  sxxeg=seg(:,1); sxyeg=seg(:,2);
  syxeg=seg(:,3); syyeg=seg(:,4);
elseif nsd==3
  [seg]=...
    computeStress('yes','no','no',Parameters(iD),Xem,...
    [Fxxeg,Fxyeg,Fxzeg,Fyxeg,Fyyeg,Fyzeg,Fzxeg,Fzyeg,Fzzeg],Sizes(iD));
  sxxeg=seg(:,1); sxyeg=seg(:,2); sxzeg=seg(:,3);
  syxeg=seg(:,4); syyeg=seg(:,5); syzeg=seg(:,6);
  szxeg=seg(:,7); szyeg=seg(:,8); szzeg=seg(:,9);
end

% Map stress to element nodes
if nsd==2
  sxxe=pinvNe*(sxxeg); sxye=pinvNe*(sxyeg);
  syxe=pinvNe*(syxeg); syye=pinvNe*(syyeg);
elseif nsd==3
  sxxe=pinvNe*(sxxeg); sxye=pinvNe*(sxyeg);  sxze=pinvNe*(sxzeg);
  syxe=pinvNe*(syxeg); syye=pinvNe*(syyeg);  syze=pinvNe*(syzeg);
  szxe=pinvNe*(szxeg); szye=pinvNe*(szyeg);  szze=pinvNe*(szzeg);
end

% Store stress
if nsd==2
  Stress=[sxxe,sxye,syxe,syye];
elseif nsd==3
  Stress=[sxxe,sxye,sxze,syxe,syye,syze,szxe,szye,szze];
end

end