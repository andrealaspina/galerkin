classdef WeaklyCompressibleFlowDM_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1+NumSpaceDim;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2+1+NumSpaceDim;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) NumSpaceDim;

    % Discretization type
    DiscretizationType='HDG';

    % Time derivative order
    TimeDerOrder=1;
    
    % Time/frequency domain
    Domain='Time';
    
  end
  
  methods
    
    %% Initialize unknowns
    function [Block]=initializeUnknowns(~,iD,Block,Parameters,Time,Sizes)
      Block(iD,iD).SolutionGlobal=zeros(Sizes(iD).NumGlobalNodes*Sizes(iD).NumGlobalComp,1);
      Block(iD,iD).SolutionLocal=zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumLocalComp,1);
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionOld=zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumLocalComp,...
          Parameters(iD).TimeDerOrder);
      end
      
      % Consider reference density instead of zero
      for iFace=1:Sizes(iD).NumFaces
        Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp.*Sizes(iD).NumFaceNodes+...
          (1:Sizes(iD).NumFaceNodes))=Parameters(iD).ReferenceDensity;
      end
      Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+1)=Parameters(iD).ReferenceDensity;
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
          Parameters(iD).Density(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).Momentum(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
              Parameters(iD).Density(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).Momentum(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,Simulation,Parameters,~,Faces,Time,...
        RefElement,Sizes)
      NodesElem=Elements(iD1).Nodes;
      NodesInitialElem=Elements(iD1).NodesInitial;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionLocalElem=Elements(iD1).SolutionLocal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      MatLocal=cell(Sizes(iD1).NumElements,1);
      VecLocal=cell(Sizes(iD1).NumElements,1);
      
      % Extract coupling data (coinciding elements)
      SolutionGlobalSameElemCoupled=cell(Sizes(iD1).NumElements,1);
      SolutionOldSameElemCoupled=cell(Sizes(iD1).NumElements,1);
      if strcmp(Simulation.Problem,'FluidALE') || ...
         strcmp(Simulation.Problem,'FluidStructureInteraction')
        iD2=find(contains({Parameters.Problem},'Mesh'));
        SolutionGlobalSameElemCoupled=Elements(iD2).SolutionGlobal;
        SolutionOldSameElemCoupled=Elements(iD2).SolutionOld;
      end
      
      % Extract coupling data
      NodesElemCoupled=double.empty(Sizes.NumElements,0);
      SolutionGlobalElemCoupled=double.empty(Sizes.NumElements,0);
      SolutionOldElemCoupled=double.empty(Sizes.NumElements,0);
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        if strcmp(Simulation.Problem,'FluidStructureInteraction')
          iD2=find(contains({Parameters.Problem},'Structural'));
        end
        NodesElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        SolutionGlobalElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        SolutionOldElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        iEF=sub2ind([Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces],...
          Faces(iD1,iD2).Interface(:,1),Faces(iD1,iD2).Interface(:,2));
        NodesElemCoupled(iEF)=Elements(iD2).Nodes(Faces(iD1,iD2).Interface(:,3));
        SolutionGlobalElemCoupled(iEF)=Elements(iD2).SolutionGlobal(Faces(iD1,iD2).Interface(:,3));
        SolutionOldElemCoupled(iEF)=Elements(iD2).SolutionOld(Faces(iD1,iD2).Interface(:,3));
      end
      
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem,MatLocalElem,VecLocalElem]=...
          buildBlockElement(iD1,NodesElem{iElem},NodesInitialElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          SolutionGlobalSameElemCoupled{iElem},...
          SolutionOldSameElemCoupled{iElem},...
          NodesElemCoupled(iElem,:),SolutionGlobalElemCoupled(iElem,:),...
          SolutionOldElemCoupled(iElem,:),Simulation,Parameters,Time,RefElement.Value,Sizes); %#ok
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
        MatLocal{iElem}=MatLocalElem;
        VecLocal{iElem}=VecLocalElem;
      end
      Block(iD1,iD1).LhsGlobal=fsparse(Block(iD1,iD1).LhsRowIndices,...
                                       Block(iD1,iD1).LhsColIndices,LhsCoef(:));
      Block(iD1,iD1).RhsGlobal=fsparse(Block(iD1,iD1).RhsRowIndices,1,RhsCoef(:));
      Elements(iD1).MatLocal=MatLocal;
      Elements(iD1).VecLocal=VecLocal;
    end
    
    %% Do post-process
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,Time,RefElement,Sizes)
      NodesElem=Elements.Nodes;
      FacesElem=Elements.Faces;
      SolutionGlobalElem=Elements.SolutionGlobal;
      SolutionLocalElem=Elements.SolutionLocal;
      LhsPost=cell(Sizes.NumElements,1);
      RhsPost=cell(Sizes.NumElements,1);
      parfor iElem=1:Sizes.NumElements
        [LhsPostElem,RhsPostElem]=...
          doPostProcessElement(NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},...
          Parameters,Time,RefElement,Sizes);
        LhsPost{iElem}=LhsPostElem;
        RhsPost{iElem}=RhsPostElem;
      end
      Elements.LhsPost=LhsPost;
      Elements.RhsPost=RhsPost;
    end
    
    %% Do coupling
    function [Block]=doCoupling(~,iD1,iD2,Block,~,Simulation,Parameters,Mesh,Faces,Time,...
        RefElement,Sizes)
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface)) && ...
         (strcmp(Simulation.Problem,'FluidStructureInteraction') && ...
          strcmp(Parameters(iD2).Problem,'Structural'))
        LhsCoupCoef=zeros(Sizes(iD1).NumElementLhsCoupCoef(iD2),Sizes(iD1).NumFacesInterface(iD2));
        for iFaceInterface=1:Sizes(iD1).NumFacesInterface(iD2)
          [LhsCoupElem]=...
            doCouplingElement(iFaceInterface,iD1,iD2,Block,...
            Parameters,Mesh,Faces,Time,RefElement.Value,Sizes);
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
    function [Results]=storeResults(~,iD,iST,Results,Block,Simulation,Parameters,Mesh,Time,...
        RefElement,Sizes)
      if iST==1
        Results(iD).Time=[];
        Results(iD).ScaledStrainRate=[];
        Results(iD).Density=[];
        Results(iD).Momentum=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledStrainRate(:,:,iST)=Block(iD,iD).SolutionLocal(:,1:Sizes(iD).NumVoigtComp);
      Results(iD).Density(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+1);
      Results(iD).Momentum(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+1+...
        (1:Sizes(iD).NumSpaceDim));
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).VelocityPost=Block(iD,iD).SolutionPost;
      end
      if strcmp(Parameters(iD).ArbitraryLagrangianEulerian,'yes')
        for iElem=1:Sizes(iD).NumElements
          if strcmp(Simulation.Problem,'Fluid')
            Xe=Mesh(iD).Nodes(:,Mesh(iD).Elements(:,iElem)')';
            ue=Parameters(iD).Displacement(Xe(:,1),Xe(:,2),Xe(:,3),Time.Time);
          elseif strcmp(Simulation.Problem,'FluidALE') || ...
                 strcmp(Simulation.Problem,'FluidStructureInteraction')
            iD2=find(contains({Parameters.Problem},'Mesh'));
            u2e=Block(iD2,iD2).SolutionGlobal(Mesh(iD2).Elements(:,iElem)',:);
            ue=RefElement(iD,iD2).ShapeFunctionsElem\(RefElement(iD2,iD).ShapeFunctionsElem*u2e);
          end
          Results(iD).Displacement(...
            (iElem-1)*Sizes(iD).NumElementNodes+(1:Sizes(iD).NumElementNodes),:,iST)=ue;
        end
      end
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  iD1,Nodes,NodesInitial,Faces,SolutionGlobal,SolutionLocal,SolutionOld,...
  SolutionGlobalSameCoupled,SolutionOldSameCoupled,...
  NodesCoupled,SolutionGlobalCoupled,SolutionOldCoupled,...
  Simulation,Parameters,Time,RefElement,Sizes)

% Get general parameters
isFluid=strcmp(Simulation.Problem,'Fluid');
isFluidALE=strcmp(Simulation.Problem,'FluidALE');
isFSI=strcmp(Simulation.Problem,'FluidStructureInteraction');
isConvectiveFlow=strcmp(Parameters(iD1).ConvectiveFlow,'yes');
isArbitraryLagrangianEulerian=strcmp(Parameters(iD1).ArbitraryLagrangianEulerian,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes(iD1).NumSpaceDim;
msd=nsd*(nsd+1)/2;
NumElementNodes=Sizes(iD1).NumElementNodes;
NumElementFaces=Sizes(iD1).NumElementFaces;
NumFaceNodes=Sizes(iD1).NumFaceNodes;
mu=Parameters(iD1).DynamicViscosity;
lambda=-2/3*Parameters(iD1).DynamicViscosity;
p=Parameters(iD1).Pressure;
dpdr=Parameters(iD1).DPressureDDensity;
rD=Parameters(iD1).Density;
wD=Parameters(iD1).Momentum;
tN=Parameters(iD1).Traction;
Rc=Parameters(iD1).ResidualContinuity;
f=Parameters(iD1).Force;
Xe=Nodes';
X0e=NodesInitial';
tauR=Parameters(iD1).StabDensity;
tauW=Parameters(iD1).StabMomentum;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
re=reshape(SolutionLocal(:,msd+1),[],1);
we=reshape(SolutionLocal(:,msd+1+(1:nsd)),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  rolde=reshape(SolutionOld(:,msd+1,:),[],BDFo);
  wolde=reshape(SolutionOld(:,msd+1+(1:nsd),:),[],BDFo);
end
if isArbitraryLagrangianEulerian
  if isFluid
    ue=reshape(Parameters(iD1).Displacement(X0e(:,1),X0e(:,2),X0e(:,3),t),[],1);
    if isTimeDependent
      uolde=zeros(nsd*NumElementNodes,BDFo);
      for iBDF=1:BDFo
        uolde(:,iBDF)=reshape(Parameters(iD1).Displacement(X0e(:,1),X0e(:,2),X0e(:,3),t-iBDF*dt),...
          [],1);
      end
    end
  elseif isFluidALE || isFSI
    iD2=find(contains({Parameters.Problem},'Mesh'));
    ue=reshape(SolutionGlobalSameCoupled,[],1);
    if isTimeDependent
      uolde=reshape(SolutionOldSameCoupled(:,:,1:BDFo),[],BDFo);
    end
  end
end

% Initialize lhs
KLL=zeros(msd*NumElementNodes,msd*NumElementNodes);
KLr=zeros(msd*NumElementNodes,NumElementNodes);
KLw=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KLR=zeros(msd*NumElementNodes,NumElementFaces*NumFaceNodes);
KLW=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Krr=zeros(NumElementNodes,NumElementNodes);
Krw=zeros(NumElementNodes,nsd*NumElementNodes);
KrR=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KrW=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KwL=zeros(nsd*NumElementNodes,msd*NumElementNodes);
Kwr=zeros(nsd*NumElementNodes,NumElementNodes);
Kww=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KwR=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
KwW=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KRr=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KRR=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KWL=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KWw=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KWR=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KWW=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);

% Initialize rhs
fL=zeros(msd*NumElementNodes,1);
fr=zeros(NumElementNodes,1);
fw=zeros(nsd*NumElementNodes,1);
fR=zeros(NumElementFaces*NumFaceNodes,1);
fW=zeros(nsd*NumElementFaces*NumFaceNodes,1);

% Update nodes coordinates
if isArbitraryLagrangianEulerian
  if isFluid
    Xe(:,1:nsd)=X0e(:,1:nsd)+reshape(ue,[],nsd);
  elseif isFluidALE || isFSI
    N210e=RefElement(iD2,iD1).ShapeFunctionsElem;
    pinvN102e=RefElement(iD1,iD2).PseudoinverseShapeFunctionsElem;
    Xe(:,1:nsd)=X0e(:,1:nsd)+pinvN102e*(N210e*reshape(ue,[],nsd));
  end
end

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions('Element',RefElement(iD1,iD1),RefElement(iD1,iD1),Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;
ne4=ne3+NumElementNodes;
ne5=ne4+NumElementNodes;
ne6=ne5+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Lxxe=Le(ne1);
  Lyye=Le(ne2);
  Lxye=Le(ne3);
elseif nsd==3
  Lxxe=Le(ne1);
  Lyye=Le(ne2);
  Lzze=Le(ne3);
  Lxye=Le(ne4);
  Lxze=Le(ne5);
  Lyze=Le(ne6);
end
wxe=we(ne1);
wye=we(ne2);
if nsd==3
  wze=we(ne3);
end
if isTimeDependent
  woldxe=wolde(ne1,:);
  woldye=wolde(ne2,:);
  if nsd==3
    woldze=wolde(ne3,:);
  end
end
pe=p(re);
if isArbitraryLagrangianEulerian
  if isFluid
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
  elseif isFluidALE || isFSI
    n2e1=1:Sizes(iD2).NumElementNodes;
    n2e2=n2e1+Sizes(iD2).NumElementNodes;
    n2e3=n2e2+Sizes(iD2).NumElementNodes;
    uxe=ue(n2e1);
    uye=ue(n2e2);
    if nsd==3
      uze=ue(n2e3);
    end
    if isTimeDependent
      uoldxe=uolde(n2e1,:);
      uoldye=uolde(n2e2,:);
      if nsd==3
        uoldze=uolde(n2e3,:);
      end
    end
  end
end
if isArbitraryLagrangianEulerian && isTimeDependent
  if isFluid
    axe=1/dt*uxe*alpha(1)+1/dt*uoldxe*alpha(2:BDFo+1,1);
    aye=1/dt*uye*alpha(1)+1/dt*uoldye*alpha(2:BDFo+1,1);
    if nsd==3
      aze=1/dt*uze*alpha(1)+1/dt*uoldze*alpha(2:BDFo+1,1);
    end
  elseif isFluidALE || isFSI
    axe=pinvN102e*(N210e*(1/dt*uxe*alpha(1)+1/dt*uoldxe*alpha(2:BDFo+1,1)));
    aye=pinvN102e*(N210e*(1/dt*uye*alpha(1)+1/dt*uoldye*alpha(2:BDFo+1,1)));
    if nsd==3
      aze=pinvN102e*(N210e*(1/dt*uze*alpha(1)+1/dt*uoldze*alpha(2:BDFo+1,1)));
    end
  end
else
  axe=wxe*0;
  aye=wye*0;
  if nsd==3
    aze=wze*0;
  end
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
if nsd==2
  Lxxeg=Ne*Lxxe;
  Lyyeg=Ne*Lyye;
  Lxyeg=Ne*Lxye;
elseif nsd==3
  Lxxeg=Ne*Lxxe;
  Lyyeg=Ne*Lyye;
  Lzzeg=Ne*Lzze;
  Lxyeg=Ne*Lxye;
  Lxzeg=Ne*Lxze;
  Lyzeg=Ne*Lyze;
end
reg=Ne*re;
wxeg=Ne*wxe;
wyeg=Ne*wye;
if nsd==3
  wzeg=Ne*wze;
end
if isTimeDependent
  roldeg=Ne*rolde;
  woldxeg=Ne*woldxe;
  woldyeg=Ne*woldye;
  if nsd==3
    woldzeg=Ne*woldze;
  end
end
Rceg=Rc(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
feg=f(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
fxeg=feg(:,1);
fyeg=feg(:,2);
if nsd==3
  fzeg=feg(:,3);
end
dpdreg=dpdr(reg);
axeg=Ne*axe;
ayeg=Ne*aye;
if nsd==3
  azeg=Ne*aze;
end
if isArbitraryLagrangianEulerian
  if nsd==2
    divaeg=Nex*axe+Ney*aye;
  else
    divaeg=Nex*axe+Ney*aye+Nez*aze;
  end
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
Cxe=NwexT*Ne;
Cye=NweyT*Ne;
if nsd==3
  Cze=NwezT*Ne;
end
CxeT=Cxe';
CyeT=Cye';
if nsd==3
  CzeT=Cze';
end

% Compute linearization of viscous stress
if nsd==2
  Voigt1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  Voigt2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  Voigt3=sqrt(mu);
elseif nsd==3
  Voigt1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  Voigt2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  Voigt3=sqrt(mu);
end

% Compute lhs
KLL(ne1,ne1)=-Me;
KLL(ne2,ne2)=-Me;
KLL(ne3,ne3)=-Me;
if nsd==3
  KLL(ne4,ne4)=-Me;
  KLL(ne5,ne5)=-Me;
  KLL(ne6,ne6)=-Me;
end

if nsd==2
  KLr(ne1,ne1)=-Voigt1*NwexT*((wxeg./reg.^2).*Ne)...
               -Voigt2*NweyT*((wyeg./reg.^2).*Ne);
  KLr(ne2,ne1)=-Voigt2*NwexT*((wxeg./reg.^2).*Ne)...
               -Voigt1*NweyT*((wyeg./reg.^2).*Ne);
  KLr(ne3,ne1)=-Voigt3*NweyT*((wxeg./reg.^2).*Ne)...
               -Voigt3*NwexT*((wyeg./reg.^2).*Ne);
elseif nsd==3
  KLr(ne1,ne1)=-Voigt1*NwexT*((wxeg./reg.^2).*Ne)...
               -Voigt2*NweyT*((wyeg./reg.^2).*Ne)...
               -Voigt2*NwezT*((wzeg./reg.^2).*Ne);
  KLr(ne2,ne1)=-Voigt2*NwexT*((wxeg./reg.^2).*Ne)...
               -Voigt1*NweyT*((wyeg./reg.^2).*Ne)...
               -Voigt2*NwezT*((wzeg./reg.^2).*Ne);
  KLr(ne3,ne1)=-Voigt2*NwexT*((wxeg./reg.^2).*Ne)...
               -Voigt2*NweyT*((wyeg./reg.^2).*Ne)...
               -Voigt1*NwezT*((wzeg./reg.^2).*Ne);
  KLr(ne4,ne1)=-Voigt3*NweyT*((wxeg./reg.^2).*Ne)...
               -Voigt3*NwexT*((wyeg./reg.^2).*Ne);
  KLr(ne5,ne1)=-Voigt3*NwezT*((wxeg./reg.^2).*Ne)...
               -Voigt3*NwexT*((wzeg./reg.^2).*Ne);
  KLr(ne6,ne1)=-Voigt3*NwezT*((wyeg./reg.^2).*Ne)...
               -Voigt3*NweyT*((wzeg./reg.^2).*Ne);
end

if nsd==2
  KLw(ne1,ne1)=Voigt1*NwexT*((1./reg).*Ne);
  KLw(ne2,ne1)=Voigt2*NwexT*((1./reg).*Ne);
  KLw(ne3,ne1)=Voigt3*NweyT*((1./reg).*Ne);
  KLw(ne1,ne2)=Voigt2*NweyT*((1./reg).*Ne);
  KLw(ne2,ne2)=Voigt1*NweyT*((1./reg).*Ne);
  KLw(ne3,ne2)=Voigt3*NwexT*((1./reg).*Ne);
elseif nsd==3
  KLw(ne1,ne1)=Voigt1*NwexT*((1./reg).*Ne);
  KLw(ne2,ne1)=Voigt2*NwexT*((1./reg).*Ne);
  KLw(ne3,ne1)=Voigt2*NwexT*((1./reg).*Ne);
  KLw(ne4,ne1)=Voigt3*NweyT*((1./reg).*Ne);
  KLw(ne5,ne1)=Voigt3*NwezT*((1./reg).*Ne);
  KLw(ne1,ne2)=Voigt2*NweyT*((1./reg).*Ne);
  KLw(ne2,ne2)=Voigt1*NweyT*((1./reg).*Ne);
  KLw(ne3,ne2)=Voigt2*NweyT*((1./reg).*Ne);
  KLw(ne4,ne2)=Voigt3*NwexT*((1./reg).*Ne);
  KLw(ne6,ne2)=Voigt3*NwezT*((1./reg).*Ne);
  KLw(ne1,ne3)=Voigt2*NwezT*((1./reg).*Ne);
  KLw(ne2,ne3)=Voigt2*NwezT*((1./reg).*Ne);
  KLw(ne3,ne3)=Voigt1*NwezT*((1./reg).*Ne);
  KLw(ne5,ne3)=Voigt3*NwexT*((1./reg).*Ne);
  KLw(ne6,ne3)=Voigt3*NweyT*((1./reg).*Ne);
end

if isTimeDependent
  Krr(ne1,ne1)=alpha(1)/dt*Me;
end

if isArbitraryLagrangianEulerian
  Krr(ne1,ne1)=Krr(ne1,ne1)+NweT*((divaeg).*Ne);
end

Krr(ne1,ne1)=Krr(ne1,ne1)+NwexT*((axeg).*Ne)...
                         +NweyT*((ayeg).*Ne);
if nsd==3
  Krr(ne1,ne1)=Krr(ne1,ne1)+NwezT*((azeg).*Ne);
end

Krw(ne1,ne1)=-Cxe;
Krw(ne1,ne2)=-Cye;
if nsd==3
  Krw(ne1,ne3)=-Cze;
end

if isTimeDependent
  Kww(ne1,ne1)=alpha(1)/dt*Me;
  Kww(ne2,ne2)=alpha(1)/dt*Me;
  if nsd==3
    Kww(ne3,ne3)=alpha(1)/dt*Me;
  end
end

if isArbitraryLagrangianEulerian
  Kww(ne1,ne1)=Kww(ne1,ne1)+NweT*((divaeg).*Ne);
  Kww(ne2,ne2)=Kww(ne2,ne2)+NweT*((divaeg).*Ne);
  if nsd==3
    Kww(ne3,ne3)=Kww(ne3,ne3)+NweT*((divaeg).*Ne);
  end
end

if isConvectiveFlow
  if nsd==2
    Kwr(ne1,ne1)=+NwexT*((wxeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wxeg.*wyeg./reg.^2).*Ne);
    Kwr(ne2,ne1)=+NwexT*((wyeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wyeg.*wyeg./reg.^2).*Ne);
  elseif nsd==3
    Kwr(ne1,ne1)=+NwexT*((wxeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wxeg.*wyeg./reg.^2).*Ne)...
                 +NwezT*((wxeg.*wzeg./reg.^2).*Ne);
    Kwr(ne2,ne1)=+NwexT*((wyeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wyeg.*wyeg./reg.^2).*Ne)...
                 +NwezT*((wyeg.*wzeg./reg.^2).*Ne);
    Kwr(ne3,ne1)=+NwexT*((wzeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wzeg.*wyeg./reg.^2).*Ne)...
                 +NwezT*((wzeg.*wzeg./reg.^2).*Ne);
  end
  
  if nsd==2
    Kww(ne1,ne1)=Kww(ne1,ne1)-NwexT*((2*wxeg./reg-axeg).*Ne)...
                             -NweyT*((  wyeg./reg-ayeg).*Ne);
    Kww(ne2,ne1)=Kww(ne2,ne1)-NwexT*((wyeg./reg).*Ne);
    Kww(ne1,ne2)=Kww(ne1,ne2)-NweyT*((wxeg./reg).*Ne);
    Kww(ne2,ne2)=Kww(ne2,ne2)-NwexT*((  wxeg./reg-axeg).*Ne)...
                             -NweyT*((2*wyeg./reg-ayeg).*Ne);
  elseif nsd==3
    Kww(ne1,ne1)=Kww(ne1,ne1)-NwexT*((2*wxeg./reg-axeg).*Ne)...
                             -NweyT*((  wyeg./reg-ayeg).*Ne)...
                             -NwezT*((  wzeg./reg-azeg).*Ne);
    Kww(ne2,ne1)=Kww(ne2,ne1)-NwexT*((wyeg./reg).*Ne);
    Kww(ne3,ne1)=Kww(ne3,ne1)-NwezT*((wzeg./reg).*Ne);
    Kww(ne1,ne2)=Kww(ne1,ne2)-NweyT*((wxeg./reg).*Ne);
    Kww(ne2,ne2)=Kww(ne2,ne2)-NwexT*((  wxeg./reg-axeg).*Ne)...
                             -NweyT*((2*wyeg./reg-ayeg).*Ne)...
                             -NwezT*((  wzeg./reg-azeg).*Ne);
    Kww(ne3,ne2)=Kww(ne3,ne2)-NwezT*((wzeg./reg).*Ne);
    Kww(ne1,ne3)=Kww(ne1,ne3)-NwezT*((wxeg./reg).*Ne);
    Kww(ne2,ne3)=Kww(ne2,ne3)-NwezT*((wyeg./reg).*Ne);
    Kww(ne3,ne3)=Kww(ne3,ne3)-NwexT*((  wxeg./reg-axeg).*Ne)...
                             -NweyT*((  wyeg./reg-ayeg).*Ne)...
                             -NwezT*((2*wzeg./reg-azeg).*Ne);
  end
end

if nsd==2
  KwL(ne1,ne1)=Voigt1*CxeT;
  KwL(ne2,ne1)=Voigt2*CyeT;
  KwL(ne1,ne2)=Voigt2*CxeT;
  KwL(ne2,ne2)=Voigt1*CyeT;
  KwL(ne1,ne3)=Voigt3*CyeT;
  KwL(ne2,ne3)=Voigt3*CxeT;
elseif nsd==3
  KwL(ne1,ne1)=Voigt1*CxeT;
  KwL(ne2,ne1)=Voigt2*CyeT;
  KwL(ne3,ne1)=Voigt2*CzeT;
  KwL(ne1,ne2)=Voigt2*CxeT;
  KwL(ne2,ne2)=Voigt1*CyeT;
  KwL(ne3,ne2)=Voigt2*CzeT;
  KwL(ne1,ne3)=Voigt2*CxeT;
  KwL(ne2,ne3)=Voigt2*CyeT;
  KwL(ne3,ne3)=Voigt1*CzeT;
  KwL(ne1,ne4)=Voigt3*CyeT;
  KwL(ne2,ne4)=Voigt3*CxeT;
  KwL(ne1,ne5)=Voigt3*CzeT;
  KwL(ne3,ne5)=Voigt3*CxeT;
  KwL(ne2,ne6)=Voigt3*CzeT;
  KwL(ne3,ne6)=Voigt3*CyeT;
end

Kwr(ne1,ne1)=Kwr(ne1,ne1)+NweT*((dpdreg).*Nex);
Kwr(ne2,ne1)=Kwr(ne2,ne1)+NweT*((dpdreg).*Ney);
if nsd==3
  Kwr(ne3,ne1)=Kwr(ne3,ne1)+NweT*((dpdreg).*Nez);
end

% Compute rhs
if nsd==2
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(Voigt1*wxeg./reg)...
            -NweyT*(Voigt2*wyeg./reg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(Voigt2*wxeg./reg)...
            -NweyT*(Voigt1*wyeg./reg);
  fL(ne3,1)=+NweT*(Lxyeg)...
            -NweyT*(Voigt3*wxeg./reg)...
            -NwexT*(Voigt3*wyeg./reg);
elseif nsd==3
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(Voigt1*wxeg./reg)...
            -NweyT*(Voigt2*wyeg./reg)...
            -NwezT*(Voigt2*wzeg./reg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(Voigt2*wxeg./reg)...
            -NweyT*(Voigt1*wyeg./reg)...
            -NwezT*(Voigt2*wzeg./reg);
  fL(ne3,1)=+NweT*(Lzzeg)...
            -NwexT*(Voigt2*wxeg./reg)...
            -NweyT*(Voigt2*wyeg./reg)...
            -NwezT*(Voigt1*wzeg./reg);
  fL(ne4,1)=+NweT*(Lxyeg)...
            -NwexT*(Voigt3*wyeg./reg)...
            -NweyT*(Voigt3*wxeg./reg);
  fL(ne5,1)=+NweT*(Lxzeg)...
            -NwexT*(Voigt3*wzeg./reg)...
            -NwezT*(Voigt3*wxeg./reg);
  fL(ne6,1)=+NweT*(Lyzeg)...
            -NweyT*(Voigt3*wzeg./reg)...
            -NwezT*(Voigt3*wyeg./reg);
end

if isTimeDependent
  fr(ne1,1)=-NweT*(1/dt*reg*alpha(1)...
                  +1/dt*roldeg*alpha(2:BDFo+1,1));
end

if isArbitraryLagrangianEulerian
  fr(ne1,1)=fr(ne1,1)-NweT*(reg.*divaeg);
end

fr(ne1,1)=fr(ne1,1)+NwexT*(wxeg-reg.*axeg)...
                   +NweyT*(wyeg-reg.*ayeg);
if nsd==3
  fr(ne1,1)=fr(ne1,1)+NwezT*(wzeg-reg.*azeg);
end

fr(ne1,1)=fr(ne1,1)+NweT*(Rceg);

if isTimeDependent
  fw(ne1,1)=-NweT*(1/dt*wxeg*alpha(1)...
                  +1/dt*woldxeg*alpha(2:BDFo+1,1));
  fw(ne2,1)=-NweT*(1/dt*wyeg*alpha(1)...
                  +1/dt*woldyeg*alpha(2:BDFo+1,1));
  if nsd==3
    fw(ne3,1)=-NweT*(1/dt*wzeg*alpha(1)...
                    +1/dt*woldzeg*alpha(2:BDFo+1,1));
  end
end

if isArbitraryLagrangianEulerian
  fw(ne1,1)=fw(ne1,1)-NweT*(wxeg.*divaeg);
  fw(ne2,1)=fw(ne2,1)-NweT*(wyeg.*divaeg);
  if nsd==3
    fw(ne3,1)=fw(ne3,1)-NweT*(wzeg.*divaeg);
  end
end

if isConvectiveFlow
  if nsd==2
    fw(ne1,1)=fw(ne1,1)+NwexT*(wxeg.*(wxeg./reg-axeg))...
                       +NweyT*(wxeg.*(wyeg./reg-ayeg));
    fw(ne2,1)=fw(ne2,1)+NwexT*(wyeg.*(wxeg./reg-axeg))...
                       +NweyT*(wyeg.*(wyeg./reg-ayeg));
  elseif nsd==3
    fw(ne1,1)=fw(ne1,1)+NwexT*(wxeg.*(wxeg./reg-axeg))...
                       +NweyT*(wxeg.*(wyeg./reg-ayeg))...
                       +NwezT*(wxeg.*(wzeg./reg-azeg));
    fw(ne2,1)=fw(ne2,1)+NwexT*(wyeg.*(wxeg./reg-axeg))...
                       +NweyT*(wyeg.*(wyeg./reg-ayeg))...
                       +NwezT*(wyeg.*(wzeg./reg-azeg));
    fw(ne3,1)=fw(ne3,1)+NwexT*(wzeg.*(wxeg./reg-axeg))...
                       +NweyT*(wzeg.*(wyeg./reg-ayeg))...
                       +NwezT*(wzeg.*(wzeg./reg-azeg));
  end
end

if nsd==2
  fw(ne1,1)=fw(ne1,1)-NweT*(+Voigt1*(Nex*Lxxe)+Voigt2*(Nex*Lyye)...
                            +Voigt3*(Ney*Lxye)...
                            +(Nex*pe));
  fw(ne2,1)=fw(ne2,1)-NweT*(+Voigt2*(Ney*Lxxe)+Voigt1*(Ney*Lyye)...
                            +Voigt3*(Nex*Lxye)...
                            +(Ney*pe));
elseif nsd==3
  fw(ne1,1)=fw(ne1,1)-NweT*(+Voigt1*(Nex*Lxxe)+Voigt2*(Nex*Lyye)+Voigt2*(Nex*Lzze)...
                            +Voigt3*(Ney*Lxye)+Voigt3*(Nez*Lxze)...
                            +(Nex*pe));
  fw(ne2,1)=fw(ne2,1)-NweT*(+Voigt2*(Ney*Lxxe)+Voigt1*(Ney*Lyye)+Voigt2*(Ney*Lzze)...
                            +Voigt3*(Nex*Lxye)+Voigt3*(Nez*Lyze)...
                            +(Ney*pe));
  fw(ne3,1)=fw(ne3,1)-NweT*(+Voigt2*(Nez*Lxxe)+Voigt2*(Nez*Lyye)+Voigt1*(Nez*Lzze)...
                            +Voigt3*(Nex*Lxze)+Voigt3*(Ney*Lyze)...
                            +(Nez*pe));
end

fw(ne1,1)=fw(ne1,1)+NweT*(fxeg);
fw(ne2,1)=fw(ne2,1)+NweT*(fyeg);
if nsd==3
  fw(ne3,1)=fw(ne3,1)+NweT*(fzeg);
end

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=true;
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement(iD1,iD1).FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions('Face',RefElement(iD1,iD1),RefElement(iD1,iD1),Xf,nsd);
    
    % Check boundary
    isExterior=Faces.Exterior(iFace);
    isDirichlet_r=Faces.Dirichlet_r(iFace);
    isDirichlet_w_x=Faces.Dirichlet_w_x(iFace);
    isDirichlet_w_y=Faces.Dirichlet_w_y(iFace);
    if nsd==3; isDirichlet_w_z=Faces.Dirichlet_w_z(iFace); end
    isNeumann_t_x=Faces.Neumann_t_x(iFace);
    isNeumann_t_y=Faces.Neumann_t_y(iFace);
    if nsd==3; isNeumann_t_z=Faces.Neumann_t_z(iFace); end
    isDirichlet_w=isDirichlet_w_x || isDirichlet_w_y || (nsd==3 && isDirichlet_w_z);
    isNeumann_t=isNeumann_t_x || isNeumann_t_y || (nsd==3 && isNeumann_t_z);
    if isFSI
      isInterface=Faces.Interface(1,iFace);
    end
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nf4=nf3+NumElementNodes;
    nf5=nf4+NumElementNodes;
    nf6=nf5+NumElementNodes;
    nefU1=(iFace-1)*(1+nsd)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
    nefU4=nefU3+NumFaceNodes;
    nefR1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefW1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefW2=nefW1+NumFaceNodes;
    nefW3=nefW2+NumFaceNodes;
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters(iD1).Degree,Node2Match1stNode1);
      nefU1=nefU1(order);
      nefU2=nefU2(order);
      nefU3=nefU3(order);
      nefU4=nefU4(order);
      nefR1=nefR1(order);
      nefW1=nefW1(order);
      nefW2=nefW2(order);
      nefW3=nefW3(order);
    end
    
    % Compute variables at nodes
    if nsd==2
      Lxxf=Lxxe(nf1);
      Lyyf=Lyye(nf1);
      Lxyf=Lxye(nf1);
    elseif nsd==3
      Lxxf=Lxxe(nf1);
      Lyyf=Lyye(nf1);
      Lzzf=Lzze(nf1);
      Lxyf=Lxye(nf1);
      Lxzf=Lxze(nf1);
      Lyzf=Lyze(nf1);
    end
    rf=re(nf1);
    wxf=wxe(nf1);
    wyf=wye(nf1);
    if nsd==3
      wzf=wze(nf1);
    end
    axf=axe(nf1);
    ayf=aye(nf1);
    if nsd==3
      azf=aze(nf1);
    end
    Rf=Ue(nefU1);
    Wxf=Ue(nefU2);
    Wyf=Ue(nefU3);
    if nsd==3
      Wzf=Ue(nefU4);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    if nsd==2
      Lxxfg=Nf*Lxxf;
      Lyyfg=Nf*Lyyf;
      Lxyfg=Nf*Lxyf;
    elseif nsd==3
      Lxxfg=Nf*Lxxf;
      Lyyfg=Nf*Lyyf;
      Lzzfg=Nf*Lzzf;
      Lxyfg=Nf*Lxyf;
      Lxzfg=Nf*Lxzf;
      Lyzfg=Nf*Lyzf;
    end
    rfg=Nf*rf;
    wxfg=Nf*wxf;
    wyfg=Nf*wyf;
    if nsd==3
      wzfg=Nf*wzf;
    end
    axfg=Nf*axf;
    ayfg=Nf*ayf;
    if nsd==3
      azfg=Nf*azf;
    end
    Rfg=Nf*Rf;
    Wxfg=Nf*Wxf;
    Wyfg=Nf*Wyf;
    if nsd==3
      Wzfg=Nf*Wzf;
    end
    if isDirichlet_r
      rDfg=rD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isDirichlet_w
      wDfg=wD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      wDxfg=wDfg(:,1);
      wDyfg=wDfg(:,2);
      if nsd==3
        wDzfg=wDfg(:,3);
      end
    end
    if isNeumann_t
      tNfg=tN(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
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
    Mfnx=NwfT*(nx.*Nf);
    Mfny=NwfT*(ny.*Nf);
    if nsd==3
      Mfnz=NwfT*(nz.*Nf);
    end
    
    % Get quantities for coupling ------------------------------------------------------------------
    if isFSI && isInterface
      % Get general parameters
      iD2=find(contains({Parameters.Problem},'Structural'));
      iFace2=Faces.Interface(2,iFace);
      NumElementNodes2=Sizes(iD2).NumElementNodes;
      NumElementFaces2=Sizes(iD2).NumElementFaces;
      X2e=NodesCoupled{iFace}';
      X2em=sum(X2e(1:NumElementFaces2,:),1)/NumElementFaces2;
      gamma=Parameters(iD2).NitschePenalty;
      
      % Get solution
      u2e=reshape(SolutionGlobalCoupled{iFace},[],1);
      if isTimeDependent
        u2olde=reshape(SolutionOldCoupled{iFace},[],BDFo+1);
      end
      
      % Compute weights at Gauss points
      [~,N2ex,N2ey,N2ez,~,~,pinvN2e]=mapShapeFunctions('Element',RefElement(iD2,iD2),...
                                                                 RefElement(iD2,iD2),X2e,nsd);
      
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
      if isTimeDependent
        u2oldxe=u2olde(n2e1,:);
        u2oldye=u2olde(n2e2,:);
        if nsd==3
          u2oldze=u2olde(n2e3,:);
        end
      end
      if nsd==2
        F2xxe=pinvN2e*(N2ex*u2xe+1); F2xye=pinvN2e*(N2ey*u2xe);
        F2yxe=pinvN2e*(N2ex*u2ye);   F2yye=pinvN2e*(N2ey*u2ye+1);
      elseif nsd==3
        F2xxe=pinvN2e*(N2ex*u2xe+1); F2xye=pinvN2e*(N2ey*u2xe);   F2xze=pinvN2e*(N2ez*u2xe);
        F2yxe=pinvN2e*(N2ex*u2ye);   F2yye=pinvN2e*(N2ey*u2ye+1); F2yze=pinvN2e*(N2ez*u2ye);
        F2zxe=pinvN2e*(N2ex*u2ze);   F2zye=pinvN2e*(N2ey*u2ze);   F2zze=pinvN2e*(N2ez*u2ze+1);
      end
      
      % Compute weights at Gauss points
      FaceNodes2=RefElement(iD2,iD2).FaceNodesElem;
      X2f=X2e(FaceNodes2(iFace2,:),:);
      [N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions('Face',RefElement(iD1,iD2),...
                                                           RefElement(iD1,iD2),Xf,nsd);
      N21f=RefElement(iD2,iD1).ShapeFunctionsFace;
      [~,~,~,~,w2fg]=mapShapeFunctions('Face',RefElement(iD2,iD2),...
                                              RefElement(iD2,iD2),X2f,nsd);
      
      % Compute characteristic element size
      h=sum(w2fg);
      
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
      if isTimeDependent
        u2oldxf=u2oldxe(n2f1,:);
        u2oldyf=u2oldye(n2f1,:);
        if nsd==3
          u2oldzf=u2oldze(n2f1,:);
        end
      end
      if nsd==2
        F2xxf=F2xxe(n2f1); F2xyf=F2xye(n2f1);
        F2yxf=F2yxe(n2f1); F2yyf=F2yye(n2f1);
      elseif nsd==3
        F2xxf=F2xxe(n2f1); F2xyf=F2xye(n2f1); F2xzf=F2xze(n2f1);
        F2yxf=F2yxe(n2f1); F2yyf=F2yye(n2f1); F2yzf=F2yze(n2f1);
        F2zxf=F2zxe(n2f1); F2zyf=F2zye(n2f1); F2zzf=F2zze(n2f1);
      end
      
      % Compute variables at Gauss points
      u2xfg=N21f*u2xf;
      u2yfg=N21f*u2yf;
      if nsd==3
        u2zfg=N21f*u2zf;
      end
      if isTimeDependent
        u2oldxfg=N21f*u2oldxf;
        u2oldyfg=N21f*u2oldyf;
        if nsd==3
          u2oldzfg=N21f*u2oldzf;
        end
      end
      if nsd==2
        F2xxfg=N21f*F2xxf; F2xyfg=N21f*F2xyf;
        F2yxfg=N21f*F2yxf; F2yyfg=N21f*F2yyf;
      elseif nsd==3
        F2xxfg=N21f*F2xxf; F2xyfg=N21f*F2xyf; F2xzfg=N21f*F2xzf;
        F2yxfg=N21f*F2yxf; F2yyfg=N21f*F2yyf; F2yzfg=N21f*F2yzf;
        F2zxfg=N21f*F2zxf; F2zyfg=N21f*F2zyf; F2zzfg=N21f*F2zzf;
      end
      
      % Compute stress
      if nsd==2
        [s2fg]=...
          computeStress('yes','no','push-forward',Parameters(iD2),X2em,...
          [F2xxfg,F2xyfg,F2yxfg,F2yyfg],Sizes(iD2));
        s2xxfg=s2fg(:,1); s2xyfg=s2fg(:,2);
        s2yxfg=s2fg(:,3); s2yyfg=s2fg(:,4);
      elseif nsd==3
        [s2fg]=...
          computeStress('yes','no','push-forward',Parameters(iD2),X2em,...
          [F2xxfg,F2xyfg,F2xzfg,F2yxfg,F2yyfg,F2yzfg,F2zxfg,F2zyfg,F2zzfg],Sizes(iD2));
        s2xxfg=s2fg(:,1); s2xyfg=s2fg(:,2); s2xzfg=s2fg(:,3);
        s2yxfg=s2fg(:,4); s2yyfg=s2fg(:,5); s2yzfg=s2fg(:,6);
        s2zxfg=s2fg(:,7); s2zyfg=s2fg(:,8); s2zzfg=s2fg(:,9);
      end
      
      % Compute basic matrices
      Nw12fT=(w12fg.*N12f)';
    end
    % ----------------------------------------------------------------------------------------------
    
    % Compute common terms
    if isDirichlet_r
      Rfg=rDfg;
    end
    if isDirichlet_w_x
      Wxfg=wDxfg;
    end
    if isDirichlet_w_y
      Wyfg=wDyfg;
    end
    if nsd==3 && isDirichlet_w_z
      Wzfg=wDzfg;
    end
    pfg=p(Rfg);
    dpdRfg=dpdr(Rfg);
    
    % Compute lhs
    Krr(nf1,nf1)=Krr(nf1,nf1)+tauR*Mf;
    
    Kww(nf1,nf1)=Kww(nf1,nf1)+tauW*Mf;
    Kww(nf2,nf2)=Kww(nf2,nf2)+tauW*Mf;
    if nsd==3
      Kww(nf3,nf3)=Kww(nf3,nf3)+tauW*Mf;
    end
    
    if not(isDirichlet_r)
      if nsd==2
        KLR(nf1,nefR1)=KLR(nf1,nefR1)+Voigt1*NwfT*((Wxfg./Rfg.^2.*nx).*Nf)...
                                     +Voigt2*NwfT*((Wyfg./Rfg.^2.*ny).*Nf);
        KLR(nf2,nefR1)=KLR(nf2,nefR1)+Voigt2*NwfT*((Wxfg./Rfg.^2.*nx).*Nf)...
                                     +Voigt1*NwfT*((Wyfg./Rfg.^2.*ny).*Nf);
        KLR(nf3,nefR1)=KLR(nf3,nefR1)+Voigt3*NwfT*((Wxfg./Rfg.^2.*ny).*Nf)...
                                     +Voigt3*NwfT*((Wyfg./Rfg.^2.*nx).*Nf);
      elseif nsd==3
        KLR(nf1,nefR1)=KLR(nf1,nefR1)+Voigt1*NwfT*((Wxfg./Rfg.^2.*nx).*Nf)...
                                     +Voigt2*NwfT*((Wyfg./Rfg.^2.*ny).*Nf)...
                                     +Voigt2*NwfT*((Wzfg./Rfg.^2.*nz).*Nf);
        KLR(nf2,nefR1)=KLR(nf2,nefR1)+Voigt2*NwfT*((Wxfg./Rfg.^2.*nx).*Nf)...
                                     +Voigt1*NwfT*((Wyfg./Rfg.^2.*ny).*Nf)...
                                     +Voigt2*NwfT*((Wzfg./Rfg.^2.*nz).*Nf);
        KLR(nf3,nefR1)=KLR(nf3,nefR1)+Voigt2*NwfT*((Wxfg./Rfg.^2.*nx).*Nf)...
                                     +Voigt2*NwfT*((Wyfg./Rfg.^2.*ny).*Nf)...
                                     +Voigt1*NwfT*((Wzfg./Rfg.^2.*nz).*Nf);
        KLR(nf4,nefR1)=KLR(nf4,nefR1)+Voigt3*NwfT*((Wxfg./Rfg.^2.*ny).*Nf)...
                                     +Voigt3*NwfT*((Wyfg./Rfg.^2.*nx).*Nf);
        KLR(nf5,nefR1)=KLR(nf5,nefR1)+Voigt3*NwfT*((Wxfg./Rfg.^2.*nz).*Nf)...
                                     +Voigt3*NwfT*((Wzfg./Rfg.^2.*nx).*Nf);
        KLR(nf6,nefR1)=KLR(nf6,nefR1)+Voigt3*NwfT*((Wyfg./Rfg.^2.*nz).*Nf)...
                                     +Voigt3*NwfT*((Wzfg./Rfg.^2.*ny).*Nf);
      end
    end
    
    if not(isDirichlet_w_x)
      if nsd==2
        KLW(nf1,nefW1)=KLW(nf1,nefW1)-Voigt1*NwfT*((1./Rfg.*nx).*Nf);
        KLW(nf2,nefW1)=KLW(nf2,nefW1)-Voigt2*NwfT*((1./Rfg.*nx).*Nf);
        KLW(nf3,nefW1)=KLW(nf3,nefW1)-Voigt3*NwfT*((1./Rfg.*ny).*Nf);
      elseif nsd==3
        KLW(nf1,nefW1)=KLW(nf1,nefW1)-Voigt1*NwfT*((1./Rfg.*nx).*Nf);
        KLW(nf2,nefW1)=KLW(nf2,nefW1)-Voigt2*NwfT*((1./Rfg.*nx).*Nf);
        KLW(nf3,nefW1)=KLW(nf3,nefW1)-Voigt2*NwfT*((1./Rfg.*nx).*Nf);
        KLW(nf4,nefW1)=KLW(nf4,nefW1)-Voigt3*NwfT*((1./Rfg.*ny).*Nf);
        KLW(nf5,nefW1)=KLW(nf5,nefW1)-Voigt3*NwfT*((1./Rfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet_w_y)
      if nsd==2
        KLW(nf1,nefW2)=KLW(nf1,nefW2)-Voigt2*NwfT*((1./Rfg.*ny).*Nf);
        KLW(nf2,nefW2)=KLW(nf2,nefW2)-Voigt1*NwfT*((1./Rfg.*ny).*Nf);
        KLW(nf3,nefW2)=KLW(nf3,nefW2)-Voigt3*NwfT*((1./Rfg.*nx).*Nf);
      elseif nsd==3
        KLW(nf1,nefW2)=KLW(nf1,nefW2)-Voigt2*NwfT*((1./Rfg.*ny).*Nf);
        KLW(nf2,nefW2)=KLW(nf2,nefW2)-Voigt1*NwfT*((1./Rfg.*ny).*Nf);
        KLW(nf3,nefW2)=KLW(nf3,nefW2)-Voigt2*NwfT*((1./Rfg.*ny).*Nf);
        KLW(nf4,nefW2)=KLW(nf4,nefW2)-Voigt3*NwfT*((1./Rfg.*nx).*Nf);
        KLW(nf6,nefW2)=KLW(nf6,nefW2)-Voigt3*NwfT*((1./Rfg.*nz).*Nf);
      end
    end
    
    if nsd==3 && not(isDirichlet_w_z)
      KLW(nf1,nefW3)=KLW(nf1,nefW3)-Voigt2*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf2,nefW3)=KLW(nf2,nefW3)-Voigt2*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf3,nefW3)=KLW(nf3,nefW3)-Voigt1*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf5,nefW3)=KLW(nf5,nefW3)-Voigt3*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf6,nefW3)=KLW(nf6,nefW3)-Voigt3*NwfT*((1./Rfg.*ny).*Nf);
    end
    
    if not(isDirichlet_r)
      if nsd==2
        KrR(nf1,nefR1)=KrR(nf1,nefR1)-NwfT*((axfg.*nx+ayfg.*ny).*Nf);
      elseif nsd==3
        KrR(nf1,nefR1)=KrR(nf1,nefR1)-NwfT*((axfg.*nx+ayfg.*ny+azfg.*nz).*Nf);
      end
    end
      
    if not(isDirichlet_r)
      KrR(nf1,nefR1)=KrR(nf1,nefR1)-tauR*Mf;
    end
    
    if not(isDirichlet_w_x)
      KrW(nf1,nefW1)=KrW(nf1,nefW1)+Mfnx;
    end
    
    if not(isDirichlet_w_y)
      KrW(nf1,nefW2)=KrW(nf1,nefW2)+Mfny;
    end
    
    if nsd==3 && not(isDirichlet_w_z)
        KrW(nf1,nefW3)=KrW(nf1,nefW3)+Mfnz;
    end
    
    if not(isDirichlet_r) && isConvectiveFlow
      if nsd==2
        KwR(nf1,nefR1)=KwR(nf1,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Wxfg).*Nf);
        KwR(nf2,nefR1)=KwR(nf2,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Wyfg).*Nf);
      elseif nsd==3
        KwR(nf1,nefR1)=KwR(nf1,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wxfg).*Nf);
        KwR(nf2,nefR1)=KwR(nf2,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wyfg).*Nf);
        KwR(nf3,nefR1)=KwR(nf3,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wzfg).*Nf);
      end
    end
    
    if not(isDirichlet_w_x) && isConvectiveFlow
      if nsd==2
        KwW(nf1,nefW1)=KwW(nf1,nefW1)+NwfT*(((2*Wxfg./Rfg-axfg).*nx...
                                            +(  Wyfg./Rfg-ayfg).*ny).*Nf);
        KwW(nf2,nefW1)=KwW(nf2,nefW1)+NwfT*((Wyfg./Rfg.*nx).*Nf);
      elseif nsd==3
        KwW(nf1,nefW1)=KwW(nf1,nefW1)+NwfT*(((2*Wxfg./Rfg-axfg).*nx...
                                            +(  Wyfg./Rfg-ayfg).*ny...
                                            +(  Wzfg./Rfg-azfg).*nz).*Nf);
        KwW(nf2,nefW1)=KwW(nf2,nefW1)+NwfT*((Wyfg./Rfg.*nx).*Nf);
        KwW(nf3,nefW1)=KwW(nf3,nefW1)+NwfT*((Wzfg./Rfg.*nx).*Nf);
      end
    end
    
    if not(isDirichlet_w_y) && isConvectiveFlow
      if nsd==2
        KwW(nf1,nefW2)=KwW(nf1,nefW2)+NwfT*((Wxfg./Rfg.*ny).*Nf);
        KwW(nf2,nefW2)=KwW(nf2,nefW2)+NwfT*(((  Wxfg./Rfg-axfg).*nx...
                                            +(2*Wyfg./Rfg-ayfg).*ny).*Nf);
      elseif nsd==3
        KwW(nf1,nefW2)=KwW(nf1,nefW2)+NwfT*((Wxfg./Rfg.*ny).*Nf);
        KwW(nf2,nefW2)=KwW(nf2,nefW2)+NwfT*(((  Wxfg./Rfg-axfg).*nx...
                                            +(2*Wyfg./Rfg-ayfg).*ny...
                                            +(  Wzfg./Rfg-azfg).*nz).*Nf);
        KwW(nf3,nefW2)=KwW(nf3,nefW2)+NwfT*((Wzfg./Rfg.*ny).*Nf);
      end
    end
    
    if nsd==3 && not(isDirichlet_w_z) && isConvectiveFlow
      KwW(nf1,nefW3)=KwW(nf1,nefW3)+NwfT*((Wxfg./Rfg.*nz).*Nf);
      KwW(nf2,nefW3)=KwW(nf2,nefW3)+NwfT*((Wyfg./Rfg.*nz).*Nf);
      KwW(nf3,nefW3)=KwW(nf3,nefW3)+NwfT*(((  Wxfg./Rfg-axfg).*nx...
                                          +(  Wyfg./Rfg-ayfg).*ny...
                                          +(2*Wzfg./Rfg-azfg).*nz).*Nf);
    end
    
    if not(isDirichlet_w_x)
      KwW(nf1,nefW1)=KwW(nf1,nefW1)-tauW*Mf;
    end
    
    if not(isDirichlet_w_y)
      KwW(nf2,nefW2)=KwW(nf2,nefW2)-tauW*Mf;
    end
    
    if nsd==3 && not(isDirichlet_w_z)
      KwW(nf3,nefW3)=KwW(nf3,nefW3)-tauW*Mf;
    end
    
    if not(isDirichlet_r)
      KRr(nefR1,nf1)=KRr(nefR1,nf1)+tauR*Mf;
    end
    
    if not(isDirichlet_r)
      KRR(nefR1,nefR1)=KRR(nefR1,nefR1)-tauR*Mf;
    end
    
    if not(isExterior) || isNeumann_t_x || (isFSI && isInterface)
      if nsd==2
        KWL(nefW1,nf1)=KWL(nefW1,nf1)-Voigt1*Mfnx;
        KWL(nefW1,nf2)=KWL(nefW1,nf2)-Voigt2*Mfnx;
        KWL(nefW1,nf3)=KWL(nefW1,nf3)-Voigt3*Mfny;
      elseif nsd==3
        KWL(nefW1,nf1)=KWL(nefW1,nf1)-Voigt1*Mfnx;
        KWL(nefW1,nf2)=KWL(nefW1,nf2)-Voigt2*Mfnx;
        KWL(nefW1,nf3)=KWL(nefW1,nf3)-Voigt2*Mfnx;
        KWL(nefW1,nf4)=KWL(nefW1,nf4)-Voigt3*Mfny;
        KWL(nefW1,nf5)=KWL(nefW1,nf5)-Voigt3*Mfnz;
      end
    end
    
    if not(isExterior) || isNeumann_t_y || (isFSI && isInterface)
      if nsd==2
        KWL(nefW2,nf1)=KWL(nefW2,nf1)-Voigt2*Mfny;
        KWL(nefW2,nf2)=KWL(nefW2,nf2)-Voigt1*Mfny;
        KWL(nefW2,nf3)=KWL(nefW2,nf3)-Voigt3*Mfnx;
      elseif nsd==3
        KWL(nefW2,nf1)=KWL(nefW2,nf1)-Voigt2*Mfny;
        KWL(nefW2,nf2)=KWL(nefW2,nf2)-Voigt1*Mfny;
        KWL(nefW2,nf3)=KWL(nefW2,nf3)-Voigt2*Mfny;
        KWL(nefW2,nf4)=KWL(nefW2,nf4)-Voigt3*Mfnx;
        KWL(nefW2,nf6)=KWL(nefW2,nf6)-Voigt3*Mfnz;
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z || (isFSI && isInterface))
      KWL(nefW3,nf1)=KWL(nefW3,nf1)-Voigt2*Mfnz;
      KWL(nefW3,nf2)=KWL(nefW3,nf2)-Voigt2*Mfnz;
      KWL(nefW3,nf3)=KWL(nefW3,nf3)-Voigt1*Mfnz;
      KWL(nefW3,nf5)=KWL(nefW3,nf5)-Voigt3*Mfnx;
      KWL(nefW3,nf6)=KWL(nefW3,nf6)-Voigt3*Mfny;
    end
    
    if (not(isExterior) || isNeumann_t_x || (isFSI && isInterface)) && not(isDirichlet_r)
      KWR(nefW1,nefR1)=KWR(nefW1,nefR1)-NwfT*((dpdRfg.*nx).*Nf);
    end
    
    if (not(isExterior) || isNeumann_t_y || (isFSI && isInterface)) && not(isDirichlet_r)
      KWR(nefW2,nefR1)=KWR(nefW2,nefR1)-NwfT*((dpdRfg.*ny).*Nf);
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z || (isFSI && isInterface)) && not(isDirichlet_r)
      KWR(nefW3,nefR1)=KWR(nefW3,nefR1)-NwfT*((dpdRfg.*nz).*Nf);
    end
    
    if not(isDirichlet_w_x)
      KWw(nefW1,nf1)=KWw(nefW1,nf1)-tauW*Mf;
    end
    
    if not(isDirichlet_w_y)
      KWw(nefW2,nf2)=KWw(nefW2,nf2)-tauW*Mf;
    end
    
    if nsd==3 && not(isDirichlet_w_z)
      KWw(nefW3,nf3)=KWw(nefW3,nf3)-tauW*Mf;
    end
    
    if not(isDirichlet_w_x)
      KWW(nefW1,nefW1)=KWW(nefW1,nefW1)+tauW*Mf;
    end
    
    if not(isDirichlet_w_y)
      KWW(nefW2,nefW2)=KWW(nefW2,nefW2)+tauW*Mf;
    end
    
    if nsd==3 && not(isDirichlet_w_z)
      KWW(nefW3,nefW3)=KWW(nefW3,nefW3)+tauW*Mf;
    end
    
    if isFSI && isInterface
      KWR(nefW1,nefR1)=KWR(nefW1,nefR1)-NwfT*((gamma/h*Wxfg./Rfg.^2).*Nf);
      KWR(nefW2,nefR1)=KWR(nefW2,nefR1)-NwfT*((gamma/h*Wyfg./Rfg.^2).*Nf);
      if nsd==3
        KWR(nefW3,nefR1)=KWR(nefW3,nefR1)-NwfT*((gamma/h*Wzfg./Rfg.^2).*Nf);
      end
      
      KWW(nefW1,nefW1)=KWW(nefW1,nefW1)+NwfT*((gamma/h*1./Rfg).*Nf);
      KWW(nefW2,nefW2)=KWW(nefW2,nefW2)+NwfT*((gamma/h*1./Rfg).*Nf);
      if nsd==3
        KWW(nefW3,nefW3)=KWW(nefW3,nefW3)+NwfT*((gamma/h*1./Rfg).*Nf);
      end
    end
    
    % Compute rhs
    fr(nf1,1)=fr(nf1,1)-NwfT*(tauR*rfg);
    
    fw(nf1,1)=fw(nf1,1)-NwfT*(tauW*wxfg);
    fw(nf2,1)=fw(nf2,1)-NwfT*(tauW*wyfg);
    if nsd==3
      fw(nf3,1)=fw(nf3,1)-NwfT*(tauW*wzfg);
    end
    
    if nsd==2
      fL(nf1,1)=fL(nf1,1)+NwfT*(+Voigt1*nx.*Wxfg./Rfg...
                                +Voigt2*ny.*Wyfg./Rfg);
      fL(nf2,1)=fL(nf2,1)+NwfT*(+Voigt2*nx.*Wxfg./Rfg...
                                +Voigt1*ny.*Wyfg./Rfg);
      fL(nf3,1)=fL(nf3,1)+NwfT*(+Voigt3*ny.*Wxfg./Rfg...
                                +Voigt3*nx.*Wyfg./Rfg);
    elseif nsd==3
      fL(nf1,1)=fL(nf1,1)+NwfT*(+Voigt1*nx.*Wxfg./Rfg...
                                +Voigt2*ny.*Wyfg./Rfg...
                                +Voigt2*nz.*Wzfg./Rfg);
      fL(nf2,1)=fL(nf2,1)+NwfT*(+Voigt2*nx.*Wxfg./Rfg...
                                +Voigt1*ny.*Wyfg./Rfg...
                                +Voigt2*nz.*Wzfg./Rfg);
      fL(nf3,1)=fL(nf3,1)+NwfT*(+Voigt2*nx.*Wxfg./Rfg...
                                +Voigt2*ny.*Wyfg./Rfg...
                                +Voigt1*nz.*Wzfg./Rfg);
      fL(nf4,1)=fL(nf4,1)+NwfT*(+Voigt3*nx.*Wyfg./Rfg...
                                +Voigt3*ny.*Wxfg./Rfg);
      fL(nf5,1)=fL(nf5,1)+NwfT*(+Voigt3*nx.*Wzfg./Rfg...
                                +Voigt3*nz.*Wxfg./Rfg);
      fL(nf6,1)=fL(nf6,1)+NwfT*(+Voigt3*ny.*Wzfg./Rfg...
                                +Voigt3*nz.*Wyfg./Rfg);
    end
    
    fr(nf1,1)=fr(nf1,1)-NwfT*((Wxfg-Rfg.*axfg).*nx+(Wyfg-Rfg.*ayfg).*ny-tauR*Rfg);
    if nsd==3
      fr(nf1,1)=fr(nf1,1)-NwfT*((Wzfg-Rfg.*azfg).*nz);
    end
    
    if isConvectiveFlow
      if nsd==2
        fw(nf1,1)=fw(nf1,1)-NwfT*(Wxfg.*((Wxfg./Rfg-axfg).*nx...
                                        +(Wyfg./Rfg-ayfg).*ny));
        fw(nf2,1)=fw(nf2,1)-NwfT*(Wyfg.*((Wxfg./Rfg-axfg).*nx...
                                        +(Wyfg./Rfg-ayfg).*ny));
      elseif nsd==3
        fw(nf1,1)=fw(nf1,1)-NwfT*(Wxfg.*((Wxfg./Rfg-axfg).*nx...
                                        +(Wyfg./Rfg-ayfg).*ny...
                                        +(Wzfg./Rfg-azfg).*nz));
        fw(nf2,1)=fw(nf2,1)-NwfT*(Wyfg.*((Wxfg./Rfg-axfg).*nx...
                                        +(Wyfg./Rfg-ayfg).*ny...
                                        +(Wzfg./Rfg-azfg).*nz));
        fw(nf3,1)=fw(nf3,1)-NwfT*(Wzfg.*((Wxfg./Rfg-axfg).*nx...
                                        +(Wyfg./Rfg-ayfg).*ny...
                                        +(Wzfg./Rfg-azfg).*nz));
      end
    end
    
    fw(nf1,1)=fw(nf1,1)+NwfT*(tauW*Wxfg);
    fw(nf2,1)=fw(nf2,1)+NwfT*(tauW*Wyfg);
    if nsd==3
      fw(nf3,1)=fw(nf3,1)+NwfT*(tauW*Wzfg);
    end
    
    if not(isDirichlet_r)
      fR(nefR1,1)=fR(nefR1,1)-NwfT*(tauR*(rfg-Rfg));
    end
    
    if not(isExterior) || isNeumann_t_x || (isFSI && isInterface)
      if nsd==2
        fW(nefW1,1)=fW(nefW1,1)+NwfT*(+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg...
                                      +Voigt3*ny.*Lxyfg...
                                      +pfg.*nx);
      elseif nsd==3
        fW(nefW1,1)=fW(nefW1,1)+NwfT*(+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg+Voigt2*nx.*Lzzfg...
                                      +Voigt3*ny.*Lxyfg+Voigt3*nz.*Lxzfg...
                                      +pfg.*nx);
      end
    end
    
    if not(isExterior) || isNeumann_t_y || (isFSI && isInterface)
      if nsd==2
        fW(nefW2,1)=fW(nefW2,1)+NwfT*(+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg...
                                      +Voigt3*nx.*Lxyfg...
                                      +pfg.*ny);
      elseif nsd==3
        fW(nefW2,1)=fW(nefW2,1)+NwfT*(+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg+Voigt2*ny.*Lzzfg...
                                      +Voigt3*nx.*Lxyfg+Voigt3*nz.*Lyzfg...
                                      +pfg.*ny);
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z || (isFSI && isInterface))
      fW(nefW3,1)=fW(nefW3,1)+NwfT*(+Voigt2*nz.*Lxxfg+Voigt2*nz.*Lyyfg+Voigt1*nz.*Lzzfg...
                                    +Voigt3*nx.*Lxzfg+Voigt3*ny.*Lyzfg...
                                    +pfg.*nz);
    end
    
    if not(isDirichlet_w_x)
      fW(nefW1,1)=fW(nefW1,1)+NwfT*(+tauW*(wxfg-Wxfg));
    end
    
    if not(isDirichlet_w_y)
      fW(nefW2,1)=fW(nefW2,1)+NwfT*(+tauW*(wyfg-Wyfg));
    end
    
    if nsd==3 && not(isDirichlet_w_z)
      fW(nefW3,1)=fW(nefW3,1)+NwfT*(+tauW*(wzfg-Wzfg));
    end
    
    if isFSI && isInterface
      if nsd==2
        fW(nefW1,1)=fW(nefW1,1)-Nw12fT*(+s2xxfg.*(-n12x)+s2xyfg.*(-n12y));
        fW(nefW2,1)=fW(nefW2,1)-Nw12fT*(+s2yxfg.*(-n12x)+s2yyfg.*(-n12y));
      elseif nsd==3
        fW(nefW1,1)=fW(nefW1,1)-Nw12fT*(+s2xxfg.*(-n12x)+s2xyfg.*(-n12y)+s2xzfg.*(-n12z));
        fW(nefW2,1)=fW(nefW2,1)-Nw12fT*(+s2yxfg.*(-n12x)+s2yyfg.*(-n12y)+s2yzfg.*(-n12z));
        fW(nefW3,1)=fW(nefW3,1)-Nw12fT*(+s2zxfg.*(-n12x)+s2zyfg.*(-n12y)+s2zzfg.*(-n12z));
      end
      
      if isTimeDependent
        fW(nefW1,1)=fW(nefW1,1)+Nw12fT*(gamma/h*(1/dt*u2xfg*alpha(1)...
                                                +1/dt*u2oldxfg(:,1:BDFo)*alpha(2:BDFo+1,1)));
        fW(nefW2,1)=fW(nefW2,1)+Nw12fT*(gamma/h*(1/dt*u2yfg*alpha(1)...
                                                +1/dt*u2oldyfg(:,1:BDFo)*alpha(2:BDFo+1,1)));
        if nsd==3
          fW(nefW3,1)=fW(nefW3,1)+Nw12fT*(gamma/h*(1/dt*u2zfg*alpha(1)...
                                                  +1/dt*u2oldzfg(:,1:BDFo)*alpha(2:BDFo+1,1)));
        end
      end
      
      fW(nefW1,1)=fW(nefW1,1)-NwfT*(gamma/h*Wxfg./Rfg);
      fW(nefW2,1)=fW(nefW2,1)-NwfT*(gamma/h*Wyfg./Rfg);
      if nsd==3
        fW(nefW3,1)=fW(nefW3,1)-NwfT*(gamma/h*Wzfg./Rfg);
      end
    end
    
    if isNeumann_t_x
      fW(nefW1,1)=fW(nefW1,1)+NwfT*(tNxfg);
    end
    
    if isNeumann_t_y
      fW(nefW2,1)=fW(nefW2,1)+NwfT*(tNyfg);
    end
    
    if nsd==3 && isNeumann_t_z
      fW(nefW3,1)=fW(nefW3,1)+NwfT*(tNzfg);
    end
    
    % Remove undetermination
    if isDirichlet_r
      KRR(nefR1,nefR1)=eye(NumFaceNodes);
    end
    if isDirichlet_w_x
      KWW(nefW1,nefW1)=eye(NumFaceNodes);
    end
    if isDirichlet_w_y
      KWW(nefW2,nefW2)=eye(NumFaceNodes);
    end
    if nsd==3 && isDirichlet_w_z
      KWW(nefW3,nefW3)=eye(NumFaceNodes);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iL=1:msd*NumElementNodes;
ir=iL(end)+(1:NumElementNodes);
iw=ir(end)+(1:nsd*NumElementNodes);
iR=reshape((0:NumElementFaces-1)*(1+nsd)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iW=reshape((0:NumElementFaces-1)*(1+nsd)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[])+NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+1+nsd)*NumElementNodes,(msd+1+nsd)*NumElementNodes);
LhsLG=zeros((msd+1+nsd)*NumElementNodes,(1+nsd)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((1+nsd)*NumElementFaces*NumFaceNodes,(msd+1+nsd)*NumElementNodes);
LhsGG=zeros((1+nsd)*NumElementFaces*NumFaceNodes,(1+nsd)*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+1+nsd)*NumElementNodes,1);
RhsG=zeros((1+nsd)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iL,iL)=KLL;
LhsLL(iL,ir)=KLr;
LhsLL(iL,iw)=KLw;
LhsLL(ir,ir)=Krr;
LhsLL(ir,iw)=Krw;
LhsLL(iw,iL)=KwL;
LhsLL(iw,ir)=Kwr;
LhsLL(iw,iw)=Kww;

% Lhs local-global
LhsLG(iL,iR)=KLR;
LhsLG(iL,iW)=KLW;
LhsLG(ir,iR)=KrR;
LhsLG(ir,iW)=KrW;
LhsLG(iw,iR)=KwR;
LhsLG(iw,iW)=KwW;

% Rhs local
RhsL(iL,1)=fL;
RhsL(ir,1)=fr;
RhsL(iw,1)=fw;

% Lhs global-local
LhsGL(iR,ir)=KRr;
LhsGL(iW,iL)=KWL;
LhsGL(iW,iw)=KWw;

% Lhs global-global
LhsGG(iR,iR)=KRR;
LhsGG(iW,iR)=KWR;
LhsGG(iW,iW)=KWW;

% Rhs global
RhsG(iR,1)=fR;
RhsG(iW,1)=fW;

% Matrix and vector for local problem
MatVecLocal=LhsLL\[LhsLG,RhsL];

% Extract matrix for local problem
MatLocal=MatVecLocal(:,1:end-1);

% Extract vector for local problem
VecLocal=MatVecLocal(:,end);

% Lhs for global problem
LhsGlobal=LhsGG-LhsGL*MatLocal;

% Rhs for global problem
RhsGlobal=RhsG-LhsGL*VecLocal;

end

%% Do coupling element
function [LhsCoup]=doCouplingElement(...
  iFaceInterface,iD1,iD2,Block,Parameters,Mesh,Faces,Time,RefElement,Sizes)

% Get general parameters
isArbitraryLagrangianEulerian=strcmp(Parameters(iD1).ArbitraryLagrangianEulerian,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
iElem1=Faces(iD1,iD2).Interface(iFaceInterface,1);
iElem2=Faces(iD1,iD2).Interface(iFaceInterface,3);
iFace1=Faces(iD1,iD2).Interface(iFaceInterface,2);
iFace2=Faces(iD1,iD2).Interface(iFaceInterface,4);
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes2=Sizes(iD2).NumElementNodes;
NumElementFaces1=Sizes(iD1).NumElementFaces;
NumElementFaces2=Sizes(iD2).NumElementFaces;
NumFaceNodes1=Sizes(iD1).NumFaceNodes;
C1e=Mesh(iD1).Elements(:,iElem1)';
C2e=Mesh(iD2).Elements(:,iElem2)';
X1e=Mesh(iD1).Nodes(:,C1e)';
X10e=Mesh(iD1).NodesInitial(:,C1e)';
X2e=Mesh(iD2).Nodes(:,C2e)';
X2em=sum(X2e(1:NumElementFaces2,:),1)/NumElementFaces2;
gamma=Parameters(iD2).NitschePenalty;
if isTimeDependent
  dt=Time.TimeStepSize;
  alpha=Time.BDF1stDerEff;
end

% Get solution
u2e=reshape(Block(iD2,iD2).SolutionGlobal(C2e,:),[],1);
if isArbitraryLagrangianEulerian
  iD3=find(contains({Parameters.Problem},'Mesh'));
  C3e=Mesh(iD3).Elements(:,iElem1)';
  u3e=reshape(Block(iD3,iD3).SolutionGlobal(C3e,:),[],1);
end

% Initialize lhs
KW1u2=zeros(nsd*NumElementFaces1*NumFaceNodes1,nsd*NumElementNodes2);

% Update nodes coordinates
if isArbitraryLagrangianEulerian
  N310e=RefElement(iD3,iD1).ShapeFunctionsElem;
  pinvN103e=RefElement(iD1,iD3).PseudoinverseShapeFunctionsElem;
  X1e(:,1:nsd)=X10e(:,1:nsd)+pinvN103e*(N310e*reshape(u3e,[],nsd));
end

% Compute weights at Gauss points
[~,N2ex,N2ey,N2ez,~,~,pinvN2e]=mapShapeFunctions('Element',RefElement(iD2,iD2),...
                                                           RefElement(iD2,iD2),X2e,nsd);

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
if nsd==2
  F2xxe=pinvN2e*(N2ex*u2xe+1); F2xye=pinvN2e*(N2ey*u2xe);
  F2yxe=pinvN2e*(N2ex*u2ye);   F2yye=pinvN2e*(N2ey*u2ye+1);
elseif nsd==3
  F2xxe=pinvN2e*(N2ex*u2xe+1); F2xye=pinvN2e*(N2ey*u2xe);   F2xze=pinvN2e*(N2ez*u2xe);
  F2yxe=pinvN2e*(N2ex*u2ye);   F2yye=pinvN2e*(N2ey*u2ye+1); F2yze=pinvN2e*(N2ez*u2ye);
  F2zxe=pinvN2e*(N2ex*u2ze);   F2zye=pinvN2e*(N2ey*u2ze);   F2zze=pinvN2e*(N2ez*u2ze+1);
end

% Compute weights at Gauss points
FaceNodes1=RefElement(iD1,iD1).FaceNodesElem;
FaceNodes2=RefElement(iD2,iD2).FaceNodesElem;
X1f=X1e(FaceNodes1(iFace1,:),:);
X2f=X2e(FaceNodes2(iFace2,:),:);
[N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions('Face',RefElement(iD1,iD2),...
                                                     RefElement(iD1,iD2),X1f,nsd);
N21f=RefElement(iD2,iD1).ShapeFunctionsFace;
[~,~,~,~,w2fg]=mapShapeFunctions('Face',RefElement(iD2,iD2),...
                                        RefElement(iD2,iD2),X2f,nsd);

% Compute characteristic element size
h=sum(w2fg);

% Indices
n1efW1=(iFace1-1)*nsd*NumFaceNodes1+(1:NumFaceNodes1);
n1efW2=n1efW1+NumFaceNodes1;
n1efW3=n1efW2+NumFaceNodes1;
n2f1=FaceNodes2(iFace2,:);
n2f2=n2f1+NumElementNodes2;
n2f3=n2f2+NumElementNodes2;

% Flip face
Node2Match1stNode1=Faces(iD1,iD2).Interface(iFaceInterface,5);
order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
n2f1=n2f1(order);
n2f2=n2f2(order);
n2f3=n2f3(order);

% Compute derivatives of shape functions
N2xe=pinvN2e*N2ex;
N2ye=pinvN2e*N2ey;
if nsd==3
  N2ze=pinvN2e*N2ez;
end
N21xf=N21f*N2xe(n2f1,:);
N21yf=N21f*N2ye(n2f1,:);
if nsd==3
  N21zf=N21f*N2ze(n2f1,:);
end

% Compute variables at nodes
if nsd==2
  F2xxf=F2xxe(n2f1); F2xyf=F2xye(n2f1);
  F2yxf=F2yxe(n2f1); F2yyf=F2yye(n2f1);
elseif nsd==3
  F2xxf=F2xxe(n2f1); F2xyf=F2xye(n2f1); F2xzf=F2xze(n2f1);
  F2yxf=F2yxe(n2f1); F2yyf=F2yye(n2f1); F2yzf=F2yze(n2f1);
  F2zxf=F2zxe(n2f1); F2zyf=F2zye(n2f1); F2zzf=F2zze(n2f1);
end

% Compute variables at Gauss points
if nsd==2
  F2xxfg=N21f*F2xxf; F2xyfg=N21f*F2xyf;
  F2yxfg=N21f*F2yxf; F2yyfg=N21f*F2yyf;
elseif nsd==3
  F2xxfg=N21f*F2xxf; F2xyfg=N21f*F2xyf; F2xzfg=N21f*F2xzf;
  F2yxfg=N21f*F2yxf; F2yyfg=N21f*F2yyf; F2yzfg=N21f*F2yzf;
  F2zxfg=N21f*F2zxf; F2zyfg=N21f*F2zyf; F2zzfg=N21f*F2zzf;
end

% Compute basic matrices
Nw12fT=(w12fg.*N12f)';
M12f=Nw12fT*N21f;

% Compute stress linearization
if nsd==2
  [~,ds2dF2fg]=...
    computeStress('no','yes','push-forward',Parameters(iD2),X2em,...
    [F2xxfg,F2xyfg,F2yxfg,F2yyfg],Sizes(iD2));
  ds2xxdF2xxfg=ds2dF2fg(:,1);  ds2xxdF2xyfg=ds2dF2fg(:,2);  ds2xxdF2yxfg=ds2dF2fg(:,3);  ds2xxdF2yyfg=ds2dF2fg(:,4);
  ds2xydF2xxfg=ds2dF2fg(:,5);  ds2xydF2xyfg=ds2dF2fg(:,6);  ds2xydF2yxfg=ds2dF2fg(:,7);  ds2xydF2yyfg=ds2dF2fg(:,8);
  ds2yxdF2xxfg=ds2dF2fg(:,9);  ds2yxdF2xyfg=ds2dF2fg(:,10); ds2yxdF2yxfg=ds2dF2fg(:,11); ds2yxdF2yyfg=ds2dF2fg(:,12);
  ds2yydF2xxfg=ds2dF2fg(:,13); ds2yydF2xyfg=ds2dF2fg(:,14); ds2yydF2yxfg=ds2dF2fg(:,15); ds2yydF2yyfg=ds2dF2fg(:,16);
elseif nsd==3
  [~,ds2dF2fg]=...
    computeStress('no','yes','push-forward',Parameters(iD2),X2em,...
    [F2xxfg,F2xyfg,F2xzfg,F2yxfg,F2yyfg,F2yzfg,F2zxfg,F2zyfg,F2zzfg],Sizes(iD2));
  ds2xxdF2xxfg=ds2dF2fg(:,1);  ds2xxdF2xyfg=ds2dF2fg(:,2);  ds2xxdF2xzfg=ds2dF2fg(:,3);  ds2xxdF2yxfg=ds2dF2fg(:,4);  ds2xxdF2yyfg=ds2dF2fg(:,5);  ds2xxdF2yzfg=ds2dF2fg(:,6);  ds2xxdF2zxfg=ds2dF2fg(:,7);  ds2xxdF2zyfg=ds2dF2fg(:,8);  ds2xxdF2zzfg=ds2dF2fg(:,9);
  ds2xydF2xxfg=ds2dF2fg(:,10); ds2xydF2xyfg=ds2dF2fg(:,11); ds2xydF2xzfg=ds2dF2fg(:,12); ds2xydF2yxfg=ds2dF2fg(:,13); ds2xydF2yyfg=ds2dF2fg(:,14); ds2xydF2yzfg=ds2dF2fg(:,15); ds2xydF2zxfg=ds2dF2fg(:,16); ds2xydF2zyfg=ds2dF2fg(:,17); ds2xydF2zzfg=ds2dF2fg(:,18);
  ds2xzdF2xxfg=ds2dF2fg(:,19); ds2xzdF2xyfg=ds2dF2fg(:,20); ds2xzdF2xzfg=ds2dF2fg(:,21); ds2xzdF2yxfg=ds2dF2fg(:,22); ds2xzdF2yyfg=ds2dF2fg(:,23); ds2xzdF2yzfg=ds2dF2fg(:,24); ds2xzdF2zxfg=ds2dF2fg(:,25); ds2xzdF2zyfg=ds2dF2fg(:,26); ds2xzdF2zzfg=ds2dF2fg(:,27);
  ds2yxdF2xxfg=ds2dF2fg(:,28); ds2yxdF2xyfg=ds2dF2fg(:,29); ds2yxdF2xzfg=ds2dF2fg(:,30); ds2yxdF2yxfg=ds2dF2fg(:,31); ds2yxdF2yyfg=ds2dF2fg(:,32); ds2yxdF2yzfg=ds2dF2fg(:,33); ds2yxdF2zxfg=ds2dF2fg(:,34); ds2yxdF2zyfg=ds2dF2fg(:,35); ds2yxdF2zzfg=ds2dF2fg(:,36);
  ds2yydF2xxfg=ds2dF2fg(:,37); ds2yydF2xyfg=ds2dF2fg(:,38); ds2yydF2xzfg=ds2dF2fg(:,39); ds2yydF2yxfg=ds2dF2fg(:,40); ds2yydF2yyfg=ds2dF2fg(:,41); ds2yydF2yzfg=ds2dF2fg(:,42); ds2yydF2zxfg=ds2dF2fg(:,43); ds2yydF2zyfg=ds2dF2fg(:,44); ds2yydF2zzfg=ds2dF2fg(:,45);
  ds2yzdF2xxfg=ds2dF2fg(:,46); ds2yzdF2xyfg=ds2dF2fg(:,47); ds2yzdF2xzfg=ds2dF2fg(:,48); ds2yzdF2yxfg=ds2dF2fg(:,49); ds2yzdF2yyfg=ds2dF2fg(:,50); ds2yzdF2yzfg=ds2dF2fg(:,51); ds2yzdF2zxfg=ds2dF2fg(:,52); ds2yzdF2zyfg=ds2dF2fg(:,53); ds2yzdF2zzfg=ds2dF2fg(:,54);
  ds2zxdF2xxfg=ds2dF2fg(:,55); ds2zxdF2xyfg=ds2dF2fg(:,56); ds2zxdF2xzfg=ds2dF2fg(:,57); ds2zxdF2yxfg=ds2dF2fg(:,58); ds2zxdF2yyfg=ds2dF2fg(:,59); ds2zxdF2yzfg=ds2dF2fg(:,60); ds2zxdF2zxfg=ds2dF2fg(:,61); ds2zxdF2zyfg=ds2dF2fg(:,62); ds2zxdF2zzfg=ds2dF2fg(:,63);
  ds2zydF2xxfg=ds2dF2fg(:,64); ds2zydF2xyfg=ds2dF2fg(:,65); ds2zydF2xzfg=ds2dF2fg(:,66); ds2zydF2yxfg=ds2dF2fg(:,67); ds2zydF2yyfg=ds2dF2fg(:,68); ds2zydF2yzfg=ds2dF2fg(:,69); ds2zydF2zxfg=ds2dF2fg(:,70); ds2zydF2zyfg=ds2dF2fg(:,71); ds2zydF2zzfg=ds2dF2fg(:,72);
  ds2zzdF2xxfg=ds2dF2fg(:,73); ds2zzdF2xyfg=ds2dF2fg(:,74); ds2zzdF2xzfg=ds2dF2fg(:,75); ds2zzdF2yxfg=ds2dF2fg(:,76); ds2zzdF2yyfg=ds2dF2fg(:,77); ds2zzdF2yzfg=ds2dF2fg(:,78); ds2zzdF2zxfg=ds2dF2fg(:,79); ds2zzdF2zyfg=ds2dF2fg(:,80); ds2zzdF2zzfg=ds2dF2fg(:,81);
end

% Compute lhs
if nsd==2
  KW1u2(n1efW1,n2e1)=KW1u2(n1efW1,n2e1)...
            +Nw12fT*((+ds2xxdF2xxfg.*(-n12x)+ds2xydF2xxfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2xxdF2xyfg.*(-n12x)+ds2xydF2xyfg.*(-n12y)).*N21yf);
  KW1u2(n1efW1,n2e2)=KW1u2(n1efW1,n2e2)...
            +Nw12fT*((+ds2xxdF2yxfg.*(-n12x)+ds2xydF2yxfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2xxdF2yyfg.*(-n12x)+ds2xydF2yyfg.*(-n12y)).*N21yf);
  KW1u2(n1efW2,n2e1)=KW1u2(n1efW2,n2e1)...
            +Nw12fT*((+ds2yxdF2xxfg.*(-n12x)+ds2yydF2xxfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2yxdF2xyfg.*(-n12x)+ds2yydF2xyfg.*(-n12y)).*N21yf);
  KW1u2(n1efW2,n2e2)=KW1u2(n1efW2,n2e2)...
            +Nw12fT*((+ds2yxdF2yxfg.*(-n12x)+ds2yydF2yxfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2yxdF2yyfg.*(-n12x)+ds2yydF2yyfg.*(-n12y)).*N21yf);
elseif nsd==3
  KW1u2(n1efW1,n2e1)=KW1u2(n1efW1,n2e1)...
            +Nw12fT*((+ds2xxdF2xxfg.*(-n12x)+ds2xydF2xxfg.*(-n12y)+ds2xzdF2xxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xxdF2xyfg.*(-n12x)+ds2xydF2xyfg.*(-n12y)+ds2xzdF2xyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xxdF2xzfg.*(-n12x)+ds2xydF2xzfg.*(-n12y)+ds2xzdF2xzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW1,n2e2)=KW1u2(n1efW1,n2e2)...
            +Nw12fT*((+ds2xxdF2yxfg.*(-n12x)+ds2xydF2yxfg.*(-n12y)+ds2xzdF2yxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xxdF2yyfg.*(-n12x)+ds2xydF2yyfg.*(-n12y)+ds2xzdF2yyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xxdF2yzfg.*(-n12x)+ds2xydF2yzfg.*(-n12y)+ds2xzdF2yzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW1,n2e3)=KW1u2(n1efW1,n2e3)...
            +Nw12fT*((+ds2xxdF2zxfg.*(-n12x)+ds2xydF2zxfg.*(-n12y)+ds2xzdF2zxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xxdF2zyfg.*(-n12x)+ds2xydF2zyfg.*(-n12y)+ds2xzdF2zyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xxdF2zzfg.*(-n12x)+ds2xydF2zzfg.*(-n12y)+ds2xzdF2zzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW2,n2e1)=KW1u2(n1efW2,n2e1)...
            +Nw12fT*((+ds2yxdF2xxfg.*(-n12x)+ds2yydF2xxfg.*(-n12y)+ds2yzdF2xxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2yxdF2xyfg.*(-n12x)+ds2yydF2xyfg.*(-n12y)+ds2yzdF2xyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2yxdF2xzfg.*(-n12x)+ds2yydF2xzfg.*(-n12y)+ds2yzdF2xzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW2,n2e2)=KW1u2(n1efW2,n2e2)...
            +Nw12fT*((+ds2yxdF2yxfg.*(-n12x)+ds2yydF2yxfg.*(-n12y)+ds2yzdF2yxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2yxdF2yyfg.*(-n12x)+ds2yydF2yyfg.*(-n12y)+ds2yzdF2yyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2yxdF2yzfg.*(-n12x)+ds2yydF2yzfg.*(-n12y)+ds2yzdF2yzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW2,n2e3)=KW1u2(n1efW2,n2e3)...
            +Nw12fT*((+ds2yxdF2zxfg.*(-n12x)+ds2yydF2zxfg.*(-n12y)+ds2yzdF2zxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2yxdF2zyfg.*(-n12x)+ds2yydF2zyfg.*(-n12y)+ds2yzdF2zyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2yxdF2zzfg.*(-n12x)+ds2yydF2zzfg.*(-n12y)+ds2yzdF2zzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW3,n2e1)=KW1u2(n1efW3,n2e1)...
            +Nw12fT*((+ds2zxdF2xxfg.*(-n12x)+ds2zydF2xxfg.*(-n12y)+ds2zzdF2xxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2zxdF2xyfg.*(-n12x)+ds2zydF2xyfg.*(-n12y)+ds2zzdF2xyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2zxdF2xzfg.*(-n12x)+ds2zydF2xzfg.*(-n12y)+ds2zzdF2xzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW3,n2e2)=KW1u2(n1efW3,n2e2)...
            +Nw12fT*((+ds2zxdF2yxfg.*(-n12x)+ds2zydF2yxfg.*(-n12y)+ds2zzdF2yxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2zxdF2yyfg.*(-n12x)+ds2zydF2yyfg.*(-n12y)+ds2zzdF2yyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2zxdF2yzfg.*(-n12x)+ds2zydF2yzfg.*(-n12y)+ds2zzdF2yzfg.*(-n12z)).*N21zf);
  KW1u2(n1efW3,n2e3)=KW1u2(n1efW3,n2e3)...
            +Nw12fT*((+ds2zxdF2zxfg.*(-n12x)+ds2zydF2zxfg.*(-n12y)+ds2zzdF2zxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2zxdF2zyfg.*(-n12x)+ds2zydF2zyfg.*(-n12y)+ds2zzdF2zyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2zxdF2zzfg.*(-n12x)+ds2zydF2zzfg.*(-n12y)+ds2zzdF2zzfg.*(-n12z)).*N21zf);
end

if isTimeDependent
  KW1u2(n1efW1,n2f1)=KW1u2(n1efW1,n2f1)-gamma/h*alpha(1)/dt*M12f;
  KW1u2(n1efW2,n2f2)=KW1u2(n1efW2,n2f2)-gamma/h*alpha(1)/dt*M12f;
  if nsd==3
    KW1u2(n1efW3,n2f3)=KW1u2(n1efW3,n2f3)-gamma/h*alpha(1)/dt*M12f;
  end
end

% Indices
iW1=reshape((0:NumElementFaces1-1)*(1+nsd)*NumFaceNodes1+repmat((1:nsd*NumFaceNodes1)',...
  1,NumElementFaces1),1,[])+NumFaceNodes1;
iu2=1:nsd*NumElementNodes2;

% Initialization of lhs and rhs
LhsCoup=zeros((1+nsd)*NumElementFaces1*NumFaceNodes1,nsd*NumElementNodes2);

% Compute elemental contributions to lhs
LhsCoup(iW1,iu2)=KW1u2;

end

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,Parameters,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.DynamicViscosity;
lambda=-2/3*Parameters.DynamicViscosity;
rD=Parameters.Density;
wD=Parameters.Momentum;
Xe=Nodes';
t=Time.Time;

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
re=reshape(SolutionLocal(:,msd+1),[],1);
we=reshape(SolutionLocal(:,msd+1+(1:nsd)),[],1);
Ue=SolutionGlobal;

% Initialize lhs
Kpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
Ktp=zeros(nsd,nsd*NumElementNodesPost);
Krp=zeros(qsd,nsd*NumElementNodesPost);

% Initialize rhs
fp=zeros(nsd*NumElementNodesPost,1);
ft=zeros(nsd,1);
fr=zeros(qsd,1);

% Get reference data
FaceNodes=RefElement.FaceNodesElem;

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions('Element',RefElement.PostLow,RefElement.Post,Xe,nsd);
N1e=ones(length(weg),1);

% Indices
ne1=1:NumElementNodesPost;
ne2=ne1+NumElementNodesPost;
ne3=ne2+NumElementNodesPost;
nle1=1:NumElementNodes;
nle2=nle1+NumElementNodes;
nle3=nle2+NumElementNodes;
nle4=nle3+NumElementNodes;
nle5=nle4+NumElementNodes;
nle6=nle5+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Lxxe=Le(nle1);
  Lyye=Le(nle2);
  Lxye=Le(nle3);
elseif nsd==3
  Lxxe=Le(nle1);
  Lyye=Le(nle2);
  Lzze=Le(nle3);
  Lxye=Le(nle4);
  Lxze=Le(nle5);
  Lyze=Le(nle6);
end
wxe=we(nle1);
wye=we(nle2);
if nsd==3
  wze=we(nle3);
end

% Compute variables at Gauss points
if nsd==2
  Lxxeg=Nle*Lxxe;
  Lyyeg=Nle*Lyye;
  Lxyeg=Nle*Lxye;
elseif nsd==3
  Lxxeg=Nle*Lxxe;
  Lyyeg=Nle*Lyye;
  Lzzeg=Nle*Lzze;
  Lxyeg=Nle*Lxye;
  Lxzeg=Nle*Lxze;
  Lyzeg=Nle*Lyze;
end
reg=Nle*re;
wxeg=Nle*wxe;
wyeg=Nle*wye;
if nsd==3
  wzeg=Nle*wze;
end

% Compute basic matrices
Nw1eT=(weg.*N1e)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
if nsd==2
  Kxxe=NwexT*Nex; Kxye=NwexT*Ney;
  Kyxe=NweyT*Nex; Kyye=NweyT*Ney;
elseif nsd==3
  Kxxe=NwexT*Nex; Kxye=NwexT*Ney; Kxze=NwexT*Nez;
  Kyxe=NweyT*Nex; Kyye=NweyT*Ney; Kyze=NweyT*Nez;
  Kzxe=NwezT*Nex; Kzye=NwezT*Ney; Kzze=NwezT*Nez;
end

% Compute linearization of viscous stress
if nsd==2
  Voigt1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  Voigt2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  Voigt3=sqrt(mu);
elseif nsd==3
  Voigt1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  Voigt2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  Voigt3=sqrt(mu);
end

% Compute lhs
if nsd==2
  Kpp(ne1,ne1)=-Voigt1*Kxxe-Voigt3*Kyye;
  Kpp(ne1,ne2)=-Voigt2*Kxye-Voigt3*Kyxe;
  Kpp(ne2,ne1)=-Voigt2*Kyxe-Voigt3*Kxye;
  Kpp(ne2,ne2)=-Voigt1*Kyye-Voigt3*Kxxe;
elseif nsd==3
  Kpp(ne1,ne1)=-Voigt1*Kxxe-Voigt3*Kyye-Voigt3*Kzze;
  Kpp(ne1,ne2)=-Voigt2*Kxye-Voigt3*Kyxe;
  Kpp(ne1,ne3)=-Voigt2*Kxze-Voigt3*Kzxe;
  Kpp(ne2,ne1)=-Voigt2*Kyxe-Voigt3*Kxye;
  Kpp(ne2,ne2)=-Voigt1*Kyye-Voigt3*Kxxe-Voigt3*Kzze;
  Kpp(ne2,ne3)=-Voigt2*Kyze-Voigt3*Kzye;
  Kpp(ne3,ne1)=-Voigt2*Kzxe-Voigt3*Kxze;
  Kpp(ne3,ne2)=-Voigt2*Kzye-Voigt3*Kyze;
  Kpp(ne3,ne3)=-Voigt1*Kzze-Voigt3*Kxxe-Voigt3*Kyye;
end

Ktp(1,ne1)=Nw1eT*(Ne);
Ktp(2,ne2)=Nw1eT*(Ne);
if nsd==3
  Ktp(3,ne3)=Nw1eT*(Ne);
end

if nsd==2
  Krp(1,ne1)=-Nw1eT*(Ney);
  Krp(1,ne2)=+Nw1eT*(Nex);
elseif nsd==3
  Krp(1,ne2)=-Nw1eT*(Nez);
  Krp(1,ne3)=+Nw1eT*(Ney);
  Krp(2,ne1)=+Nw1eT*(Nez);
  Krp(2,ne3)=-Nw1eT*(Nex);
  Krp(3,ne1)=-Nw1eT*(Ney);
  Krp(3,ne2)=+Nw1eT*(Nex);
end

% Compute rhs
if nsd==2
  fp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg);
  fp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg);
elseif nsd==3
  fp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg)+NwezT*(Lxzeg);
  fp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg)+NwezT*(Lyzeg);
  fp(ne3,1)=NwexT*(Lxzeg)+NweyT*(Lyzeg)+NwezT*(Lzzeg);
end

ft(1,1)=Nw1eT*(wxeg./reg);
ft(2,1)=Nw1eT*(wyeg./reg);
if nsd==3
  ft(3,1)=Nw1eT*(wzeg./reg);
end

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions('Face',RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  Xfg=Nlf*Xf;
  
  % Check boundary
  isDirichlet_r=Faces.Dirichlet_r(iFace);
  isDirichlet_w_x=Faces.Dirichlet_w_x(iFace);
  isDirichlet_w_y=Faces.Dirichlet_w_y(iFace);
  if nsd==3; isDirichlet_w_z=Faces.Dirichlet_w_z(iFace); end
  isDirichlet_w=isDirichlet_w_x || isDirichlet_w_y || (nsd==3 && isDirichlet_w_z);
  
  % Indices
  nlefU1=(iFace-1)*(1+nsd)*NumFaceNodes+(1:NumFaceNodes);
  nlefU2=nlefU1+NumFaceNodes;
  nlefU3=nlefU2+NumFaceNodes;
  nlefU4=nlefU3+NumFaceNodes;
  
  % Flip face
  Node2Match1stNode1=Faces.Interior(2,iFace);
  FlipFace=max(Node2Match1stNode1);
  if FlipFace
    order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
    nlefU1=nlefU1(order);
    nlefU2=nlefU2(order);
    nlefU3=nlefU3(order);
    nlefU4=nlefU4(order);
  end
  
  % Compute variables at nodes
  Rf=Ue(nlefU1);
  Wxf=Ue(nlefU2);
  Wyf=Ue(nlefU3);
  if nsd==3
    Wzf=Ue(nlefU4);
  end
  
  % Compute variables at Gauss points
  Rfg=Nlf*Rf;
  Wxfg=Nlf*Wxf;
  Wyfg=Nlf*Wyf;
  if nsd==3
    Wzfg=Nlf*Wzf;
  end
  if isDirichlet_r
    rDfg=rD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
  end
  if isDirichlet_w
    wDfg=wD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    wDxfg=wDfg(:,1);
    wDyfg=wDfg(:,2);
    if nsd==3
      wDzfg=wDfg(:,3);
    end
  end
  
  % Compute common terms
  if isDirichlet_r
    Rfg=rDfg;
  end
  if isDirichlet_w_x
    Wxfg=wDxfg;
  end
  if isDirichlet_w_y
    Wyfg=wDyfg;
  end
  if nsd==3 && isDirichlet_w_z
    Wzfg=wDzfg;
  end
  
  % Compute basic matrices
  Nw1fT=(wfg.*N1f)';
  
  % Compute rhs
  if nsd==2
    fr(1)=fr(1)+Nw1fT*(-Wxfg./Rfg.*ny+Wyfg./Rfg.*nx);
  elseif nsd==3
    fr(1)=fr(1)+Nw1fT*(-Wyfg./Rfg.*nz+Wzfg./Rfg.*ny);
    fr(2)=fr(2)+Nw1fT*(+Wxfg./Rfg.*nz-Wzfg./Rfg.*nx);
    fr(3)=fr(3)+Nw1fT*(-Wxfg./Rfg.*ny+Wyfg./Rfg.*nx);
  end
end

% Indices
ip=1:nsd*NumElementNodesPost;
it=ip(end)+(1:nsd);
ir=it(end)+(1:qsd);

% Initialization of lhs and rhs
LhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd,nsd*NumElementNodesPost);
RhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd,1);

% Lhs for post-processing
LhsPost(ip,ip)=Kpp;
LhsPost(it,ip)=Ktp;
LhsPost(ir,ip)=Krp;

% Rhs for post-processing
RhsPost(ip,1)=fp;
RhsPost(it,1)=ft;
RhsPost(ir,1)=fr;

end