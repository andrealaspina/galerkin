classdef MagnetohydrodynamicsCURLCURL_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim+1+...
                                 NumSpaceDim-1+1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2+NumSpaceDim+1+...
                                NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim+1;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) NumSpaceDim+...
                               NumSpaceDim;

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
      
      % Consider reference pressure instead of zero
      for iFace=1:Sizes(iD).NumFaces
        Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp.*Sizes(iD).NumFaceNodes+...
          Sizes(iD).NumSpaceDim*Sizes(iD).NumFaceNodes+(1:Sizes(iD).NumFaceNodes))=...
          Parameters(iD).ReferencePressure;
      end
      Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim)=...
        Parameters(iD).ReferencePressure;
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
          Parameters(iD).Velocity(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).Pressure(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2),...
          Parameters(iD).MagneticInduction(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          zeros(Sizes(iD).NumLocalNodes,1)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
              Parameters(iD).Velocity(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).Pressure(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2),...
              Parameters(iD).MagneticInduction(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              zeros(Sizes(iD).NumLocalNodes,1)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD,Block,Elements,~,Parameters,~,~,Time,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      FacesElem=Elements(iD).Faces;
      SolutionGlobalElem=Elements(iD).SolutionGlobal;
      SolutionLocalElem=Elements(iD).SolutionLocal;
      SolutionOldElem=Elements(iD).SolutionOld;
      LhsCoef=zeros(Sizes(iD).NumElementLhsCoef,Sizes(iD).NumElements);
      RhsCoef=zeros(Sizes(iD).NumElementRhsCoef,Sizes(iD).NumElements);
      MatLocal=cell(Sizes(iD).NumElements,1);
      VecLocal=cell(Sizes(iD).NumElements,1);
      parfor iElem=1:Sizes(iD).NumElements
        [LhsGlobalElem,RhsGlobalElem,MatLocalElem,VecLocalElem]=...
          buildBlockElement(NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          Parameters,Time,RefElement.Value,Sizes); %#ok
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
        MatLocal{iElem}=MatLocalElem;
        VecLocal{iElem}=VecLocalElem;
      end
      Block(iD,iD).LhsGlobal=fsparse(Block(iD,iD).LhsRowIndices,...
                                     Block(iD,iD).LhsColIndices,LhsCoef(:));
      Block(iD,iD).RhsGlobal=fsparse(Block(iD,iD).RhsRowIndices,1,RhsCoef(:));
      Elements(iD).MatLocal=MatLocal;
      Elements(iD).VecLocal=VecLocal;
    end
    
    %% Do post-process
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,~,RefElement,Sizes)
      RefElement.PostLowHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree  ,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+1,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHighHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+2,Parameters.Degree+2,Parameters.NodesDistribution);
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
          Parameters,RefElement,Sizes);
        LhsPost{iElem}=LhsPostElem;
        RhsPost{iElem}=RhsPostElem;
      end
      Elements.LhsPost=LhsPost;
      Elements.RhsPost=RhsPost;
    end
    
    %% Store results
    function [Results]=storeResults(~,iD,iST,Results,Block,~,Parameters,~,Time,~,Sizes)
      if iST==1
        Results(iD).Time=[];
        Results(iD).ScaledStrainRate=[];
        Results(iD).Velocity=[];
        Results(iD).Pressure=[];
        Results(iD).ScaledMagneticCurl=[];
        Results(iD).MagneticInduction=[];
        Results(iD).LagrangeMultiplier=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledStrainRate(:,:,iST)=Block(iD,iD).SolutionLocal(:,1:Sizes(iD).NumVoigtComp);
      Results(iD).Velocity(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        (1:Sizes(iD).NumSpaceDim));
      Results(iD).Pressure(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        Sizes(iD).NumSpaceDim+1);
      Results(iD).ScaledMagneticCurl(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1+...
        (1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2));
      Results(iD).MagneticInduction(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        Sizes(iD).NumSpaceDim+1+Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+...
        (1:Sizes(iD).NumSpaceDim));
      Results(iD).LagrangeMultiplier(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1+...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim+1);
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).VelocityPost=Block(iD,iD).SolutionPost(:,1:Sizes(iD).NumSpaceDim);
        Results(iD).MagneticInductionPost=Block(iD,iD).SolutionPost(:,Sizes(iD).NumSpaceDim+...
          (1:Sizes(iD).NumSpaceDim));
      end
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isConvectiveFlow=strcmp(Parameters.ConvectiveFlow,'yes');
isCouplingTerms=strcmp(Parameters.CouplingTerms,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.DynamicViscosity;
lambda=-2/3*Parameters.DynamicViscosity;
eta=Parameters.MagneticDiffusivity;
mu_m=Parameters.MagneticPermeability;
r=Parameters.Density;
drdp=Parameters.DDensityDPressure;
vD=Parameters.Velocity;
pD=Parameters.Pressure;
tN=Parameters.Traction;
bD=Parameters.MagneticInduction;
qD=Parameters.LagrangeMultiplier;
eNt=Parameters.TangentialElectricField;
bNn=Parameters.NormalMagneticInduction;
Rc=Parameters.ResidualContinuity;
f=Parameters.Force;
g=Parameters.Source;
Xe=Nodes';
tauV=Parameters.StabVelocity;
tauP=Parameters.StabPressure;
tauB=Parameters.StabMagneticInduction;
tauQ=Parameters.StabLagrangeMultiplier;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
ve=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
pe=reshape(SolutionLocal(:,msd+nsd+1),[],1);
Je=reshape(SolutionLocal(:,msd+nsd+1+(1:qsd)),[],1);
be=reshape(SolutionLocal(:,msd+nsd+1+qsd+(1:nsd)),[],1);
qe=reshape(SolutionLocal(:,msd+nsd+1+qsd+nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  volde=reshape(SolutionOld(:,msd+(1:nsd),:),[],BDFo);
  polde=reshape(SolutionOld(:,msd+nsd+1,:),[],BDFo);
  bolde=reshape(SolutionOld(:,msd+nsd+1+qsd+(1:nsd),:),[],BDFo);
end

% Initialize lhs
KLL=zeros(msd*NumElementNodes,msd*NumElementNodes);
KLv=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KLV=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kvv=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kvp=zeros(nsd*NumElementNodes,NumElementNodes);
KvJ=zeros(nsd*NumElementNodes,qsd*NumElementNodes);
Kvb=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KvV=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KvP=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Kpv=zeros(NumElementNodes,nsd*NumElementNodes);
Kpp=zeros(NumElementNodes,NumElementNodes);
KpV=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KpP=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KJJ=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KJb=zeros(qsd*NumElementNodes,nsd*NumElementNodes);
KJB=zeros(qsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
Kbv=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KbJ=zeros(nsd*NumElementNodes,qsd*NumElementNodes);
Kbb=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kbq=zeros(nsd*NumElementNodes,NumElementNodes);
KbV=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KbB=zeros(nsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KbQ=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Kqb=zeros(NumElementNodes,nsd*NumElementNodes);
Kqq=zeros(NumElementNodes,NumElementNodes);
KqQ=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KVL=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KVv=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KVp=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementNodes);
KVV=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KPp=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KPP=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KBJ=zeros((nsd-1)*NumElementFaces*NumFaceNodes,qsd*NumElementNodes);
KBb=zeros((nsd-1)*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KBV=zeros((nsd-1)*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KBB=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KBQ=zeros((nsd-1)*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KQb=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KQq=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KQQ=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fL=zeros(msd*NumElementNodes,1);
fv=zeros(nsd*NumElementNodes,1);
fp=zeros(NumElementNodes,1);
fJ=zeros(qsd*NumElementNodes,1);
fb=zeros(nsd*NumElementNodes,1);
fq=zeros(NumElementNodes,1);
fV=zeros(nsd*NumElementFaces*NumFaceNodes,1);
fP=zeros(NumElementFaces*NumFaceNodes,1);
fB=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);
fQ=zeros(NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions('Element',RefElement,RefElement,Xe,nsd);

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
vxe=ve(ne1);
vye=ve(ne2);
if nsd==3
  vze=ve(ne3);
end
if nsd==2
  Jze=Je(ne1);
elseif nsd==3
  Jxe=Je(ne1);
  Jye=Je(ne2);
  Jze=Je(ne3);
end
bxe=be(ne1);
bye=be(ne2);
if nsd==3
  bze=be(ne3);
end
if isTimeDependent
  voldxe=volde(ne1,:);
  voldye=volde(ne2,:);
  if nsd==3
    voldze=volde(ne3,:);
  end
  boldxe=bolde(ne1,:);
  boldye=bolde(ne2,:);
  if nsd==3
    boldze=bolde(ne3,:);
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
vxeg=Ne*vxe;
vyeg=Ne*vye;
if nsd==3
  vzeg=Ne*vze;
end
peg=Ne*pe;
if nsd==2
  Jzeg=Ne*Jze;
elseif nsd==3
  Jxeg=Ne*Jxe;
  Jyeg=Ne*Jye;
  Jzeg=Ne*Jze;
end
bxeg=Ne*bxe;
byeg=Ne*bye;
if nsd==3
  bzeg=Ne*bze;
end
qeg=Ne*qe;
if isTimeDependent
  voldxeg=Ne*voldxe;
  voldyeg=Ne*voldye;
  if nsd==3
    voldzeg=Ne*voldze;
  end
  poldeg=Ne*polde;
  boldxeg=Ne*boldxe;
  boldyeg=Ne*boldye;
  if nsd==3
    boldzeg=Ne*boldze;
  end
end
Rceg=Rc(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
feg=f(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
fxeg=feg(:,1);
fyeg=feg(:,2);
if nsd==3
  fzeg=feg(:,3);
end
geg=g(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
gxeg=geg(:,1);
gyeg=geg(:,2);
if nsd==3
  gzeg=geg(:,3);
end
reg=r(peg);
drdpeg=drdp(peg);

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
  VoigtV1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  VoigtV2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  VoigtV3=sqrt(mu);
elseif nsd==3
  VoigtV1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  VoigtV2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  VoigtV3=sqrt(mu);
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
  KLv(ne1,ne1)=VoigtV1*Cxe;
  KLv(ne2,ne1)=VoigtV2*Cxe;
  KLv(ne3,ne1)=VoigtV3*Cye;
  KLv(ne1,ne2)=VoigtV2*Cye;
  KLv(ne2,ne2)=VoigtV1*Cye;
  KLv(ne3,ne2)=VoigtV3*Cxe;
elseif nsd==3
  KLv(ne1,ne1)=VoigtV1*Cxe;
  KLv(ne2,ne1)=VoigtV2*Cxe;
  KLv(ne3,ne1)=VoigtV2*Cxe;
  KLv(ne4,ne1)=VoigtV3*Cye;
  KLv(ne5,ne1)=VoigtV3*Cze;
  KLv(ne1,ne2)=VoigtV2*Cye;
  KLv(ne2,ne2)=VoigtV1*Cye;
  KLv(ne3,ne2)=VoigtV2*Cye;
  KLv(ne4,ne2)=VoigtV3*Cxe;
  KLv(ne6,ne2)=VoigtV3*Cze;
  KLv(ne1,ne3)=VoigtV2*Cze;
  KLv(ne2,ne3)=VoigtV2*Cze;
  KLv(ne3,ne3)=VoigtV1*Cze;
  KLv(ne5,ne3)=VoigtV3*Cxe;
  KLv(ne6,ne3)=VoigtV3*Cye;
end

if isTimeDependent
  Kvv(ne1,ne1)=alpha(1)/dt*NweT*((reg).*Ne);
  Kvv(ne2,ne2)=alpha(1)/dt*NweT*((reg).*Ne);
  if nsd==3
    Kvv(ne3,ne3)=alpha(1)/dt*NweT*((reg).*Ne);
  end
  
  Kvv(ne1,ne1)=Kvv(ne1,ne1)+NweT*((drdpeg.*(1/dt*peg*alpha(1)+1/dt*poldeg*alpha(2:BDFo+1,1))).*Ne);
  Kvv(ne2,ne2)=Kvv(ne2,ne2)+NweT*((drdpeg.*(1/dt*peg*alpha(1)+1/dt*poldeg*alpha(2:BDFo+1,1))).*Ne);
  if nsd==3
    Kvv(ne3,ne3)=Kvv(ne3,ne3)+NweT*((drdpeg.*(1/dt*peg*alpha(1)+1/dt*poldeg*alpha(2:BDFo+1,1))).*Ne);
  end
  
  Kvp(ne1,ne1)=NweT*((drdpeg.*(1/dt*vxeg*alpha(1)+1/dt*voldxeg*alpha(2:BDFo+1,1))).*Ne);
  Kvp(ne2,ne1)=NweT*((drdpeg.*(1/dt*vyeg*alpha(1)+1/dt*voldyeg*alpha(2:BDFo+1,1))).*Ne);
  if nsd==3
    Kvp(ne3,ne1)=NweT*((drdpeg*(1/dt*vzeg*alpha(1)+1/dt*voldzeg*alpha(2:BDFo+1,1))).*Ne);
  end
  
  Kvp(ne1,ne1)=Kvp(ne1,ne1)+alpha(1)/dt*NweT*((drdpeg.*vxeg).*Ne);
  Kvp(ne2,ne1)=Kvp(ne2,ne1)+alpha(1)/dt*NweT*((drdpeg.*vyeg).*Ne);
  if nsd==3
    Kvp(ne3,ne1)=Kvp(ne3,ne1)+alpha(1)/dt*NweT*((drdpeg.*vzeg).*Ne);
  end
end

if isConvectiveFlow
  if nsd==2
    Kvv(ne1,ne1)=Kvv(ne1,ne1)-NwexT*((2*reg.*vxeg).*Ne)...
                             -NweyT*((reg.*vyeg).*Ne);
    Kvv(ne2,ne1)=Kvv(ne2,ne1)-NwexT*((reg.*vyeg).*Ne);
    Kvv(ne1,ne2)=Kvv(ne1,ne2)-NweyT*((reg.*vxeg).*Ne);
    Kvv(ne2,ne2)=Kvv(ne2,ne2)-NwexT*((reg.*vxeg).*Ne)...
                             -NweyT*((2*reg.*vyeg).*Ne);
  elseif nsd==3
    Kvv(ne1,ne1)=Kvv(ne1,ne1)-NwexT*((2*reg.*vxeg).*Ne)...
                             -NweyT*((reg.*vyeg).*Ne)...
                             -NwezT*((reg.*vzeg).*Ne);
    Kvv(ne2,ne1)=Kvv(ne2,ne1)-NwexT*((reg.*vyeg).*Ne);
    Kvv(ne3,ne1)=Kvv(ne3,ne1)-NwezT*((reg.*vzeg).*Ne);
    Kvv(ne1,ne2)=Kvv(ne1,ne2)-NweyT*((reg.*vxeg).*Ne);
    Kvv(ne2,ne2)=Kvv(ne2,ne2)-NwexT*((reg.*vxeg).*Ne)...
                             -NweyT*((2*reg.*vyeg).*Ne)...
                             -NwezT*((reg.*vzeg).*Ne);
    Kvv(ne3,ne2)=Kvv(ne3,ne2)-NwezT*((reg.*vzeg).*Ne);
    Kvv(ne1,ne3)=Kvv(ne1,ne3)-NwezT*((reg.*vxeg).*Ne);
    Kvv(ne2,ne3)=Kvv(ne2,ne3)-NwezT*((reg.*vyeg).*Ne);
    Kvv(ne3,ne3)=Kvv(ne3,ne3)-NwexT*((reg.*vxeg).*Ne)...
                             -NweyT*((reg.*vyeg).*Ne)...
                             -NwezT*((2*reg.*vzeg).*Ne);
  end
  
  if nsd==2
    Kvp(ne1,ne1)=Kvp(ne1,ne1)-NwexT*((drdpeg.*vxeg.*vxeg).*Ne)...
                             -NweyT*((drdpeg.*vxeg.*vyeg).*Ne);
    Kvp(ne2,ne1)=Kvp(ne2,ne1)-NwexT*((drdpeg.*vyeg.*vxeg).*Ne)...
                             -NweyT*((drdpeg.*vyeg.*vyeg).*Ne);
  elseif nsd==3
    Kvp(ne1,ne1)=Kvp(ne1,ne1)-NwexT*((drdpeg.*vxeg.*vxeg).*Ne)...
                             -NweyT*((drdpeg.*vxeg.*vyeg).*Ne)...
                             -NwezT*((drdpeg.*vxeg.*vzeg).*Ne);
    Kvp(ne2,ne1)=Kvp(ne2,ne1)-NwexT*((drdpeg.*vyeg.*vxeg).*Ne)...
                             -NweyT*((drdpeg.*vyeg.*vyeg).*Ne)...
                             -NwezT*((drdpeg.*vyeg.*vzeg).*Ne);
    Kvp(ne3,ne1)=Kvp(ne3,ne1)-NwexT*((drdpeg.*vzeg.*vxeg).*Ne)...
                             -NweyT*((drdpeg.*vzeg.*vyeg).*Ne)...
                             -NwezT*((drdpeg.*vzeg.*vzeg).*Ne);
  end
end

Kvp(ne1,ne1)=Kvp(ne1,ne1)+CxeT;
Kvp(ne2,ne1)=Kvp(ne2,ne1)+CyeT;
if nsd==3
  Kvp(ne3,ne1)=Kvp(ne3,ne1)+CzeT;
end

if isCouplingTerms
  if nsd==2
    KvJ(ne1,ne1)=-1/(mu_m*sqrt(eta))*NweT*((-byeg).*Ne);
    KvJ(ne2,ne1)=-1/(mu_m*sqrt(eta))*NweT*((+bxeg).*Ne);
  elseif nsd==3
    KvJ(ne2,ne1)=-1/(mu_m*sqrt(eta))*NweT*((-bzeg).*Ne);
    KvJ(ne3,ne1)=-1/(mu_m*sqrt(eta))*NweT*((+byeg).*Ne);
    KvJ(ne1,ne2)=-1/(mu_m*sqrt(eta))*NweT*((+bzeg).*Ne);
    KvJ(ne3,ne2)=-1/(mu_m*sqrt(eta))*NweT*((-bxeg).*Ne);
    KvJ(ne1,ne3)=-1/(mu_m*sqrt(eta))*NweT*((-byeg).*Ne);
    KvJ(ne2,ne3)=-1/(mu_m*sqrt(eta))*NweT*((+bxeg).*Ne);
  end

  if nsd==2
    Kvb(ne2,ne1)=-1/(mu_m*sqrt(eta))*NweT*((+Jzeg).*Ne);
    Kvb(ne1,ne2)=-1/(mu_m*sqrt(eta))*NweT*((-Jzeg).*Ne);
  elseif nsd==3
    Kvb(ne2,ne1)=-1/(mu_m*sqrt(eta))*NweT*((+Jzeg).*Ne);
    Kvb(ne3,ne1)=-1/(mu_m*sqrt(eta))*NweT*((-Jyeg).*Ne);
    Kvb(ne1,ne2)=-1/(mu_m*sqrt(eta))*NweT*((-Jzeg).*Ne);
    Kvb(ne3,ne2)=-1/(mu_m*sqrt(eta))*NweT*((+Jxeg).*Ne);
    Kvb(ne1,ne3)=-1/(mu_m*sqrt(eta))*NweT*((+Jyeg).*Ne);
    Kvb(ne2,ne3)=-1/(mu_m*sqrt(eta))*NweT*((-Jxeg).*Ne);
  end
end

if isTimeDependent
  Kpp(ne1,ne1)=alpha(1)/dt*NweT*((drdpeg).*Ne);
end

Kpv(ne1,ne1)=-NwexT*((reg).*Ne);
Kpv(ne1,ne2)=-NweyT*((reg).*Ne);
if nsd==3
  Kpv(ne1,ne3)=-NwezT*((reg).*Ne);
end

Kpp(ne1,ne1)=Kpp(ne1,ne1)-NwexT*((drdpeg.*vxeg).*Ne)...
                         -NweyT*((drdpeg.*vyeg).*Ne);
if nsd==3
  Kpp(ne1,ne1)=Kpp(ne1,ne1)-NwezT*((drdpeg.*vzeg).*Ne);
end

KJJ(ne1,ne1)=-Me;
if nsd==3
  KJJ(ne2,ne2)=-Me;
  KJJ(ne3,ne3)=-Me;
end

if nsd==2
  KJb(ne1,ne1)=+sqrt(eta)*Cye;
  KJb(ne1,ne2)=-sqrt(eta)*Cxe;
elseif nsd==3
  KJb(ne2,ne1)=-sqrt(eta)*Cze;
  KJb(ne3,ne1)=+sqrt(eta)*Cye;
  KJb(ne1,ne2)=+sqrt(eta)*Cze;
  KJb(ne3,ne2)=-sqrt(eta)*Cxe;
  KJb(ne1,ne3)=-sqrt(eta)*Cye;
  KJb(ne2,ne3)=+sqrt(eta)*Cxe;
end

if isTimeDependent
  Kbb(ne1,ne1)=alpha(1)/dt*Me;
  Kbb(ne2,ne2)=alpha(1)/dt*Me;
  if nsd==3
    Kbb(ne3,ne3)=alpha(1)/dt*Me;
  end
end

if isCouplingTerms
  if nsd==2
    Kbb(ne1,ne1)=Kbb(ne1,ne1)-NweyT*((+vyeg).*Ne);
    Kbb(ne2,ne1)=Kbb(ne2,ne1)-NwexT*((-vyeg).*Ne);
    Kbb(ne1,ne2)=Kbb(ne1,ne2)-NweyT*((-vxeg).*Ne);
    Kbb(ne2,ne2)=Kbb(ne2,ne2)-NwexT*((+vxeg).*Ne);
  elseif nsd==3
    Kbb(ne1,ne1)=Kbb(ne1,ne1)-NweyT*((+vyeg).*Ne)...
                             -NwezT*((+vzeg).*Ne);
    Kbb(ne2,ne1)=Kbb(ne2,ne1)-NwexT*((-vyeg).*Ne);
    Kbb(ne3,ne1)=Kbb(ne3,ne1)-NwexT*((-vzeg).*Ne);
    Kbb(ne1,ne2)=Kbb(ne1,ne2)-NweyT*((-vxeg).*Ne);
    Kbb(ne2,ne2)=Kbb(ne2,ne2)-NwexT*((+vxeg).*Ne)...
                             -NwezT*((+vzeg).*Ne);
    Kbb(ne3,ne2)=Kbb(ne3,ne2)-NweyT*((-vzeg).*Ne);
    Kbb(ne1,ne3)=Kbb(ne1,ne3)-NwezT*((-vxeg).*Ne);
    Kbb(ne2,ne3)=Kbb(ne2,ne3)-NwezT*((-vyeg).*Ne);
    Kbb(ne3,ne3)=Kbb(ne3,ne3)-NwexT*((+vxeg).*Ne)...
                             -NweyT*((+vyeg).*Ne);
  end
  
  if nsd==2
    Kbv(ne1,ne1)=-NweyT*((-byeg).*Ne);
    Kbv(ne2,ne1)=-NwexT*((+byeg).*Ne);
    Kbv(ne1,ne2)=-NweyT*((+bxeg).*Ne);
    Kbv(ne2,ne2)=-NwexT*((-bxeg).*Ne);
  elseif nsd==3
    Kbv(ne1,ne1)=-NweyT*((-byeg).*Ne)...
                 -NwezT*((-bzeg).*Ne);
    Kbv(ne2,ne1)=-NwexT*((+byeg).*Ne);
    Kbv(ne3,ne1)=-NwexT*((+bzeg).*Ne);
    Kbv(ne1,ne2)=-NweyT*((+bxeg).*Ne);
    Kbv(ne2,ne2)=-NwexT*((-bxeg).*Ne)...
                 -NwezT*((-bzeg).*Ne);
    Kbv(ne3,ne2)=-NweyT*((+bzeg).*Ne);
    Kbv(ne1,ne3)=-NwezT*((+bxeg).*Ne);
    Kbv(ne2,ne3)=-NwezT*((+byeg).*Ne);
    Kbv(ne3,ne3)=-NwexT*((-bxeg).*Ne)...
                 -NweyT*((-byeg).*Ne);
  end
end

if nsd==2
  KbJ(ne1,ne1)=-sqrt(eta)*Cye;
  KbJ(ne2,ne1)=+sqrt(eta)*Cxe;
elseif nsd==3
  KbJ(ne2,ne1)=-sqrt(eta)*Cze;
  KbJ(ne3,ne1)=+sqrt(eta)*Cye;
  KbJ(ne1,ne2)=+sqrt(eta)*Cze;
  KbJ(ne3,ne2)=-sqrt(eta)*Cxe;
  KbJ(ne1,ne3)=-sqrt(eta)*Cye;
  KbJ(ne2,ne3)=+sqrt(eta)*Cxe;
end

Kbq(ne1,ne1)=-Cxe;
Kbq(ne2,ne1)=-Cye;
if nsd==3
  Kbq(ne3,ne1)=-Cze;
end

Kqb(ne1,ne1)=Cxe;
Kqb(ne1,ne2)=Cye;
if nsd==3
  Kqb(ne1,ne3)=Cze;
end

% Compute rhs
if nsd==2
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(VoigtV1*vxeg)...
            -NweyT*(VoigtV2*vyeg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(VoigtV2*vxeg)...
            -NweyT*(VoigtV1*vyeg);
  fL(ne3,1)=+NweT*(Lxyeg)...
            -NweyT*(VoigtV3*vxeg)...
            -NwexT*(VoigtV3*vyeg);
elseif nsd==3
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(VoigtV1*vxeg)...
            -NweyT*(VoigtV2*vyeg)...
            -NwezT*(VoigtV2*vzeg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(VoigtV2*vxeg)...
            -NweyT*(VoigtV1*vyeg)...
            -NwezT*(VoigtV2*vzeg);
  fL(ne3,1)=+NweT*(Lzzeg)...
            -NwexT*(VoigtV2*vxeg)...
            -NweyT*(VoigtV2*vyeg)...
            -NwezT*(VoigtV1*vzeg);
  fL(ne4,1)=+NweT*(Lxyeg)...
            -NwexT*(VoigtV3*vyeg)...
            -NweyT*(VoigtV3*vxeg);
  fL(ne5,1)=+NweT*(Lxzeg)...
            -NwexT*(VoigtV3*vzeg)...
            -NwezT*(VoigtV3*vxeg);
  fL(ne6,1)=+NweT*(Lyzeg)...
            -NweyT*(VoigtV3*vzeg)...
            -NwezT*(VoigtV3*vyeg);
end

if isTimeDependent
  fv(ne1,1)=-NweT*(r(peg).*(1/dt*vxeg*alpha(1)...
                           +1/dt*voldxeg*alpha(2:BDFo+1,1))...
               +drdp(peg).*(1/dt*peg*alpha(1)...
                           +1/dt*poldeg*alpha(2:BDFo+1,1)).*vxeg);
  fv(ne2,1)=-NweT*(r(peg).*(1/dt*vyeg*alpha(1)...
                           +1/dt*voldyeg*alpha(2:BDFo+1,1))...
               +drdp(peg).*(1/dt*peg*alpha(1)...
                           +1/dt*poldeg*alpha(2:BDFo+1,1)).*vyeg);
  if nsd==3
    fv(ne3,1)=-NweT*(r(peg).*(1/dt*vzeg*alpha(1)...
                             +1/dt*voldzeg*alpha(2:BDFo+1,1))...
                 +drdp(peg).*(1/dt*peg*alpha(1)...
                             +1/dt*poldeg*alpha(2:BDFo+1,1)).*vzeg);
  end
end

if isConvectiveFlow
  if nsd==2
    fv(ne1,1)=fv(ne1,1)+NwexT*(r(peg).*vxeg.*vxeg)...
                       +NweyT*(r(peg).*vxeg.*vyeg);
    fv(ne2,1)=fv(ne2,1)+NwexT*(r(peg).*vyeg.*vxeg)...
                       +NweyT*(r(peg).*vyeg.*vyeg);
  elseif nsd==3
    fv(ne1,1)=fv(ne1,1)+NwexT*(r(peg).*vxeg.*vxeg)...
                       +NweyT*(r(peg).*vxeg.*vyeg)...
                       +NwezT*(r(peg).*vxeg.*vzeg);
    fv(ne2,1)=fv(ne2,1)+NwexT*(r(peg).*vyeg.*vxeg)...
                       +NweyT*(r(peg).*vyeg.*vyeg)...
                       +NwezT*(r(peg).*vyeg.*vzeg);
    fv(ne3,1)=fv(ne3,1)+NwexT*(r(peg).*vzeg.*vxeg)...
                       +NweyT*(r(peg).*vzeg.*vyeg)...
                       +NwezT*(r(peg).*vzeg.*vzeg);
  end
end

if nsd==2
  fv(ne1,1)=fv(ne1,1)-NweT*(+VoigtV1*(Nex*Lxxe)+VoigtV2*(Nex*Lyye)...
                            +VoigtV3*(Ney*Lxye)...
                            +(Nex*pe));
  fv(ne2,1)=fv(ne2,1)-NweT*(+VoigtV2*(Ney*Lxxe)+VoigtV1*(Ney*Lyye)...
                            +VoigtV3*(Nex*Lxye)...
                            +(Ney*pe));
elseif nsd==3
  fv(ne1,1)=fv(ne1,1)-NweT*(+VoigtV1*(Nex*Lxxe)+VoigtV2*(Nex*Lyye)+VoigtV2*(Nex*Lzze)...
                            +VoigtV3*(Ney*Lxye)+VoigtV3*(Nez*Lxze)...
                            +(Nex*pe));
  fv(ne2,1)=fv(ne2,1)-NweT*(+VoigtV2*(Ney*Lxxe)+VoigtV1*(Ney*Lyye)+VoigtV2*(Ney*Lzze)...
                            +VoigtV3*(Nex*Lxye)+VoigtV3*(Nez*Lyze)...
                            +(Ney*pe));
  fv(ne3,1)=fv(ne3,1)-NweT*(+VoigtV2*(Nez*Lxxe)+VoigtV2*(Nez*Lyye)+VoigtV1*(Nez*Lzze)...
                            +VoigtV3*(Nex*Lxze)+VoigtV3*(Ney*Lyze)...
                            +(Nez*pe));
end

if isCouplingTerms
  if nsd==2
    fv(ne1,1)=fv(ne1,1)+NweT*(1/(mu_m*sqrt(eta))*(-byeg.*Jzeg));
    fv(ne2,1)=fv(ne2,1)+NweT*(1/(mu_m*sqrt(eta))*(+bxeg.*Jzeg));
  elseif nsd==3
    fv(ne1,1)=fv(ne1,1)+NweT*(1/(mu_m*sqrt(eta))*(+bzeg.*Jyeg...
                                                  -byeg.*Jzeg));
    fv(ne2,1)=fv(ne2,1)+NweT*(1/(mu_m*sqrt(eta))*(-bzeg.*Jxeg...
                                                  +bxeg.*Jzeg));
    fv(ne3,1)=fv(ne3,1)+NweT*(1/(mu_m*sqrt(eta))*(+byeg.*Jxeg...
                                                  -bxeg.*Jyeg));
  end
end

fv(ne1,1)=fv(ne1,1)+NweT*(fxeg);
fv(ne2,1)=fv(ne2,1)+NweT*(fyeg);
if nsd==3
  fv(ne3,1)=fv(ne3,1)+NweT*(fzeg);
end

if isTimeDependent
  fp(ne1,1)=-NweT*(drdp(peg).*(1/dt*peg*alpha(1)...
                              +1/dt*poldeg*alpha(2:BDFo+1,1)));
end

fp(ne1,1)=fp(ne1,1)+NwexT*(r(peg).*vxeg)...
                   +NweyT*(r(peg).*vyeg);
if nsd==3
  fp(ne1,1)=fp(ne1,1)+NwezT*(r(peg).*vzeg);
end

fp(ne1,1)=fp(ne1,1)+NweT*(Rceg);

if nsd==2
  fJ(ne1,1)=+NweT*(Jzeg);
elseif nsd==3
  fJ(ne1,1)=+NweT*(Jxeg);
  fJ(ne2,1)=+NweT*(Jyeg);
  fJ(ne3,1)=+NweT*(Jzeg);
end

if nsd==2
  fJ(ne1,1)=fJ(ne1,1)+NwexT*(sqrt(eta)*byeg)...
                     -NweyT*(sqrt(eta)*bxeg);
elseif nsd==3
  fJ(ne1,1)=fJ(ne1,1)+NweyT*(sqrt(eta)*bzeg)...
                     -NwezT*(sqrt(eta)*byeg);
  fJ(ne2,1)=fJ(ne2,1)+NwezT*(sqrt(eta)*bxeg)...
                     -NwexT*(sqrt(eta)*bzeg);
  fJ(ne3,1)=fJ(ne3,1)+NwexT*(sqrt(eta)*byeg)...
                     -NweyT*(sqrt(eta)*bxeg);
end

if isTimeDependent
  fb(ne1,1)=-NweT*(1/dt*bxeg*alpha(1)...
                  +1/dt*boldxeg*alpha(2:BDFo+1,1));
  fb(ne2,1)=-NweT*(1/dt*byeg*alpha(1)...
                  +1/dt*boldyeg*alpha(2:BDFo+1,1));
  if nsd==3
    fb(ne3,1)=-NweT*(1/dt*bzeg*alpha(1)...
                    +1/dt*boldzeg*alpha(2:BDFo+1,1));
  end
end

if isCouplingTerms
  if nsd==2
    fb(ne1,1)=fb(ne1,1)+NweyT*(bxeg.*vyeg-vxeg.*byeg);
    fb(ne2,1)=fb(ne2,1)+NwexT*(byeg.*vxeg-vyeg.*bxeg);
  elseif nsd==3
    fb(ne1,1)=fb(ne1,1)+NweyT*(bxeg.*vyeg-vxeg.*byeg)...
                       +NwezT*(bxeg.*vzeg-vxeg.*bzeg);
    fb(ne2,1)=fb(ne2,1)+NwexT*(byeg.*vxeg-vyeg.*bxeg)...
                       +NwezT*(byeg.*vzeg-vyeg.*bzeg);
    fb(ne3,1)=fb(ne3,1)+NwexT*(bzeg.*vxeg-vzeg.*bxeg)...
                       +NweyT*(bzeg.*vyeg-vzeg.*byeg);
  end
end

if nsd==2
  fb(ne1,1)=fb(ne1,1)+NweyT*(sqrt(eta)*Jzeg);
  fb(ne2,1)=fb(ne2,1)-NwexT*(sqrt(eta)*Jzeg);
elseif nsd==3
  fb(ne1,1)=fb(ne1,1)-NwezT*(sqrt(eta)*Jyeg)...
                     +NweyT*(sqrt(eta)*Jzeg);
  fb(ne2,1)=fb(ne2,1)-NwexT*(sqrt(eta)*Jzeg)...
                     +NwezT*(sqrt(eta)*Jxeg);
  fb(ne3,1)=fb(ne3,1)-NweyT*(sqrt(eta)*Jxeg)...
                     +NwexT*(sqrt(eta)*Jyeg);
end

fb(ne1,1)=fb(ne1,1)+NwexT*(qeg);
fb(ne2,1)=fb(ne2,1)+NweyT*(qeg);
if nsd==3
  fb(ne3,1)=fb(ne3,1)+NwezT*(qeg);
end

fb(ne1,1)=fb(ne1,1)+NweT*(gxeg);
fb(ne2,1)=fb(ne2,1)+NweT*(gyeg);
if nsd==3
  fb(ne3,1)=fb(ne3,1)+NweT*(gzeg);
end

fq(ne1,1)=-NwexT*(bxeg)...
          -NweyT*(byeg);
if nsd==3
  fq(ne1,1)=fq(ne1,1)-NwezT*(bzeg);
end

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=true;
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement.FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions('Face',RefElement,RefElement,Xf,nsd);
    
    % Check boundary
    isExterior=Faces.Exterior(iFace);
    isDirichlet_v_x=Faces.Dirichlet_v_x(iFace);
    isDirichlet_v_y=Faces.Dirichlet_v_y(iFace);
    if nsd==3; isDirichlet_v_z=Faces.Dirichlet_v_z(iFace); end
    isDirichlet_p=Faces.Dirichlet_p(iFace);
    isNeumann_t_x=Faces.Neumann_t_x(iFace);
    isNeumann_t_y=Faces.Neumann_t_y(iFace);
    if nsd==3; isNeumann_t_z=Faces.Neumann_t_z(iFace); end
    isDirichlet_m=Faces.Dirichlet_m(iFace);
    isNatural_m=Faces.Natural_m(iFace);
    isDirichlet_v=isDirichlet_v_x || isDirichlet_v_y || (nsd==3 && isDirichlet_v_z);
    isNeumann_t=isNeumann_t_x || isNeumann_t_y || (nsd==3 && isNeumann_t_z);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nf4=nf3+NumElementNodes;
    nf5=nf4+NumElementNodes;
    nf6=nf5+NumElementNodes;
    nefU1=(iFace-1)*(nsd+1+nsd-1+1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
    nefU4=nefU3+NumFaceNodes;
    nefU5=nefU4+NumFaceNodes;
    nefU6=nefU5+NumFaceNodes;
    nefU7=nefU6+NumFaceNodes;
    nefV1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefV2=nefV1+NumFaceNodes;
    nefV3=nefV2+NumFaceNodes;
    nefP1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefB1=(iFace-1)*(nsd-1)*NumFaceNodes+(1:NumFaceNodes);
    nefB2=nefB1+NumFaceNodes;
    nefQ1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefR=1:nsd;
    
    % Rotation operator in 2D (curved faces)
    if nsd==2
      Rxt1=-ny;
      Ryt1=+nx;
    end
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
      nefU1=nefU1(order);
      nefU2=nefU2(order);
      nefU3=nefU3(order);
      nefU4=nefU4(order);
      nefU5=nefU5(order);
      nefU6=nefU6(order);
      nefU7=nefU7(order);
      nefV1=nefV1(order);
      nefV2=nefV2(order);
      nefV3=nefV3(order);
      nefP1=nefP1(order);
      nefB1=nefB1(order);
      nefB2=nefB2(order);
      nefQ1=nefQ1(order);
      if nsd==2
        Rxt1=-Rxt1;
        Ryt1=-Ryt1;
      elseif nsd==3
        nefR=circshift(flip(nefR(1:nsd)),Node2Match1stNode1);
      end
    end
    
    % Gram-Schmidt orthonormalization in 3D (straight faces)
    if nsd==3
      V=Xf(nefR,1:nsd)';
      W=(V(:,2:end)-V(:,1))./vecnorm(V(:,2:end)-V(:,1));
      R(:,1)=W(:,1);
      R(:,2)=(W(:,2)-(R(:,1)'*W(:,2))*R(:,1))/norm(W(:,2)-(R(:,1)'*W(:,2))*R(:,1));
    end
    
    % Get rotation operator
    if nsd==3
      Rxt1=R(1,1); Rxt2=R(1,2);
      Ryt1=R(2,1); Ryt2=R(2,2);
      Rzt1=R(3,1); Rzt2=R(3,2);
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
    vxf=vxe(nf1);
    vyf=vye(nf1);
    if nsd==3
      vzf=vze(nf1);
    end
    pf=pe(nf1);
    if nsd==2
      Jzf=Jze(nf1);
    elseif nsd==3
      Jxf=Jxe(nf1);
      Jyf=Jye(nf1);
      Jzf=Jze(nf1);
    end
    bxf=bxe(nf1);
    byf=bye(nf1);
    if nsd==3
      bzf=bze(nf1);
    end
    qf=qe(nf1);
    Vxf=Ue(nefU1);
    Vyf=Ue(nefU2);
    if nsd==3
      Vzf=Ue(nefU3);
    end
    if nsd==2
      Pf=Ue(nefU3);
    elseif nsd==3
      Pf=Ue(nefU4);
    end
    if nsd==2
      Bt1f=Ue(nefU4);
    elseif nsd==3
      Bt1f=Ue(nefU5);
      Bt2f=Ue(nefU6);
    end
    if nsd==2
      Qf=Ue(nefU5);
    elseif nsd==3
      Qf=Ue(nefU7);
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
    vxfg=Nf*vxf;
    vyfg=Nf*vyf;
    if nsd==3
      vzfg=Nf*vzf;
    end
    pfg=Nf*pf;
    if nsd==2
      Jzfg=Nf*Jzf;
    elseif nsd==3
      Jxfg=Nf*Jxf;
      Jyfg=Nf*Jyf;
      Jzfg=Nf*Jzf;
    end
    bxfg=Nf*bxf;
    byfg=Nf*byf;
    if nsd==3
      bzfg=Nf*bzf;
    end
    qfg=Nf*qf;
    Vxfg=Nf*Vxf;
    Vyfg=Nf*Vyf;
    if nsd==3
      Vzfg=Nf*Vzf;
    end
    Pfg=Nf*Pf;
    Bt1fg=Nf*Bt1f;
    if nsd==3
      Bt2fg=Nf*Bt2f;
    end
    if nsd==2
      Btxfg=Rxt1.*Bt1fg;
      Btyfg=Ryt1.*Bt1fg;
    elseif nsd==3
      Btxfg=Rxt1.*Bt1fg+Rxt2.*Bt2fg;
      Btyfg=Ryt1.*Bt1fg+Ryt2.*Bt2fg;
      Btzfg=Rzt1.*Bt1fg+Rzt2.*Bt2fg;
    end
    Qfg=Nf*Qf;
    if isDirichlet_v
      vDfg=vD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      vDxfg=vDfg(:,1);
      vDyfg=vDfg(:,2);
      if nsd==3
        vDzfg=vDfg(:,3);
      end
    end
    if isDirichlet_p
      pDfg=pD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isNeumann_t
      tNfg=tN(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      tNxfg=tNfg(:,1);
      tNyfg=tNfg(:,2);
      if nsd==3
        tNzfg=tNfg(:,3);
      end
    end
    if isDirichlet_m
      bDfg=bD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      if nsd==2
        bDt1fg=Rxt1.*bDfg(:,1)+Ryt1.*bDfg(:,2);
      elseif nsd==3
        bDt1fg=Rxt1.*bDfg(:,1)+Ryt1.*bDfg(:,2)+Rzt1.*bDfg(:,3);
        bDt2fg=Rxt2.*bDfg(:,1)+Ryt2.*bDfg(:,2)+Rzt2.*bDfg(:,3);
      end
      if nsd==2
        bDtxfg=Rxt1.*bDt1fg;
        bDtyfg=Ryt1.*bDt1fg;
      elseif nsd==3
        bDtxfg=Rxt1.*bDt1fg+Rxt2.*bDt2fg;
        bDtyfg=Ryt1.*bDt1fg+Ryt2.*bDt2fg;
        bDtzfg=Rzt1.*bDt1fg+Rzt2.*bDt2fg;
      end
      qDfg=qD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isNatural_m
      eNtfg=eNt(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      eNtxfg=eNtfg(:,1);
      eNtyfg=eNtfg(:,2);
      if nsd==3
        eNtzfg=eNtfg(:,3);
      end
      bNnfg=bNn(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    rfg=r(Pfg);
    drdPfg=drdp(Pfg);
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    Mf=NwfT*Nf;
    Mf=(Mf+Mf')/2;
    Mfnx=NwfT*(nx.*Nf);
    Mfny=NwfT*(ny.*Nf);
    if nsd==3
      Mfnz=NwfT*(nz.*Nf);
    end
    if nsd==2
      Mfnxnx=NwfT*(nx.*nx.*Nf);
      Mfnxny=NwfT*(nx.*ny.*Nf);
      Mfnyny=NwfT*(ny.*ny.*Nf);
    elseif nsd==3
      Mfnxnx=NwfT*(nx.*nx.*Nf);
      Mfnxny=NwfT*(nx.*ny.*Nf);
      Mfnxnz=NwfT*(nx.*nz.*Nf);
      Mfnyny=NwfT*(ny.*ny.*Nf);
      Mfnynz=NwfT*(ny.*nz.*Nf);
      Mfnznz=NwfT*(nz.*nz.*Nf);
    end
    if nsd==2
      MfRxt1=NwfT*(Rxt1.*Nf);
      MfRyt1=NwfT*(Ryt1.*Nf);
    elseif nsd==3
      MfRxt1=NwfT*(Rxt1.*Nf);
      MfRyt1=NwfT*(Ryt1.*Nf);
      MfRzt1=NwfT*(Rzt1.*Nf);
      MfRxt2=NwfT*(Rxt2.*Nf);
      MfRyt2=NwfT*(Ryt2.*Nf);
      MfRzt2=NwfT*(Rzt2.*Nf);
    end
    if nsd==2
      MfRxt1Rxt1=NwfT*(Rxt1.*Rxt1.*Nf);
      MfRyt1Ryt1=NwfT*(Ryt1.*Ryt1.*Nf);
    elseif nsd==3
      MfRxt1Rxt1=NwfT*(Rxt1.*Rxt1.*Nf);
      MfRyt1Ryt1=NwfT*(Ryt1.*Ryt1.*Nf);
      MfRzt1Rzt1=NwfT*(Rzt1.*Rzt1.*Nf);
      MfRxt1Rxt2=NwfT*(Rxt1.*Rxt2.*Nf);
      MfRyt1Ryt2=NwfT*(Ryt1.*Ryt2.*Nf);
      MfRzt1Rzt2=NwfT*(Rzt1.*Rzt2.*Nf);
      MfRxt2Rxt2=NwfT*(Rxt2.*Rxt2.*Nf);
      MfRyt2Ryt2=NwfT*(Ryt2.*Ryt2.*Nf);
      MfRzt2Rzt2=NwfT*(Rzt2.*Rzt2.*Nf);
    end
    if nsd==2
      MfnxRxt1=NwfT*(nx.*Rxt1.*Nf);
      MfnxRyt1=NwfT*(nx.*Ryt1.*Nf);
      MfnyRxt1=NwfT*(ny.*Rxt1.*Nf);
      MfnyRyt1=NwfT*(ny.*Ryt1.*Nf);
    elseif nsd==3
      MfnxRxt1=NwfT*(nx.*Rxt1.*Nf);
      MfnxRyt1=NwfT*(nx.*Ryt1.*Nf);
      MfnxRzt1=NwfT*(nx.*Rzt1.*Nf);
      MfnyRxt1=NwfT*(ny.*Rxt1.*Nf);
      MfnyRyt1=NwfT*(ny.*Ryt1.*Nf);
      MfnyRzt1=NwfT*(ny.*Rzt1.*Nf);
      MfnzRxt1=NwfT*(nz.*Rxt1.*Nf);
      MfnzRyt1=NwfT*(nz.*Ryt1.*Nf);
      MfnzRzt1=NwfT*(nz.*Rzt1.*Nf);
      MfnxRxt2=NwfT*(nx.*Rxt2.*Nf);
      MfnxRyt2=NwfT*(nx.*Ryt2.*Nf);
      MfnxRzt2=NwfT*(nx.*Rzt2.*Nf);
      MfnyRxt2=NwfT*(ny.*Rxt2.*Nf);
      MfnyRyt2=NwfT*(ny.*Ryt2.*Nf);
      MfnyRzt2=NwfT*(ny.*Rzt2.*Nf);
      MfnzRxt2=NwfT*(nz.*Rxt2.*Nf);
      MfnzRyt2=NwfT*(nz.*Ryt2.*Nf);
      MfnzRzt2=NwfT*(nz.*Rzt2.*Nf);
    end
    
    % Compute lhs
    if nsd==2
      KLV(nf1,nefV1)=KLV(nf1,nefV1)-VoigtV1*Mfnx;
      KLV(nf2,nefV1)=KLV(nf2,nefV1)-VoigtV2*Mfnx;
      KLV(nf3,nefV1)=KLV(nf3,nefV1)-VoigtV3*Mfny;
      KLV(nf1,nefV2)=KLV(nf1,nefV2)-VoigtV2*Mfny;
      KLV(nf2,nefV2)=KLV(nf2,nefV2)-VoigtV1*Mfny;
      KLV(nf3,nefV2)=KLV(nf3,nefV2)-VoigtV3*Mfnx;
    elseif nsd==3
      KLV(nf1,nefV1)=KLV(nf1,nefV1)-VoigtV1*Mfnx;
      KLV(nf2,nefV1)=KLV(nf2,nefV1)-VoigtV2*Mfnx;
      KLV(nf3,nefV1)=KLV(nf3,nefV1)-VoigtV2*Mfnx;
      KLV(nf4,nefV1)=KLV(nf4,nefV1)-VoigtV3*Mfny;
      KLV(nf5,nefV1)=KLV(nf5,nefV1)-VoigtV3*Mfnz;
      KLV(nf1,nefV2)=KLV(nf1,nefV2)-VoigtV2*Mfny;
      KLV(nf2,nefV2)=KLV(nf2,nefV2)-VoigtV1*Mfny;
      KLV(nf3,nefV2)=KLV(nf3,nefV2)-VoigtV2*Mfny;
      KLV(nf4,nefV2)=KLV(nf4,nefV2)-VoigtV3*Mfnx;
      KLV(nf6,nefV2)=KLV(nf6,nefV2)-VoigtV3*Mfnz;
      KLV(nf1,nefV3)=KLV(nf1,nefV3)-VoigtV2*Mfnz;
      KLV(nf2,nefV3)=KLV(nf2,nefV3)-VoigtV2*Mfnz;
      KLV(nf3,nefV3)=KLV(nf3,nefV3)-VoigtV1*Mfnz;
      KLV(nf5,nefV3)=KLV(nf5,nefV3)-VoigtV3*Mfnx;
      KLV(nf6,nefV3)=KLV(nf6,nefV3)-VoigtV3*Mfny;
    end
    
    Kvv(nf1,nf1)=Kvv(nf1,nf1)+tauV*Mf;
    Kvv(nf2,nf2)=Kvv(nf2,nf2)+tauV*Mf;
    if nsd==3
      Kvv(nf3,nf3)=Kvv(nf3,nf3)+tauV*Mf;
    end
    
    if isConvectiveFlow
      if nsd==2
        KvV(nf1,nefV1)=KvV(nf1,nefV1)+NwfT*((2*rfg.*Vxfg.*nx+rfg.*Vyfg.*ny).*Nf);
        KvV(nf2,nefV1)=KvV(nf2,nefV1)+NwfT*((rfg.*Vyfg.*nx).*Nf);
        KvV(nf1,nefV2)=KvV(nf1,nefV2)+NwfT*((rfg.*Vxfg.*ny).*Nf);
        KvV(nf2,nefV2)=KvV(nf2,nefV2)+NwfT*((rfg.*Vxfg.*nx+2*rfg.*Vyfg.*ny).*Nf);
      elseif nsd==3
        KvV(nf1,nefV1)=KvV(nf1,nefV1)+NwfT*((2*rfg.*Vxfg.*nx+rfg.*Vyfg.*ny+rfg.*Vzfg.*nz).*Nf);
        KvV(nf2,nefV1)=KvV(nf2,nefV1)+NwfT*((rfg.*Vyfg.*nx).*Nf);
        KvV(nf3,nefV1)=KvV(nf3,nefV1)+NwfT*((rfg.*Vzfg.*nx).*Nf);
        KvV(nf1,nefV2)=KvV(nf1,nefV2)+NwfT*((rfg.*Vxfg.*ny).*Nf);
        KvV(nf2,nefV2)=KvV(nf2,nefV2)+NwfT*((rfg.*Vxfg.*nx+2*rfg.*Vyfg.*ny+rfg.*Vzfg.*nz).*Nf);
        KvV(nf3,nefV2)=KvV(nf3,nefV2)+NwfT*((rfg.*Vzfg.*ny).*Nf);
        KvV(nf1,nefV3)=KvV(nf1,nefV3)+NwfT*((rfg.*Vxfg.*nz).*Nf);
        KvV(nf2,nefV3)=KvV(nf2,nefV3)+NwfT*((rfg.*Vyfg.*nz).*Nf);
        KvV(nf3,nefV3)=KvV(nf3,nefV3)+NwfT*((rfg.*Vxfg.*nx+rfg.*Vyfg.*ny+2*rfg.*Vzfg.*nz).*Nf);
      end
      
      if nsd==2
        KvP(nf1,nefP1)=KvP(nf1,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny).*Vxfg).*Nf);
        KvP(nf2,nefP1)=KvP(nf2,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny).*Vyfg).*Nf);
      elseif nsd==3
        KvP(nf1,nefP1)=KvP(nf1,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny+Vzfg.*nz).*Vxfg).*Nf);
        KvP(nf2,nefP1)=KvP(nf2,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny+Vzfg.*nz).*Vyfg).*Nf);
        KvP(nf3,nefP1)=KvP(nf3,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny+Vzfg.*nz).*Vzfg).*Nf);
      end
    end
    
    KvV(nf1,nefV1)=KvV(nf1,nefV1)-tauV*Mf;
    KvV(nf2,nefV2)=KvV(nf2,nefV2)-tauV*Mf;
    if nsd==3
      KvV(nf3,nefV3)=KvV(nf3,nefV3)-tauV*Mf;
    end
    
    Kpp(nf1,nf1)=Kpp(nf1,nf1)+tauP*Mf;
    
    KpV(nf1,nefV1)=KpV(nf1,nefV1)+NwfT*((rfg.*nx).*Nf);
    KpV(nf1,nefV2)=KpV(nf1,nefV2)+NwfT*((rfg.*ny).*Nf);
    if nsd==3
      KpV(nf1,nefV3)=KpV(nf1,nefV3)+NwfT*((rfg.*nz).*Nf);
    end
    
    if nsd==2
      KpP(nf1,nefP1)=KpP(nf1,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny)).*Nf);
    elseif nsd==3
      KpP(nf1,nefP1)=KpP(nf1,nefP1)+NwfT*((drdPfg.*(Vxfg.*nx+Vyfg.*ny+Vzfg.*nz)).*Nf);
    end
    
    KpP(nf1,nefP1)=KpP(nf1,nefP1)-tauP*Mf;
    
    if nsd==2
      KJB(nf1,nefB1)=KJB(nf1,nefB1)-sqrt(eta)*MfnyRxt1...
                                   +sqrt(eta)*MfnxRyt1;
    elseif nsd==3
      KJB(nf1,nefB1)=KJB(nf1,nefB1)+sqrt(eta)*MfnyRzt1...
                                   -sqrt(eta)*MfnzRyt1;
      KJB(nf2,nefB1)=KJB(nf2,nefB1)+sqrt(eta)*MfnzRxt1...
                                   -sqrt(eta)*MfnxRzt1;
      KJB(nf3,nefB1)=KJB(nf3,nefB1)+sqrt(eta)*MfnxRyt1...
                                   -sqrt(eta)*MfnyRxt1;
      KJB(nf1,nefB2)=KJB(nf1,nefB2)+sqrt(eta)*MfnyRzt2...
                                   -sqrt(eta)*MfnzRyt2;
      KJB(nf2,nefB2)=KJB(nf2,nefB2)+sqrt(eta)*MfnzRxt2...
                                   -sqrt(eta)*MfnxRzt2;
      KJB(nf3,nefB2)=KJB(nf3,nefB2)+sqrt(eta)*MfnxRyt2...
                                   -sqrt(eta)*MfnyRxt2;
    end
    
    if nsd==2
      KbJ(nf1,nf1)=KbJ(nf1,nf1)+sqrt(eta)*Mfny;
      KbJ(nf2,nf1)=KbJ(nf2,nf1)-sqrt(eta)*Mfnx;
    elseif nsd==3
      KbJ(nf2,nf1)=KbJ(nf2,nf1)+sqrt(eta)*Mfnz;
      KbJ(nf3,nf1)=KbJ(nf3,nf1)-sqrt(eta)*Mfny;
      KbJ(nf1,nf2)=KbJ(nf1,nf2)-sqrt(eta)*Mfnz;
      KbJ(nf3,nf2)=KbJ(nf3,nf2)+sqrt(eta)*Mfnx;
      KbJ(nf1,nf3)=KbJ(nf1,nf3)+sqrt(eta)*Mfny;
      KbJ(nf2,nf3)=KbJ(nf2,nf3)-sqrt(eta)*Mfnx;
    end
    
    if nsd==2
      Kbb(nf1,nf1)=Kbb(nf1,nf1)+tauB*Mfnyny;
      Kbb(nf2,nf1)=Kbb(nf2,nf1)-tauB*Mfnxny;
      Kbb(nf1,nf2)=Kbb(nf1,nf2)-tauB*Mfnxny;
      Kbb(nf2,nf2)=Kbb(nf2,nf2)+tauB*Mfnxnx;
    elseif nsd==3
      Kbb(nf1,nf1)=Kbb(nf1,nf1)+tauB*Mfnyny+tauB*Mfnznz;
      Kbb(nf2,nf1)=Kbb(nf2,nf1)-tauB*Mfnxny;
      Kbb(nf3,nf1)=Kbb(nf3,nf1)-tauB*Mfnxnz;
      Kbb(nf1,nf2)=Kbb(nf1,nf2)-tauB*Mfnxny;
      Kbb(nf2,nf2)=Kbb(nf2,nf2)+tauB*Mfnxnx+tauB*Mfnznz;
      Kbb(nf3,nf2)=Kbb(nf3,nf2)-tauB*Mfnynz;
      Kbb(nf1,nf3)=Kbb(nf1,nf3)-tauB*Mfnxnz;
      Kbb(nf2,nf3)=Kbb(nf2,nf3)-tauB*Mfnynz;
      Kbb(nf3,nf3)=Kbb(nf3,nf3)+tauB*Mfnxnx+tauB*Mfnyny;
    end
    
    if isCouplingTerms
      if nsd==2
        Kbb(nf1,nf1)=Kbb(nf1,nf1)+NwfT*((+Vyfg.*ny).*Nf);
        Kbb(nf2,nf1)=Kbb(nf2,nf1)+NwfT*((-Vyfg.*nx).*Nf);
        Kbb(nf1,nf2)=Kbb(nf1,nf2)+NwfT*((-Vxfg.*ny).*Nf);
        Kbb(nf2,nf2)=Kbb(nf2,nf2)+NwfT*((+Vxfg.*nx).*Nf);
      elseif nsd==3
        Kbb(nf1,nf1)=Kbb(nf1,nf1)+NwfT*((+Vyfg.*ny+Vzfg.*nz).*Nf);
        Kbb(nf2,nf1)=Kbb(nf2,nf1)+NwfT*((-Vyfg.*nx).*Nf);
        Kbb(nf3,nf1)=Kbb(nf3,nf1)+NwfT*((-Vzfg.*nx).*Nf);
        Kbb(nf1,nf2)=Kbb(nf1,nf2)+NwfT*((-Vxfg.*ny).*Nf);
        Kbb(nf2,nf2)=Kbb(nf2,nf2)+NwfT*((+Vxfg.*nx+Vzfg.*nz).*Nf);
        Kbb(nf3,nf2)=Kbb(nf3,nf2)+NwfT*((-Vzfg.*ny).*Nf);
        Kbb(nf1,nf3)=Kbb(nf1,nf3)+NwfT*((-Vxfg.*nz).*Nf);
        Kbb(nf2,nf3)=Kbb(nf2,nf3)+NwfT*((-Vyfg.*nz).*Nf);
        Kbb(nf3,nf3)=Kbb(nf3,nf3)+NwfT*((+Vxfg.*nx+Vyfg.*ny).*Nf);
      end
      
      if nsd==2
        KbV(nf1,nefV1)=KbV(nf1,nefV1)+NwfT*((-byfg.*ny).*Nf);
        KbV(nf2,nefV1)=KbV(nf2,nefV1)+NwfT*((+byfg.*nx).*Nf);
        KbV(nf1,nefV2)=KbV(nf1,nefV2)+NwfT*((+bxfg.*ny).*Nf);
        KbV(nf2,nefV2)=KbV(nf2,nefV2)+NwfT*((-bxfg.*nx).*Nf);
      elseif nsd==3
        KbV(nf1,nefV1)=KbV(nf1,nefV1)+NwfT*((-byfg.*ny-bzfg.*nz).*Nf);
        KbV(nf2,nefV1)=KbV(nf2,nefV1)+NwfT*((+byfg.*nx).*Nf);
        KbV(nf3,nefV1)=KbV(nf3,nefV1)+NwfT*((+bzfg.*nx).*Nf);
        KbV(nf1,nefV2)=KbV(nf1,nefV2)+NwfT*((+bxfg.*ny).*Nf);
        KbV(nf2,nefV2)=KbV(nf2,nefV2)+NwfT*((-bxfg.*nx-bzfg.*nz).*Nf);
        KbV(nf3,nefV2)=KbV(nf3,nefV2)+NwfT*((+bzfg.*ny).*Nf);
        KbV(nf1,nefV3)=KbV(nf1,nefV3)+NwfT*((+bxfg.*nz).*Nf);
        KbV(nf2,nefV3)=KbV(nf2,nefV3)+NwfT*((+byfg.*nz).*Nf);
        KbV(nf3,nefV3)=KbV(nf3,nefV3)+NwfT*((-bxfg.*nx-byfg.*ny).*Nf);
      end
    end
    
    if nsd==2
      KbB(nf1,nefB1)=KbB(nf1,nefB1)-tauB*MfRxt1;
      KbB(nf2,nefB1)=KbB(nf2,nefB1)-tauB*MfRyt1;
    elseif nsd==3
      KbB(nf1,nefB1)=KbB(nf1,nefB1)-tauB*MfRxt1;
      KbB(nf2,nefB1)=KbB(nf2,nefB1)-tauB*MfRyt1;
      KbB(nf3,nefB1)=KbB(nf3,nefB1)-tauB*MfRzt1;
      KbB(nf1,nefB2)=KbB(nf1,nefB2)-tauB*MfRxt2;
      KbB(nf2,nefB2)=KbB(nf2,nefB2)-tauB*MfRyt2;
      KbB(nf3,nefB2)=KbB(nf3,nefB2)-tauB*MfRzt2;
    end
    
    KbQ(nf1,nefQ1)=KbQ(nf1,nefQ1)+Mfnx;
    KbQ(nf2,nefQ1)=KbQ(nf2,nefQ1)+Mfny;
    if nsd==3
      KbQ(nf3,nefQ1)=KbQ(nf3,nefQ1)+Mfnz;
    end
    
    Kqb(nf1,nf1)=Kqb(nf1,nf1)-Mfnx;
    Kqb(nf1,nf2)=Kqb(nf1,nf2)-Mfny;
    if nsd==3
      Kqb(nf1,nf3)=Kqb(nf1,nf3)-Mfnz;
    end
    
    Kqq(nf1,nf1)=Kqq(nf1,nf1)-tauQ*Mf;
    
    KqQ(nf1,nefQ1)=KqQ(nf1,nefQ1)+tauQ*Mf;
    
    if not(isExterior) || isNeumann_t_x
      KVp(nefV1,nf1)=KVp(nefV1,nf1)-Mfnx;
    end
    
    if not(isExterior) || isNeumann_t_y
      KVp(nefV2,nf1)=KVp(nefV2,nf1)-Mfny;
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z)
      KVp(nefV2,nf1)=KVp(nefV2,nf1)-Mfnz;
    end
    
    if not(isDirichlet_v_x)
      KVv(nefV1,nf1)=KVv(nefV1,nf1)-tauV*Mf;
    end
    
    if not(isDirichlet_v_y)
      KVv(nefV2,nf2)=KVv(nefV2,nf2)-tauV*Mf;
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      KVv(nefV3,nf3)=KVv(nefV3,nf3)-tauV*Mf;
    end
    
    KVV(nefV1,nefV1)=KVV(nefV1,nefV1)+tauV*Mf;
    KVV(nefV2,nefV2)=KVV(nefV2,nefV2)+tauV*Mf;
    if nsd==3
      KVV(nefV3,nefV3)=KVV(nefV3,nefV3)+tauV*Mf;
    end
    
    if not(isDirichlet_p)
      KPp(nefP1,nf1)=KPp(nefP1,nf1)+tauP*Mf;
    end
    
    KPP(nefP1,nefP1)=KPP(nefP1,nefP1)-tauP*Mf;
    
    if isCouplingTerms && not(isDirichlet_m)
      if nsd==2
        KBb(nefB1,nf1)=KBb(nefB1,nf1)-NwfT*((+Rxt1.*Vyfg.*ny...
                                             -Ryt1.*Vyfg.*nx).*Nf);
        KBb(nefB1,nf2)=KBb(nefB1,nf2)-NwfT*((-Rxt1.*Vxfg.*ny...
                                             +Ryt1.*Vxfg.*nx).*Nf);
      elseif nsd==3
        KBb(nefB1,nf1)=KBb(nefB1,nf1)-NwfT*((+Rxt1.*(Vyfg.*ny+Vzfg.*nz)...
                                             -Ryt1.*(Vyfg.*nx)...
                                             -Rzt1.*(Vzfg.*nx)).*Nf);
        KBb(nefB2,nf1)=KBb(nefB2,nf1)-NwfT*((+Rxt2.*(Vyfg.*ny+Vzfg.*nz)...
                                             -Ryt2.*(Vyfg.*nx)...
                                             -Rzt2.*(Vzfg.*nx)).*Nf);
        KBb(nefB1,nf2)=KBb(nefB1,nf2)-NwfT*((-Rxt1.*(Vxfg.*ny)...
                                             +Ryt1.*(Vxfg.*nx+Vzfg.*nz)...
                                             -Rzt1.*(Vzfg.*ny)).*Nf);
        KBb(nefB2,nf2)=KBb(nefB2,nf2)-NwfT*((-Rxt2.*(Vxfg.*ny)...
                                             +Ryt2.*(Vxfg.*nx+Vzfg.*nz)...
                                             -Rzt2.*(Vzfg.*ny)).*Nf);
        KBb(nefB1,nf3)=KBb(nefB1,nf3)-NwfT*((-Rxt1.*(Vxfg.*nz)...
                                             -Ryt1.*(Vyfg.*nz)...
                                             +Rzt1.*(Vxfg.*nx+Vyfg.*ny)).*Nf);
        KBb(nefB2,nf3)=KBb(nefB2,nf3)-NwfT*((-Rxt2.*(Vxfg.*nz)...
                                             -Ryt2.*(Vyfg.*nz)...
                                             +Rzt2.*(Vxfg.*nx+Vyfg.*ny)).*Nf);
      end
      
      if nsd==2
        KBV(nefB1,nefV1)=KBV(nefB1,nefV1)-NwfT*((-Rxt1.*byfg.*ny...
                                                 +Ryt1.*byfg.*nx).*Nf);
        KBV(nefB1,nefV2)=KBV(nefB1,nefV2)-NwfT*((+Rxt1.*bxfg.*ny...
                                                 -Ryt1.*bxfg.*nx).*Nf);
      elseif nsd==3
        KBV(nefB1,nefV1)=KBV(nefB1,nefV1)-NwfT*((-Rxt1.*(byfg.*ny+bzfg.*nz)...
                                                 +Ryt1.*(byfg.*nx)...
                                                 +Rzt1.*(bzfg.*nx)).*Nf);
        KBV(nefB2,nefV1)=KBV(nefB2,nefV1)-NwfT*((-Rxt2.*(byfg.*ny+bzfg.*nz)...
                                                 +Ryt2.*(byfg.*nx)...
                                                 +Rzt2.*(bzfg.*nx)).*Nf);
        KBV(nefB1,nefV2)=KBV(nefB1,nefV2)-NwfT*((+Rxt1.*(bxfg.*ny)...
                                                 -Ryt1.*(bxfg.*nx+bzfg.*nz)...
                                                 +Rzt1.*(bzfg.*ny)).*Nf);
        KBV(nefB2,nefV2)=KBV(nefB2,nefV2)-NwfT*((+Rxt2.*(bxfg.*ny)...
                                                 -Ryt2.*(bxfg.*nx+bzfg.*nz)...
                                                 +Rzt2.*(bzfg.*ny)).*Nf);
        KBV(nefB1,nefV3)=KBV(nefB1,nefV3)-NwfT*((+Rxt1.*(bxfg.*nz)...
                                                 +Ryt1.*(byfg.*nz)...
                                                 -Rzt1.*(bxfg.*nx+byfg.*ny)).*Nf);
        KBV(nefB2,nefV3)=KBV(nefB2,nefV3)-NwfT*((+Rxt2.*(bxfg.*nz)...
                                                 +Ryt2.*(byfg.*nz)...
                                                 -Rzt2.*(bxfg.*nx+byfg.*ny)).*Nf);
      end
    end
    
    if not(isDirichlet_m)
      if nsd==2
        KBQ(nefB1,nefQ1)=KBQ(nefB1,nefQ1)-MfnxRxt1...
                                         -MfnyRyt1;
      elseif nsd==3
        KBQ(nefB1,nefQ1)=KBQ(nefB1,nefQ1)-MfnxRxt1...
                                         -MfnyRyt1...
                                         -MfnzRzt1;
        KBQ(nefB2,nefQ1)=KBQ(nefB2,nefQ1)-MfnxRxt2...
                                         -MfnyRyt2...
                                         -MfnzRzt2;
      end
    end
    
    if not(isDirichlet_m)
      if nsd==2
        KBb(nefB1,nf1)=KBb(nefB1,nf1)-tauB*MfRxt1;
        KBb(nefB1,nf2)=KBb(nefB1,nf2)-tauB*MfRyt1;
      elseif nsd==3
        KBb(nefB1,nf1)=KBb(nefB1,nf1)-tauB*MfRxt1;
        KBb(nefB2,nf1)=KBb(nefB2,nf1)-tauB*MfRxt2;
        KBb(nefB1,nf2)=KBb(nefB1,nf2)-tauB*MfRyt1;
        KBb(nefB2,nf2)=KBb(nefB2,nf2)-tauB*MfRyt2;
        KBb(nefB1,nf3)=KBb(nefB1,nf3)-tauB*MfRzt1;
        KBb(nefB2,nf3)=KBb(nefB2,nf3)-tauB*MfRzt2;
      end
    end
    
    if nsd==2
      KBB(nefB1,nefB1)=KBB(nefB1,nefB1)+tauB*MfRxt1Rxt1...
                                       +tauB*MfRyt1Ryt1;
    elseif nsd==3
      KBB(nefB1,nefB1)=KBB(nefB1,nefB1)+tauB*MfRxt1Rxt1...
                                       +tauB*MfRyt1Ryt1...
                                       +tauB*MfRzt1Rzt1;
      KBB(nefB2,nefB1)=KBB(nefB2,nefB1)+tauB*MfRxt1Rxt2...
                                       +tauB*MfRyt1Ryt2...
                                       +tauB*MfRzt1Rzt2;
      KBB(nefB1,nefB2)=KBB(nefB1,nefB2)+tauB*MfRxt1Rxt2...
                                       +tauB*MfRyt1Ryt2...
                                       +tauB*MfRzt1Rzt2;
      KBB(nefB2,nefB2)=KBB(nefB2,nefB2)+tauB*MfRxt2Rxt2...
                                       +tauB*MfRyt2Ryt2...
                                       +tauB*MfRzt2Rzt2;
    end
    
    if not(isDirichlet_m)
      KQb(nefQ1,nf1)=KQb(nefQ1,nf1)+Mfnx;
      KQb(nefQ1,nf2)=KQb(nefQ1,nf2)+Mfny;
      if nsd==3
        KQb(nefQ1,nf3)=KQb(nefQ1,nf3)+Mfnz;
      end
    end
    
    KQQ(nefQ1,nefQ1)=KQQ(nefQ1,nefQ1)-tauQ*Mf;
    
    % Compute rhs
    if nsd==2
      fL(nf1,1)=fL(nf1,1)+NwfT*(+VoigtV1*nx.*Vxfg...
                                +VoigtV2*ny.*Vyfg);
      fL(nf2,1)=fL(nf2,1)+NwfT*(+VoigtV2*nx.*Vxfg...
                                +VoigtV1*ny.*Vyfg);
      fL(nf3,1)=fL(nf3,1)+NwfT*(+VoigtV3*ny.*Vxfg...
                                +VoigtV3*nx.*Vyfg);
    elseif nsd==3
      fL(nf1,1)=fL(nf1,1)+NwfT*(+VoigtV1*nx.*Vxfg...
                                +VoigtV2*ny.*Vyfg...
                                +VoigtV2*nz.*Vzfg);
      fL(nf2,1)=fL(nf2,1)+NwfT*(+VoigtV2*nx.*Vxfg...
                                +VoigtV1*ny.*Vyfg...
                                +VoigtV2*nz.*Vzfg);
      fL(nf3,1)=fL(nf3,1)+NwfT*(+VoigtV2*nx.*Vxfg...
                                +VoigtV2*ny.*Vyfg...
                                +VoigtV1*nz.*Vzfg);
      fL(nf4,1)=fL(nf4,1)+NwfT*(+VoigtV3*nx.*Vyfg...
                                +VoigtV3*ny.*Vxfg);
      fL(nf5,1)=fL(nf5,1)+NwfT*(+VoigtV3*nx.*Vzfg...
                                +VoigtV3*nz.*Vxfg);
      fL(nf6,1)=fL(nf6,1)+NwfT*(+VoigtV3*ny.*Vzfg...
                                +VoigtV3*nz.*Vyfg);
    end
    
    fv(nf1,1)=fv(nf1,1)-NwfT*(tauV*vxfg);
    fv(nf2,1)=fv(nf2,1)-NwfT*(tauV*vyfg);
    if nsd==3
      fv(nf3,1)=fv(nf3,1)-NwfT*(tauV*vzfg);
    end
    
    if isConvectiveFlow
      if nsd==2
        fv(nf1,1)=fv(nf1,1)-NwfT*((Vxfg.*nx+Vyfg.*ny).*r(Pfg).*Vxfg);
        fv(nf2,1)=fv(nf2,1)-NwfT*((Vxfg.*nx+Vyfg.*ny).*r(Pfg).*Vyfg);
      elseif nsd==3
        fv(nf1,1)=fv(nf1,1)-NwfT*((Vxfg.*nx+Vyfg.*ny+Vzfg.*nz).*r(Pfg).*Vxfg);
        fv(nf2,1)=fv(nf2,1)-NwfT*((Vxfg.*nx+Vyfg.*ny+Vzfg.*nz).*r(Pfg).*Vyfg);
        fv(nf3,1)=fv(nf3,1)-NwfT*((Vxfg.*nx+Vyfg.*ny+Vzfg.*nz).*r(Pfg).*Vzfg);
      end
    end
    
    fv(nf1,1)=fv(nf1,1)+NwfT*(tauV*Vxfg);
    fv(nf2,1)=fv(nf2,1)+NwfT*(tauV*Vyfg);
    if nsd==3
      fv(nf3,1)=fv(nf3,1)+NwfT*(tauV*Vzfg);
    end
    
    fp(nf1,1)=fp(nf1,1)-NwfT*(tauP*pfg);
    
    fp(nf1,1)=fp(nf1,1)-NwfT*(r(Pfg).*(Vxfg.*nx+Vyfg.*ny)-tauP*Pfg);
    if nsd==3
      fp(nf1,1)=fp(nf1,1)-NwfT*(r(Pfg).*Vzfg.*nz);
    end
    
    if nsd==2
      fJ(nf1,1)=fJ(nf1,1)-NwfT*(sqrt(eta)*(nx.*Btyfg-ny.*Btxfg));
    elseif nsd==3
      fJ(nf1,1)=fJ(nf1,1)-NwfT*(sqrt(eta)*(ny.*Btzfg-nz.*Btyfg));
      fJ(nf2,1)=fJ(nf2,1)-NwfT*(sqrt(eta)*(nz.*Btxfg-nx.*Btzfg));
      fJ(nf3,1)=fJ(nf3,1)-NwfT*(sqrt(eta)*(nx.*Btyfg-ny.*Btxfg));
    end
    
    if nsd==2
      fb(nf1,1)=fb(nf1,1)+NwfT*(sqrt(eta)*(-ny.*Jzfg));
      fb(nf2,1)=fb(nf2,1)+NwfT*(sqrt(eta)*(+nx.*Jzfg));
    elseif nsd==3
      fb(nf1,1)=fb(nf1,1)+NwfT*(sqrt(eta)*(+nz.*Jyfg-ny.*Jzfg));
      fb(nf2,1)=fb(nf2,1)+NwfT*(sqrt(eta)*(+nx.*Jzfg-nz.*Jxfg));
      fb(nf3,1)=fb(nf3,1)+NwfT*(sqrt(eta)*(+ny.*Jxfg-nx.*Jyfg));
    end
    
    if nsd==2
      fb(nf1,1)=fb(nf1,1)+NwfT*(tauB*(+ny.*(nx.*byfg-ny.*bxfg)));
      fb(nf2,1)=fb(nf2,1)+NwfT*(tauB*(-nx.*(nx.*byfg-ny.*bxfg)));
    elseif nsd==3
      fb(nf1,1)=fb(nf1,1)+NwfT*(tauB*(+ny.*(nx.*byfg-ny.*bxfg)...
                                      -nz.*(nz.*bxfg-nx.*bzfg)));
      fb(nf2,1)=fb(nf2,1)+NwfT*(tauB*(+nz.*(ny.*bzfg-nz.*byfg)...
                                      -nx.*(nx.*byfg-ny.*bxfg)));
      fb(nf3,1)=fb(nf3,1)+NwfT*(tauB*(+nx.*(nz.*bxfg-nx.*bzfg)...
                                      -ny.*(ny.*bzfg-nz.*byfg)));
    end
    
    if isCouplingTerms
      if nsd==2
        fb(nf1,1)=fb(nf1,1)-NwfT*((bxfg.*Vyfg-Vxfg.*byfg).*ny);
        fb(nf2,1)=fb(nf2,1)-NwfT*((byfg.*Vxfg-Vyfg.*bxfg).*nx);
      elseif nsd==3
        fb(nf1,1)=fb(nf1,1)-NwfT*(((bxfg.*Vyfg-Vxfg.*byfg).*ny...
                                  +(bxfg.*Vzfg-Vxfg.*bzfg).*nz));
        fb(nf2,1)=fb(nf2,1)-NwfT*(((byfg.*Vxfg-Vyfg.*bxfg).*nx...
                                  +(byfg.*Vzfg-Vyfg.*bzfg).*nz));
        fb(nf3,1)=fb(nf3,1)-NwfT*(((bzfg.*Vxfg-Vzfg.*bxfg).*nx...
                                  +(bzfg.*Vyfg-Vzfg.*byfg).*ny));
      end
    end
    
    fb(nf1,1)=fb(nf1,1)+NwfT*(tauB*Btxfg);
    fb(nf2,1)=fb(nf2,1)+NwfT*(tauB*Btyfg);
    if nsd==3
      fb(nf3,1)=fb(nf3,1)+NwfT*(tauB*Btzfg);
    end
    
    fb(nf1,1)=fb(nf1,1)-NwfT*(Qfg.*nx);
    fb(nf2,1)=fb(nf2,1)-NwfT*(Qfg.*ny);
    if nsd==3
      fb(nf3,1)=fb(nf3,1)-NwfT*(Qfg.*nz);
    end
    
    if nsd==2
      fq(nf1,1)=fq(nf1,1)+NwfT*(bxfg.*nx+byfg.*ny);
    elseif nsd==3
      fq(nf1,1)=fq(nf1,1)+NwfT*(bxfg.*nx+byfg.*ny+bzfg.*nz);
    end
    
    fq(nf1,1)=fq(nf1,1)+NwfT*(tauQ*qfg);
    
    fq(nf1,1)=fq(nf1,1)-NwfT*(tauQ*Qfg);
    
    if not(isExterior) || isNeumann_t_x
      if nsd==2
        fV(nefV1,1)=fV(nefV1,1)+NwfT*(+VoigtV1*nx.*Lxxfg+VoigtV2*nx.*Lyyfg...
                                      +VoigtV3*ny.*Lxyfg...
                                      +pfg.*nx);
      elseif nsd==3
        fV(nefV1,1)=fV(nefV1,1)+NwfT*(+VoigtV1*nx.*Lxxfg+VoigtV2*nx.*Lyyfg+VoigtV2*nx.*Lzzfg...
                                      +VoigtV3*ny.*Lxyfg+VoigtV3*nz.*Lxzfg...
                                      +pfg.*nx);
      end
    end
    
    if not(isExterior) || isNeumann_t_y
      if nsd==2
        fV(nefV2,1)=fV(nefV2,1)+NwfT*(+VoigtV2*ny.*Lxxfg+VoigtV1*ny.*Lyyfg...
                                      +VoigtV3*nx.*Lxyfg...
                                      +pfg.*ny);
      elseif nsd==3
        fV(nefV2,1)=fV(nefV2,1)+NwfT*(+VoigtV2*ny.*Lxxfg+VoigtV1*ny.*Lyyfg+VoigtV2*ny.*Lzzfg...
                                      +VoigtV3*nx.*Lxyfg+VoigtV3*nz.*Lyzfg...
                                      +pfg.*ny);
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z)
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(+VoigtV2*nz.*Lxxfg+VoigtV2*nz.*Lyyfg+VoigtV1*nz.*Lzzfg...
                                    +VoigtV3*nx.*Lxzfg+VoigtV3*ny.*Lyzfg...
                                    +pfg.*nz);
    end
    
    if not(isDirichlet_v_x)
      fV(nefV1,1)=fV(nefV1,1)+NwfT*(+tauV*vxfg);
    end
    
    if not(isDirichlet_v_y)
      fV(nefV2,1)=fV(nefV2,1)+NwfT*(+tauV*vyfg);
    end
    
    if  nsd==3 && not(isDirichlet_v_z)
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(+tauV*vzfg);
    end
    
    if nsd==2
      fV(nefV1,1)=fV(nefV1,1)-NwfT*(tauV*Vxfg);
      fV(nefV2,1)=fV(nefV2,1)-NwfT*(tauV*Vyfg);
    elseif nsd==3
      fV(nefV1,1)=fV(nefV1,1)-NwfT*(tauV*Vxfg);
      fV(nefV2,1)=fV(nefV2,1)-NwfT*(tauV*Vyfg);
      fV(nefV3,1)=fV(nefV3,1)-NwfT*(tauV*Vzfg);
    end
    
    if isDirichlet_v_x
      fV(nefV1,1)=fV(nefV1,1)+NwfT*(tauV*vDxfg);
    end
    
    if isDirichlet_v_y
      fV(nefV2,1)=fV(nefV2,1)+NwfT*(tauV*vDyfg);
    end
    
    if nsd==3 && isDirichlet_v_z
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(tauV*vDzfg);
    end
    
    if isNeumann_t_x
      fV(nefV1,1)=fV(nefV1,1)+NwfT*(tNxfg);
    end
    
    if isNeumann_t_y
      fV(nefV2,1)=fV(nefV2,1)+NwfT*(tNyfg);
    end
    
    if nsd==3 && isNeumann_t_z
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(tNzfg);
    end
    
    if not(isDirichlet_p)
      fP(nefP1,1)=fP(nefP1,1)-NwfT*(tauP*pfg);
    end
    
    fP(nefP1,1)=fP(nefP1,1)+NwfT*(tauP*Pfg);
    
    if isDirichlet_p
      fP(nefP1,1)=fP(nefP1,1)-NwfT*(tauP*pDfg);
    end
    
    if isCouplingTerms && not(isDirichlet_m)
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(bxfg.*Vyfg-Vxfg.*byfg).*ny...
                                     +Ryt1.*(byfg.*Vxfg-Vyfg.*bxfg).*nx);
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*((bxfg.*Vyfg-Vxfg.*byfg).*ny+...
                                             (bxfg.*Vzfg-Vxfg.*bzfg).*nz)...
                                     +Ryt1.*((byfg.*Vxfg-Vyfg.*bxfg).*nx+...
                                             (byfg.*Vzfg-Vyfg.*bzfg).*nz)...
                                     +Rzt1.*((bzfg.*Vxfg-Vzfg.*bxfg).*nx+...
                                             (bzfg.*Vyfg-Vzfg.*byfg).*ny));
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(Rxt2.*((bxfg.*Vyfg-Vxfg.*byfg).*ny+...
                                             (bxfg.*Vzfg-Vxfg.*bzfg).*nz)...
                                     +Ryt2.*((byfg.*Vxfg-Vyfg.*bxfg).*nx+...
                                             (byfg.*Vzfg-Vyfg.*bzfg).*nz)...
                                     +Rzt2.*((bzfg.*Vxfg-Vzfg.*bxfg).*nx+...
                                             (bzfg.*Vyfg-Vzfg.*byfg).*ny));
      end
    end
    
    if not(isDirichlet_m)
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)-NwfT*(Rxt1.*(sqrt(eta)*(-ny.*Jzfg))+...
                                      Ryt1.*(sqrt(eta)*(+nx.*Jzfg)));
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)-NwfT*(Rxt1.*(sqrt(eta)*(nz.*Jyfg-ny.*Jzfg))+...
                                      Ryt1.*(sqrt(eta)*(nx.*Jzfg-nz.*Jxfg))+...
                                      Rzt1.*(sqrt(eta)*(ny.*Jxfg-nx.*Jyfg)));
        fB(nefB2,1)=fB(nefB2,1)-NwfT*(Rxt2.*(sqrt(eta)*(nz.*Jyfg-ny.*Jzfg))+...
                                      Ryt2.*(sqrt(eta)*(nx.*Jzfg-nz.*Jxfg))+...
                                      Rzt2.*(sqrt(eta)*(ny.*Jxfg-nx.*Jyfg)));
      end
    end
    
    if not(isDirichlet_m)
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(Qfg.*nx)+...
                                      Ryt1.*(Qfg.*ny));
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(Qfg.*nx)+...
                                      Ryt1.*(Qfg.*ny)+...
                                      Rzt1.*(Qfg.*nz));
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(Rxt2.*(Qfg.*nx)+...
                                      Ryt2.*(Qfg.*ny)+...
                                      Rzt2.*(Qfg.*nz));
      end
    end
    
    if not(isDirichlet_m)
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(tauB*bxfg)+...
                                      Ryt1.*(tauB*byfg));
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(tauB*bxfg)+...
                                      Ryt1.*(tauB*byfg)+...
                                      Rzt1.*(tauB*bzfg));
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(Rxt2.*(tauB*bxfg)+...
                                      Ryt2.*(tauB*byfg)+...
                                      Rzt2.*(tauB*bzfg));
      end
    end
    
    if nsd==2
      fB(nefB1,1)=fB(nefB1,1)-NwfT*(Rxt1.*(tauB*Btxfg)+...
                                    Ryt1.*(tauB*Btyfg));
    elseif nsd==3
      fB(nefB1,1)=fB(nefB1,1)-NwfT*(Rxt1.*(tauB*Btxfg)+...
                                    Ryt1.*(tauB*Btyfg)+...
                                    Rzt1.*(tauB*Btzfg));
      fB(nefB2,1)=fB(nefB2,1)-NwfT*(Rxt2.*(tauB*Btxfg)+...
                                    Ryt2.*(tauB*Btyfg)+...
                                    Rzt2.*(tauB*Btzfg));
    end
    
    if isDirichlet_m
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(tauB*bDtxfg)+...
                                      Ryt1.*(tauB*bDtyfg));
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(tauB*bDtxfg)+...
                                      Ryt1.*(tauB*bDtyfg)+...
                                      Rzt1.*(tauB*bDtzfg));
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(Rxt2.*(tauB*bDtxfg)+...
                                      Ryt2.*(tauB*bDtyfg)+...
                                      Rzt2.*(tauB*bDtzfg));
      end
    end
    
    if isNatural_m
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*eNtxfg+...
                                      Ryt1.*eNtyfg);
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*eNtxfg+...
                                      Ryt1.*eNtyfg+...
                                      Rzt1.*eNtzfg);
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(Rxt2.*eNtxfg+...
                                      Ryt2.*eNtyfg+...
                                      Rzt2.*eNtzfg);
      end
    end
    
    if not(isDirichlet_m)
      if nsd==2
        fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(bxfg.*nx+byfg.*ny);
      elseif nsd==3
        fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(bxfg.*nx+byfg.*ny+bzfg.*nz);
      end
    end
    
    if not(isDirichlet_m)
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*qfg);
    end
    
    fQ(nefQ1,1)=fQ(nefQ1,1)+NwfT*(tauQ*Qfg);
    
    if isDirichlet_m
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*qDfg);
    end
    
    if isNatural_m
      fQ(nefQ1,1)=fQ(nefQ1,1)+NwfT*(bNnfg);
    end
    
    % Clear blocks
    if nsd==2
      KVL(nefV1,nf1)=KLV(nf1,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf1)=KLV(nf1,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV1,nf2)=KLV(nf2,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf2)=KLV(nf2,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV1,nf3)=KLV(nf3,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf3)=KLV(nf3,nefV2)'*(not(isExterior) || isNeumann_t_y);
    elseif nsd==3
      KVL(nefV1,nf1)=KLV(nf1,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf1)=KLV(nf1,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV3,nf1)=KLV(nf1,nefV3)'*(not(isExterior) || isNeumann_t_z);
      KVL(nefV1,nf2)=KLV(nf2,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf2)=KLV(nf2,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV3,nf2)=KLV(nf2,nefV3)'*(not(isExterior) || isNeumann_t_z);
      KVL(nefV1,nf3)=KLV(nf3,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf3)=KLV(nf3,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV3,nf3)=KLV(nf3,nefV3)'*(not(isExterior) || isNeumann_t_z);
      KVL(nefV1,nf4)=KLV(nf4,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV2,nf4)=KLV(nf4,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV1,nf5)=KLV(nf5,nefV1)'*(not(isExterior) || isNeumann_t_x);
      KVL(nefV3,nf5)=KLV(nf5,nefV3)'*(not(isExterior) || isNeumann_t_z);
      KVL(nefV2,nf6)=KLV(nf6,nefV2)'*(not(isExterior) || isNeumann_t_y);
      KVL(nefV3,nf6)=KLV(nf6,nefV3)'*(not(isExterior) || isNeumann_t_z);
    end
    
    if nsd==2
      KBJ(nefB1,nf1)=KJB(nf1,nefB1)'*not(isDirichlet_m);
    elseif nsd==3
      KBJ(nefB1,nf1)=KJB(nf1,nefB1)'*not(isDirichlet_m);
      KBJ(nefB1,nf2)=KJB(nf2,nefB1)'*not(isDirichlet_m);
      KBJ(nefB1,nf3)=KJB(nf3,nefB1)'*not(isDirichlet_m);
      KBJ(nefB2,nf1)=KJB(nf1,nefB2)'*not(isDirichlet_m);
      KBJ(nefB2,nf2)=KJB(nf2,nefB2)'*not(isDirichlet_m);
      KBJ(nefB2,nf3)=KJB(nf3,nefB2)'*not(isDirichlet_m);
    end
    
    KQq(nefQ1,nf1)=KqQ(nf1,nefQ1)'*not(isDirichlet_m);
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iL=1:msd*NumElementNodes;
iv=iL(end)+(1:nsd*NumElementNodes);
ip=iv(end)+(1:NumElementNodes);
iJ=ip(end)+(1:qsd*NumElementNodes);
ib=iJ(end)+(1:nsd*NumElementNodes);
iq=ib(end)+(1:NumElementNodes);
iV=reshape((0:NumElementFaces-1)*(nsd+1+nsd-1+1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iP=reshape((0:NumElementFaces-1)*(nsd+1+nsd-1+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+nsd*NumFaceNodes;
iB=reshape((0:NumElementFaces-1)*(nsd+1+nsd-1+1)*NumFaceNodes+repmat((1:(nsd-1)*NumFaceNodes)',...
  1,NumElementFaces),1,[])+(nsd+1)*NumFaceNodes;
iQ=reshape((0:NumElementFaces-1)*(nsd+1+nsd-1+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+(nsd+1+nsd-1)*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+nsd+1+qsd+nsd+1)*NumElementNodes,(msd+nsd+1+qsd+nsd+1)*NumElementNodes);
LhsLG=zeros((msd+nsd+1+qsd+nsd+1)*NumElementNodes,(nsd+1+nsd-1+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd+1+nsd-1+1)*NumElementFaces*NumFaceNodes,(msd+nsd+1+qsd+nsd+1)*NumElementNodes);
LhsGG=zeros((nsd+1+nsd-1+1)*NumElementFaces*NumFaceNodes,(nsd+1+nsd-1+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+nsd+1+qsd+nsd+1)*NumElementNodes,1);
RhsG=zeros((nsd+1+nsd-1+1)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iL,iL)=KLL;
LhsLL(iL,iv)=KLv;
LhsLL(iv,iL)=KLv';
LhsLL(iv,iv)=Kvv;
LhsLL(iv,ip)=Kvp;
LhsLL(iv,iJ)=KvJ;
LhsLL(iv,ib)=Kvb;
LhsLL(ip,iv)=Kpv;
LhsLL(ip,ip)=Kpp;
LhsLL(iJ,iJ)=KJJ;
LhsLL(iJ,ib)=KJb;
LhsLL(ib,iv)=Kbv;
LhsLL(ib,iJ)=KbJ;
LhsLL(ib,ib)=Kbb;
LhsLL(ib,iq)=Kbq;
LhsLL(iq,ib)=Kqb;
LhsLL(iq,iq)=Kqq;

% Lhs local-global
LhsLG(iL,iV)=KLV;
LhsLG(iv,iV)=KvV;
LhsLG(iv,iP)=KvP;
LhsLG(ip,iV)=KpV;
LhsLG(ip,iP)=KpP;
LhsLG(iJ,iB)=KJB;
LhsLG(ib,iV)=KbV;
LhsLG(ib,iB)=KbB;
LhsLG(ib,iQ)=KbQ;
LhsLG(iq,iQ)=KqQ;

% Rhs local
RhsL(iL,1)=fL;
RhsL(iv,1)=fv;
RhsL(ip,1)=fp;
RhsL(iJ,1)=fJ;
RhsL(ib,1)=fb;
RhsL(iq,1)=fq; 

% Lhs global-local
LhsGL(iV,iL)=KVL;
LhsGL(iV,iv)=KVv;
LhsGL(iV,ip)=KVp;
LhsGL(iP,ip)=KPp;
LhsGL(iB,iJ)=KBJ;
LhsGL(iB,ib)=KBb;
LhsGL(iQ,ib)=KQb;
LhsGL(iQ,iq)=KQq;

% Lhs global-global
LhsGG(iV,iV)=KVV;
LhsGG(iP,iP)=KPP;
LhsGG(iB,iV)=KBV;
LhsGG(iB,iB)=KBB;
LhsGG(iB,iQ)=KBQ;
LhsGG(iQ,iQ)=KQQ;

% Rhs global
RhsG(iV,1)=fV;
RhsG(iP,1)=fP;
RhsG(iB,1)=fB;
RhsG(iQ,1)=fQ;

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

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementNodesPostHigh=size(RefElement.PostHighHigh.NodesCoordElem,1);
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.DynamicViscosity;
lambda=-2/3*Parameters.DynamicViscosity;
eta=Parameters.MagneticDiffusivity;
Xe=Nodes';

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
ve=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
Je=reshape(SolutionLocal(:,msd+nsd+1+(1:qsd)),[],1);
be=reshape(SolutionLocal(:,msd+nsd+1+qsd+(1:nsd)),[],1);
Ue=SolutionGlobal;

% Initialize lhs
KVpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
KVtp=zeros(nsd,nsd*NumElementNodesPost);
KVrp=zeros(qsd,nsd*NumElementNodesPost);
KBpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
KBhp=zeros(NumElementNodesPostHigh,nsd*NumElementNodesPost);

% Initialize rhs
fVp=zeros(nsd*NumElementNodesPost,1);
fVt=zeros(nsd,1);
fVr=zeros(qsd,1);
fBp=zeros(nsd*NumElementNodesPost,1);
fBh=zeros(NumElementNodesPostHigh,1);

% Get reference data
FaceNodes=RefElement.FaceNodesElem;

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions('Element',RefElement.PostLow,RefElement.Post,Xe,nsd);
[~,Nhhex,Nhhey,Nhhez,wheg,Nlhe]=mapShapeFunctions('Element',RefElement.PostLowHigh,...
                                                            RefElement.PostHighHigh,Xe,nsd);
Nhe=RefElement.PostHigh.ShapeFunctionsElem;
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
nhe1=1:NumElementNodesPostHigh;

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
vxe=ve(nle1);
vye=ve(nle2);
if nsd==3
  vze=ve(nle3);
end
if nsd==2
  Jze=Je(nle1);
elseif nsd==3
  Jxe=Je(nle1);
  Jye=Je(nle2);
  Jze=Je(nle3);
end
bxe=be(nle1);
bye=be(nle2);
if nsd==3
  bze=be(nle3);
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
vxeg=Nle*vxe;
vyeg=Nle*vye;
if nsd==3
  vzeg=Nle*vze;
end
if nsd==2
  Jzeg=Nle*Jze;
elseif nsd==3
  Jxeg=Nle*Jxe;
  Jyeg=Nle*Jye;
  Jzeg=Nle*Jze;
end
bxeg=Nlhe*bxe;
byeg=Nlhe*bye;
if nsd==3
  bzeg=Nlhe*bze;
end

% Compute basic matrices
Nw1eT=(weg.*N1e)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
NwhhexT=(wheg.*Nhhex)';
NwhheyT=(wheg.*Nhhey)';
if nsd==3
  NwhhezT=(wheg.*Nhhez)';
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
  VoigtV1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  VoigtV2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  VoigtV3=sqrt(mu);
elseif nsd==3
  VoigtV1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  VoigtV2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  VoigtV3=sqrt(mu);
end

% Compute lhs
if nsd==2
  KVpp(ne1,ne1)=-VoigtV1*Kxxe-VoigtV3*Kyye;
  KVpp(ne1,ne2)=-VoigtV2*Kxye-VoigtV3*Kyxe;
  KVpp(ne2,ne1)=-VoigtV2*Kyxe-VoigtV3*Kxye;
  KVpp(ne2,ne2)=-VoigtV1*Kyye-VoigtV3*Kxxe;
elseif nsd==3
  KVpp(ne1,ne1)=-VoigtV1*Kxxe-VoigtV3*Kyye-VoigtV3*Kzze;
  KVpp(ne1,ne2)=-VoigtV2*Kxye-VoigtV3*Kyxe;
  KVpp(ne1,ne3)=-VoigtV2*Kxze-VoigtV3*Kzxe;
  KVpp(ne2,ne1)=-VoigtV2*Kyxe-VoigtV3*Kxye;
  KVpp(ne2,ne2)=-VoigtV1*Kyye-VoigtV3*Kxxe-VoigtV3*Kzze;
  KVpp(ne2,ne3)=-VoigtV2*Kyze-VoigtV3*Kzye;
  KVpp(ne3,ne1)=-VoigtV2*Kzxe-VoigtV3*Kxze;
  KVpp(ne3,ne2)=-VoigtV2*Kzye-VoigtV3*Kyze;
  KVpp(ne3,ne3)=-VoigtV1*Kzze-VoigtV3*Kxxe-VoigtV3*Kyye;
end

KVtp(1,ne1)=Nw1eT*(Ne);
KVtp(2,ne2)=Nw1eT*(Ne);
if nsd==3
  KVtp(3,ne3)=Nw1eT*(Ne);
end

if nsd==2
  KVrp(1,ne1)=-Nw1eT*(Ney);
  KVrp(1,ne2)=+Nw1eT*(Nex);
elseif nsd==3
  KVrp(1,ne2)=-Nw1eT*(Nez);
  KVrp(1,ne3)=+Nw1eT*(Ney);
  KVrp(2,ne1)=+Nw1eT*(Nez);
  KVrp(2,ne3)=-Nw1eT*(Nex);
  KVrp(3,ne1)=-Nw1eT*(Ney);
  KVrp(3,ne2)=+Nw1eT*(Nex);
end

if nsd==2
  KBpp(ne1,ne1)=+sqrt(eta)*Kyye;
  KBpp(ne1,ne2)=-sqrt(eta)*Kyxe;
  KBpp(ne2,ne1)=-sqrt(eta)*Kxye;
  KBpp(ne2,ne2)=+sqrt(eta)*Kxxe;
elseif nsd==3
  KBpp(ne1,ne1)=+sqrt(eta)*Kyye+sqrt(eta)*Kzze;
  KBpp(ne1,ne2)=-sqrt(eta)*Kyxe;
  KBpp(ne1,ne3)=-sqrt(eta)*Kzxe;
  KBpp(ne2,ne1)=-sqrt(eta)*Kxye;
  KBpp(ne2,ne2)=+sqrt(eta)*Kxxe+sqrt(eta)*Kzze;
  KBpp(ne2,ne3)=-sqrt(eta)*Kzye;
  KBpp(ne3,ne1)=-sqrt(eta)*Kxze;
  KBpp(ne3,ne2)=-sqrt(eta)*Kyze;
  KBpp(ne3,ne3)=+sqrt(eta)*Kxxe+sqrt(eta)*Kyye;
end

KBhp(nhe1,ne1)=NwhhexT*(Nhe);
KBhp(nhe1,ne2)=NwhheyT*(Nhe);
if nsd==3
  KBhp(nhe1,ne3)=NwhhezT*(Nhe);
end

% Compute rhs
if nsd==2
  fVp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg);
  fVp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg);
elseif nsd==3
  fVp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg)+NwezT*(Lxzeg);
  fVp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg)+NwezT*(Lyzeg);
  fVp(ne3,1)=NwexT*(Lxzeg)+NweyT*(Lyzeg)+NwezT*(Lzzeg);
end

fVt(1,1)=Nw1eT*(vxeg);
fVt(2,1)=Nw1eT*(vyeg);
if nsd==3
  fVt(3,1)=Nw1eT*(vzeg);
end

if nsd==2
  fBp(ne1,1)=-NweyT*(Jzeg);
  fBp(ne2,1)=+NwexT*(Jzeg);
elseif nsd==3
  fBp(ne1,1)=+NwezT*(Jyeg)...
             -NweyT*(Jzeg);
  fBp(ne2,1)=+NwexT*(Jzeg)...
             -NwezT*(Jxeg);
  fBp(ne3,1)=+NweyT*(Jxeg)...
             -NwexT*(Jyeg);
end

fBh(nhe1,1)=NwhhexT*(bxeg)...
           +NwhheyT*(byeg);
if nsd==3
  fBh(nhe1,1)=fBh(nhe1,1)+NwhhezT*(bzeg);
end

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions('Face',RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  
  % Indices
  nlefU1=(iFace-1)*(nsd+1+nsd-1+1)*NumFaceNodes+(1:NumFaceNodes);
  nlefU2=nlefU1+NumFaceNodes;
  nlefU3=nlefU2+NumFaceNodes;
  
  % Flip face
  Node2Match1stNode1=Faces.Interior(2,iFace);
  FlipFace=max(Node2Match1stNode1);
  if FlipFace
    order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
    nlefU1=nlefU1(order);
    nlefU2=nlefU2(order);
    nlefU3=nlefU3(order);
  end
  
  % Compute variables at nodes
  Vxf=Ue(nlefU1);
  Vyf=Ue(nlefU2);
  if nsd==3
    Vzf=Ue(nlefU3);
  end
  
  % Compute variables at Gauss points
  Vxfg=Nlf*Vxf;
  Vyfg=Nlf*Vyf;
  if nsd==3
    Vzfg=Nlf*Vzf;
  end
  
  % Compute basic matrices
  Nw1fT=(wfg.*N1f)';
  
  % Compute rhs
  if nsd==2
    fVr(1)=fVr(1)+Nw1fT*(-Vxfg.*ny+Vyfg.*nx);
  elseif nsd==3
    fVr(1)=fVr(1)+Nw1fT*(-Vyfg.*nz+Vzfg.*ny);
    fVr(2)=fVr(2)+Nw1fT*(+Vxfg.*nz-Vzfg.*nx);
    fVr(3)=fVr(3)+Nw1fT*(-Vxfg.*ny+Vyfg.*nx);
  end
end

% Indices
iVp=1:nsd*NumElementNodesPost;
iVt=iVp(end)+(1:nsd);
iVr=iVt(end)+(1:qsd);
iBp=iVr(end)+(1:nsd*NumElementNodesPost); iBp2=iVp(end)+(1:nsd*NumElementNodesPost);
iBh=iBp(end)+(1:NumElementNodesPostHigh);

% Initialization of lhs and rhs
LhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd+nsd*NumElementNodesPost+NumElementNodesPostHigh,...
              nsd*NumElementNodesPost+nsd*NumElementNodesPost);
RhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd+nsd*NumElementNodesPost+NumElementNodesPostHigh,1);

% Lhs for post-processing
LhsPost(iVp,iVp)=KVpp;
LhsPost(iVt,iVp)=KVtp;
LhsPost(iVr,iVp)=KVrp;
LhsPost(iBp,iBp2)=KBpp;
LhsPost(iBh,iBp2)=KBhp;

% Rhs for post-processing
RhsPost(iVp,1)=fVp;
RhsPost(iVt,1)=fVt;
RhsPost(iVr,1)=fVr;
RhsPost(iBp,1)=fBp;
RhsPost(iBh,1)=fBh;

end