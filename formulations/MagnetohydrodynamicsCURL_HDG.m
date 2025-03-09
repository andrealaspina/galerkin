classdef MagnetohydrodynamicsCURL_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim+1+...
                                 NumSpaceDim+1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2+NumSpaceDim+1+...
                               (NumSpaceDim*(NumSpaceDim+1)/2-NumSpaceDim)+...
                                NumSpaceDim*(NumSpaceDim+1)/2+NumSpaceDim+1;
    
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
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp-Sizes(iD).NumSpaceDim),...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
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
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp-Sizes(iD).NumSpaceDim),...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
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
        Results(iD).ScaledMagneticGradient=[];
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
        (1:(Sizes(iD).NumVoigtComp-Sizes(iD).NumSpaceDim)));
      Results(iD).ScaledMagneticGradient(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1+...
        (Sizes(iD).NumVoigtComp-Sizes(iD).NumSpaceDim)+(1:Sizes(iD).NumVoigtComp));
      Results(iD).MagneticInduction(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        Sizes(iD).NumSpaceDim+1+(Sizes(iD).NumVoigtComp-Sizes(iD).NumSpaceDim)+...
        Sizes(iD).NumVoigtComp+(1:Sizes(iD).NumSpaceDim));
      Results(iD).LagrangeMultiplier(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1+...
        (Sizes(iD).NumVoigtComp-Sizes(iD).NumSpaceDim)+Sizes(iD).NumVoigtComp+...
        Sizes(iD).NumSpaceDim+1);
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
sN=Parameters.PseudoTraction;
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
Ge=reshape(SolutionLocal(:,msd+nsd+1+qsd+(1:msd)),[],1);
be=reshape(SolutionLocal(:,msd+nsd+1+qsd+msd+(1:nsd)),[],1);
qe=reshape(SolutionLocal(:,msd+nsd+1+qsd+msd+nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  volde=reshape(SolutionOld(:,msd+(1:nsd),:),[],BDFo);
  polde=reshape(SolutionOld(:,msd+nsd+1,:),[],BDFo);
  bolde=reshape(SolutionOld(:,msd+nsd+1+qsd+msd+(1:nsd),:),[],BDFo);
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
KJB=zeros(qsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KGG=zeros(msd*NumElementNodes,msd*NumElementNodes);
KGb=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KGB=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kbv=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kbb=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kbq=zeros(nsd*NumElementNodes,NumElementNodes);
KbB=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KbV=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kqq=zeros(NumElementNodes,NumElementNodes);
KqB=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KqQ=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KVL=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KVv=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KVp=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementNodes);
KVV=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KPp=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KPP=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KBG=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KBb=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KBq=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementNodes);
KBB=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KQq=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KQQ=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fL=zeros(msd*NumElementNodes,1);
fv=zeros(nsd*NumElementNodes,1);
fp=zeros(NumElementNodes,1);
fJ=zeros(qsd*NumElementNodes,1);
fG=zeros(msd*NumElementNodes,1);
fb=zeros(nsd*NumElementNodes,1);
fq=zeros(NumElementNodes,1);
fV=zeros(nsd*NumElementFaces*NumFaceNodes,1);
fP=zeros(NumElementFaces*NumFaceNodes,1);
fB=zeros(nsd*NumElementFaces*NumFaceNodes,1);
fQ=zeros(NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

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
  Jxye=Je(ne1);
elseif nsd==3
  Jyze=Je(ne1);
  Jxze=Je(ne2);
  Jxye=Je(ne3);
end
if nsd==2
  Gxxe=Ge(ne1);
  Gyye=Ge(ne2);
  Gxye=Ge(ne3);
elseif nsd==3
  Gxxe=Ge(ne1);
  Gyye=Ge(ne2);
  Gzze=Ge(ne3);
  Gxye=Ge(ne4);
  Gxze=Ge(ne5);
  Gyze=Ge(ne6);
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
  Jxyeg=Ne*Jxye;
elseif nsd==3
  Jyzeg=Ne*Jyze;
  Jxzeg=Ne*Jxze;
  Jxyeg=Ne*Jxye;
end
if nsd==2
  Gxxeg=Ne*Gxxe;
  Gyyeg=Ne*Gyye;
  Gxyeg=Ne*Gxye;
elseif nsd==3
  Gxxeg=Ne*Gxxe;
  Gyyeg=Ne*Gyye;
  Gzzeg=Ne*Gzze;
  Gxyeg=Ne*Gxye;
  Gxzeg=Ne*Gxze;
  Gyzeg=Ne*Gyze;
end
bxeg=Ne*bxe;
byeg=Ne*bye;
if nsd==3
  bzeg=Ne*bze;
end
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

% Compute material matrix
VoigtB1=sqrt(2*eta);
VoigtB3=sqrt(eta);

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
    Kvp(ne3,ne1)=NweT*((drdpeg.*(1/dt*vzeg*alpha(1)+1/dt*voldzeg*alpha(2:BDFo+1,1))).*Ne);
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
    KvJ(ne1,ne1)=mu_m^(-1/2)*NweT*((-byeg).*Ne);
    KvJ(ne2,ne1)=mu_m^(-1/2)*NweT*((+bxeg).*Ne);
  elseif nsd==3
    KvJ(ne2,ne1)=mu_m^(-1/2)*NweT*((-bzeg).*Ne);
    KvJ(ne3,ne1)=mu_m^(-1/2)*NweT*((+byeg).*Ne);
    KvJ(ne1,ne2)=mu_m^(-1/2)*NweT*((+bzeg).*Ne);
    KvJ(ne3,ne2)=mu_m^(-1/2)*NweT*((-bxeg).*Ne);
    KvJ(ne1,ne3)=mu_m^(-1/2)*NweT*((-byeg).*Ne);
    KvJ(ne2,ne3)=mu_m^(-1/2)*NweT*((+bxeg).*Ne);
  end

  if nsd==2
    Kvb(ne2,ne1)=mu_m^(-1/2)*NweT*((+Jxyeg).*Ne);
    Kvb(ne1,ne2)=mu_m^(-1/2)*NweT*((-Jxyeg).*Ne);
  elseif nsd==3
    Kvb(ne2,ne1)=mu_m^(-1/2)*NweT*((+Jxyeg).*Ne);
    Kvb(ne3,ne1)=mu_m^(-1/2)*NweT*((-Jxyeg).*Ne);
    Kvb(ne1,ne2)=mu_m^(-1/2)*NweT*((-Jxyeg).*Ne);
    Kvb(ne3,ne2)=mu_m^(-1/2)*NweT*((+Jyzeg).*Ne);
    Kvb(ne1,ne3)=mu_m^(-1/2)*NweT*((+Jxyeg).*Ne);
    Kvb(ne2,ne3)=mu_m^(-1/2)*NweT*((-Jyzeg).*Ne);
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
  KJb(ne1,ne1)=-mu_m^(-1/2)*Cye;
  KJb(ne1,ne2)=+mu_m^(-1/2)*Cxe;
elseif nsd==3
  KJb(ne2,ne1)=+mu_m^(-1/2)*Cze;
  KJb(ne3,ne1)=-mu_m^(-1/2)*Cye;
  KJb(ne1,ne2)=-mu_m^(-1/2)*Cze;
  KJb(ne3,ne2)=+mu_m^(-1/2)*Cxe;
  KJb(ne1,ne3)=+mu_m^(-1/2)*Cye;
  KJb(ne2,ne3)=-mu_m^(-1/2)*Cxe;
end

KGG(ne1,ne1)=-Me;
KGG(ne2,ne2)=-Me;
KGG(ne3,ne3)=-Me;
if nsd==3
  KGG(ne4,ne4)=-Me;
  KGG(ne5,ne5)=-Me;
  KGG(ne6,ne6)=-Me;
end

if nsd==2
  KGb(ne1,ne1)=VoigtB1*Cxe;
  KGb(ne3,ne1)=VoigtB3*Cye;
  KGb(ne2,ne2)=VoigtB1*Cye;
  KGb(ne3,ne2)=VoigtB3*Cxe;
elseif nsd==3
  KGb(ne1,ne1)=VoigtB1*Cxe;
  KGb(ne4,ne1)=VoigtB3*Cye;
  KGb(ne5,ne1)=VoigtB3*Cze;
  KGb(ne2,ne2)=VoigtB1*Cye;
  KGb(ne4,ne2)=VoigtB3*Cxe;
  KGb(ne6,ne2)=VoigtB3*Cze;
  KGb(ne3,ne3)=VoigtB1*Cze;
  KGb(ne5,ne3)=VoigtB3*Cxe;
  KGb(ne6,ne3)=VoigtB3*Cye;
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

Kbq(ne1,ne1)=Kbq(ne1,ne1)+CxeT;
Kbq(ne2,ne1)=Kbq(ne2,ne1)+CyeT;
if nsd==3
  Kbq(ne3,ne1)=Kbq(ne3,ne1)+CzeT;
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
    fv(ne1,1)=fv(ne1,1)-NweT*(mu_m^(-1/2)*(-byeg.*Jxyeg));
    fv(ne2,1)=fv(ne2,1)-NweT*(mu_m^(-1/2)*(+bxeg.*Jxyeg));
  elseif nsd==3
    fv(ne1,1)=fv(ne1,1)-NweT*(mu_m^(-1/2)*(+bzeg.*Jxzeg...
                                           -byeg.*Jxyeg));
    fv(ne2,1)=fv(ne2,1)-NweT*(mu_m^(-1/2)*(-bzeg.*Jyzeg...
                                           +bxeg.*Jxyeg));
    fv(ne3,1)=fv(ne3,1)-NweT*(mu_m^(-1/2)*(+byeg.*Jyzeg...
                                           -bxeg.*Jxzeg));
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
  fJ(ne1,1)=+NweT*(Jxyeg)...
            -NwexT*(mu_m^(-1/2)*byeg)...
            +NweyT*(mu_m^(-1/2)*bxeg);
elseif nsd==3
  fJ(ne1,1)=+NweT*(Jyzeg)...
            -NweyT*(mu_m^(-1/2)*bzeg)...
            +NwezT*(mu_m^(-1/2)*byeg);
  fJ(ne2,1)=+NweT*(Jxzeg)...
            +NwexT*(mu_m^(-1/2)*bzeg)...
            -NwezT*(mu_m^(-1/2)*bxeg);
  fJ(ne3,1)=+NweT*(Jxyeg)...
            -NwexT*(mu_m^(-1/2)*byeg)...
            +NweyT*(mu_m^(-1/2)*bxeg);
end

if nsd==2
  fG(ne1,1)=+NweT*(Gxxeg)...
            -NwexT*(VoigtB1*bxeg);
  fG(ne2,1)=+NweT*(Gyyeg)...
            -NweyT*(VoigtB1*byeg);
  fG(ne3,1)=+NweT*(Gxyeg)...
            -NweyT*(VoigtB3*bxeg)...
            -NwexT*(VoigtB3*byeg);
elseif nsd==3
  fG(ne1,1)=+NweT*(Gxxeg)...
            -NwexT*(VoigtB1*bxeg);
  fG(ne2,1)=+NweT*(Gyyeg)...
            -NweyT*(VoigtB1*byeg);
  fG(ne3,1)=+NweT*(Gzzeg)...
            -NwezT*(VoigtB1*bzeg);
  fG(ne4,1)=+NweT*(Gxyeg)...
            -NwexT*(VoigtB3*byeg)...
            -NweyT*(VoigtB3*bxeg);
  fG(ne5,1)=+NweT*(Gxzeg)...
            -NwexT*(VoigtB3*bzeg)...
            -NwezT*(VoigtB3*bxeg);
  fG(ne6,1)=+NweT*(Gyzeg)...
            -NweyT*(VoigtB3*bzeg)...
            -NwezT*(VoigtB3*byeg);
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
  fb(ne1,1)=fb(ne1,1)-NweT*(+VoigtB1*(Nex*Gxxe)...
                            +VoigtB3*(Ney*Gxye)...
                            +(Nex*qe));
  fb(ne2,1)=fb(ne2,1)-NweT*(+VoigtB1*(Ney*Gyye)...
                            +VoigtB3*(Nex*Gxye)...
                            +(Ney*qe));
elseif nsd==3
  fb(ne1,1)=fb(ne1,1)-NweT*(+VoigtB1*(Nex*Gxxe)...
                            +VoigtB3*(Ney*Gxye)+VoigtB3*(Nez*Gxze)...
                            +(Nex*qe));
  fb(ne2,1)=fb(ne2,1)-NweT*(+VoigtB1*(Ney*Gyye)...
                            +VoigtB3*(Nex*Gxye)+VoigtB3*(Nez*Gyze)...
                            +(Ney*qe));
  fb(ne3,1)=fb(ne3,1)-NweT*(+VoigtB1*(Nez*Gzze)...
                            +VoigtB3*(Nex*Gxze)+VoigtB3*(Ney*Gyze)...
                            +(Nez*qe));
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf,nsd);
    
    % Check boundary
    isExterior=Faces.Exterior(iFace);
    isDirichlet_v_x=Faces.Dirichlet_v_x(iFace);
    isDirichlet_v_y=Faces.Dirichlet_v_y(iFace);
    if nsd==3; isDirichlet_v_z=Faces.Dirichlet_v_z(iFace); end
    isDirichlet_p=Faces.Dirichlet_p(iFace);
    isNeumann_t_x=Faces.Neumann_t_x(iFace);
    isNeumann_t_y=Faces.Neumann_t_y(iFace);
    if nsd==3; isNeumann_t_z=Faces.Neumann_t_z(iFace); end
    isDirichlet_b_x=Faces.Dirichlet_b_x(iFace);
    isDirichlet_b_y=Faces.Dirichlet_b_y(iFace);
    if nsd==3; isDirichlet_b_z=Faces.Dirichlet_b_z(iFace); end
    isDirichlet_q=Faces.Dirichlet_q(iFace);
    isNeumann_s_x=Faces.Neumann_s_x(iFace);
    isNeumann_s_y=Faces.Neumann_s_y(iFace);
    if nsd==3; isNeumann_s_z=Faces.Neumann_s_z(iFace); end
    isDirichlet_v=isDirichlet_v_x || isDirichlet_v_y || (nsd==3 && isDirichlet_v_z);
    isNeumann_t=isNeumann_t_x || isNeumann_t_y || (nsd==3 && isNeumann_t_z);
    isDirichlet_b=isDirichlet_b_x || isDirichlet_b_y || (nsd==3 && isDirichlet_b_z);
    isNeumann_s=isNeumann_s_x || isNeumann_s_y || (nsd==3 && isNeumann_s_z);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nf4=nf3+NumElementNodes;
    nf5=nf4+NumElementNodes;
    nf6=nf5+NumElementNodes;
    nefU1=(iFace-1)*(nsd+1+nsd+1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
    nefU4=nefU3+NumFaceNodes;
    nefU5=nefU4+NumFaceNodes;
    nefU6=nefU5+NumFaceNodes;
    nefU7=nefU6+NumFaceNodes;
    nefU8=nefU7+NumFaceNodes;
    nefV1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefV2=nefV1+NumFaceNodes;
    nefV3=nefV2+NumFaceNodes;
    nefP1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefB1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefB2=nefB1+NumFaceNodes;
    nefB3=nefB2+NumFaceNodes;
    nefQ1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    
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
      nefU8=nefU8(order);
      nefV1=nefV1(order);
      nefV2=nefV2(order);
      nefV3=nefV3(order);
      nefP1=nefP1(order);
      nefB1=nefB1(order);
      nefB2=nefB2(order);
      nefB3=nefB3(order);
      nefQ1=nefQ1(order);
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
      Gxxf=Gxxe(nf1);
      Gyyf=Gyye(nf1);
      Gxyf=Gxye(nf1);
    elseif nsd==3
      Gxxf=Gxxe(nf1);
      Gyyf=Gyye(nf1);
      Gzzf=Gzze(nf1);
      Gxyf=Gxye(nf1);
      Gxzf=Gxze(nf1);
      Gyzf=Gyze(nf1);
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
      Bxf=Ue(nefU4);
      Byf=Ue(nefU5);
    elseif nsd==3
      Bxf=Ue(nefU5);
      Byf=Ue(nefU6);
      Bzf=Ue(nefU7);
    end
    if nsd==2
      Qf=Ue(nefU6);
    elseif nsd==3
      Qf=Ue(nefU8);
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
      Gxxfg=Nf*Gxxf;
      Gyyfg=Nf*Gyyf;
      Gxyfg=Nf*Gxyf;
    elseif nsd==3
      Gxxfg=Nf*Gxxf;
      Gyyfg=Nf*Gyyf;
      Gzzfg=Nf*Gzzf;
      Gxyfg=Nf*Gxyf;
      Gxzfg=Nf*Gxzf;
      Gyzfg=Nf*Gyzf;
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
    Bxfg=Nf*Bxf;
    Byfg=Nf*Byf;
    if nsd==3
      Bzfg=Nf*Bzf;
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
    if isDirichlet_b
      bDfg=bD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      bDxfg=bDfg(:,1);
      bDyfg=bDfg(:,2);
      if nsd==3
        bDzfg=bDfg(:,3);
      end
    end
    if isDirichlet_q
      qDfg=qD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isNeumann_s
      sNfg=sN(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      sNxfg=sNfg(:,1);
      sNyfg=sNfg(:,2);
      if nsd==3
        sNzfg=sNfg(:,3);
      end
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
      KJB(nf1,nefB1)=KJB(nf1,nefB1)+mu_m^(-1/2)*Mfny;
      KJB(nf1,nefB2)=KJB(nf1,nefB2)-mu_m^(-1/2)*Mfnx;
    elseif nsd==3
      KJB(nf2,nefB1)=KJB(nf2,nefB1)-mu_m^(-1/2)*Mfnz;
      KJB(nf3,nefB1)=KJB(nf3,nefB1)+mu_m^(-1/2)*Mfny;
      KJB(nf1,nefB2)=KJB(nf1,nefB2)+mu_m^(-1/2)*Mfnz;
      KJB(nf3,nefB2)=KJB(nf3,nefB2)-mu_m^(-1/2)*Mfnx;
      KJB(nf1,nefB3)=KJB(nf1,nefB3)-mu_m^(-1/2)*Mfny;
      KJB(nf2,nefB3)=KJB(nf2,nefB3)+mu_m^(-1/2)*Mfnx;
    end
    
    if nsd==2
      KGB(nf1,nefB1)=KGB(nf1,nefB1)-VoigtB1*Mfnx;
      KGB(nf3,nefB1)=KGB(nf3,nefB1)-VoigtB3*Mfny;
      KGB(nf2,nefB2)=KGB(nf2,nefB2)-VoigtB1*Mfny;
      KGB(nf3,nefB2)=KGB(nf3,nefB2)-VoigtB3*Mfnx;
    elseif nsd==3
      KGB(nf1,nefB1)=KGB(nf1,nefB1)-VoigtB1*Mfnx;
      KGB(nf4,nefB1)=KGB(nf4,nefB1)-VoigtB3*Mfny;
      KGB(nf5,nefB1)=KGB(nf5,nefB1)-VoigtB3*Mfnz;
      KGB(nf2,nefB2)=KGB(nf2,nefB2)-VoigtB1*Mfny;
      KGB(nf4,nefB2)=KGB(nf4,nefB2)-VoigtB3*Mfnx;
      KGB(nf6,nefB2)=KGB(nf6,nefB2)-VoigtB3*Mfnz;
      KGB(nf3,nefB3)=KGB(nf3,nefB3)-VoigtB1*Mfnz;
      KGB(nf5,nefB3)=KGB(nf5,nefB3)-VoigtB3*Mfnx;
      KGB(nf6,nefB3)=KGB(nf6,nefB3)-VoigtB3*Mfny;
    end
    
    Kbb(nf1,nf1)=Kbb(nf1,nf1)+tauB*Mf;
    Kbb(nf2,nf2)=Kbb(nf2,nf2)+tauB*Mf;
    if nsd==3
      Kbb(nf3,nf3)=Kbb(nf3,nf3)+tauB*Mf;
    end
    
    if isCouplingTerms
      if nsd==2
        KbV(nf1,nefV1)=KbV(nf1,nefV1)+NwfT*((-Byfg.*ny).*Nf);
        KbV(nf2,nefV1)=KbV(nf2,nefV1)+NwfT*((+Byfg.*nx).*Nf);
        KbV(nf1,nefV2)=KbV(nf1,nefV2)+NwfT*((+Bxfg.*ny).*Nf);
        KbV(nf2,nefV2)=KbV(nf2,nefV2)+NwfT*((-Bxfg.*nx).*Nf);
      elseif nsd==3
        KbV(nf1,nefV1)=KbV(nf1,nefV1)+NwfT*((-Byfg.*ny-AdvMag_xx2).*Nf);
        KbV(nf2,nefV1)=KbV(nf2,nefV1)+NwfT*((+Byfg.*nx).*Nf);
        KbV(nf3,nefV1)=KbV(nf3,nefV1)+NwfT*((+Bzfg.*nx).*Nf);
        KbV(nf1,nefV2)=KbV(nf1,nefV2)+NwfT*((+Bxfg.*ny).*Nf);
        KbV(nf2,nefV2)=KbV(nf2,nefV2)+NwfT*((-Bxfg.*nx-Bzfg.*nz).*Nf);
        KbV(nf3,nefV2)=KbV(nf3,nefV2)+NwfT*((+Bzfg.*ny).*Nf);
        KbV(nf1,nefV3)=KbV(nf1,nefV3)+NwfT*((+Bxfg.*nz).*Nf);
        KbV(nf2,nefV3)=KbV(nf2,nefV3)+NwfT*((+Byfg.*nz).*Nf);
        KbV(nf3,nefV3)=KbV(nf3,nefV3)+NwfT*((-Bxfg.*nx-Byfg.*ny).*Nf);
      end
      
      if nsd==2
        KbB(nf1,nefB1)=KbB(nf1,nefB1)+NwfT*((+Vyfg.*ny).*Nf);
        KbB(nf2,nefB1)=KbB(nf2,nefB1)+NwfT*((-Vyfg.*nx).*Nf);
        KbB(nf1,nefB2)=KbB(nf1,nefB2)+NwfT*((-Vxfg.*ny).*Nf);
        KbB(nf2,nefB2)=KbB(nf2,nefB2)+NwfT*((+Vxfg.*nx).*Nf);
      elseif nsd==3
        KbB(nf1,nefB1)=KbB(nf1,nefB1)+NwfT*((+Vyfg.*ny+Vzfg.*nz).*Nf);
        KbB(nf2,nefB1)=KbB(nf2,nefB1)+NwfT*((-Vyfg.*nx).*Nf);
        KbB(nf3,nefB1)=KbB(nf3,nefB1)+NwfT*((-Vzfg.*nx).*Nf);
        KbB(nf1,nefB2)=KbB(nf1,nefB2)+NwfT*((-Vxfg.*ny).*Nf);
        KbB(nf2,nefB2)=KbB(nf2,nefB2)+NwfT*((+Vxfg.*nx+Vzfg.*nz).*Nf);
        KbB(nf3,nefB2)=KbB(nf3,nefB2)+NwfT*((-Vzfg.*ny).*Nf);
        KbB(nf1,nefB3)=KbB(nf1,nefB3)+NwfT*((-Vxfg.*nz).*Nf);
        KbB(nf2,nefB3)=KbB(nf2,nefB3)+NwfT*((-Vyfg.*nz).*Nf);
        KbB(nf3,nefB3)=KbB(nf3,nefB3)+NwfT*((+Vxfg.*nx+Vyfg.*ny).*Nf);
      end
    end
    
    KbB(nf1,nefB1)=KbB(nf1,nefB1)-tauB*Mf;
    KbB(nf2,nefB2)=KbB(nf2,nefB2)-tauB*Mf;
    if nsd==3
      KbB(nf3,nefB3)=KbB(nf3,nefB3)-tauB*Mf;
    end
    
    Kqq(nf1,nf1)=Kqq(nf1,nf1)-tauQ*Mf;
    
    KqB(nf1,nefB1)=KqB(nf1,nefB1)-Mfnx;
    KqB(nf1,nefB2)=KqB(nf1,nefB2)-Mfny;
    if nsd==3
      KqB(nf1,nefB3)=KqB(nf1,nefB3)-Mfnz;
    end
    
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
    
    if not(isDirichlet_b_x)
      KBb(nefB1,nf1)=KBb(nefB1,nf1)-tauB*Mf;
    end
    
    if not(isDirichlet_b_y)
      KBb(nefB2,nf2)=KBb(nefB2,nf2)-tauB*Mf;
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      KBb(nefB3,nf3)=KBb(nefB3,nf3)-tauB*Mf;
    end
    
    KBB(nefB1,nefB1)=KBB(nefB1,nefB1)+tauB*Mf;
    KBB(nefB2,nefB2)=KBB(nefB2,nefB2)+tauB*Mf;
    if nsd==3
      KBB(nefB3,nefB3)=KBB(nefB3,nefB3)+tauB*Mf;
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
    
    fp(nf1,1)=fp(nf1,1)-NwfT*((r(Pfg).*(Vxfg.*nx+Vyfg.*ny)-tauP*Pfg));
    if nsd==3
      fp(nf1,1)=fp(nf1,1)-NwfT*(r(Pfg).*Vzfg.*nz);
    end
    
    if nsd==2
      fJ(nf1,1)=fJ(nf1,1)+NwfT*(mu_m^(-1/2)*(nx.*Byfg-ny.*Bxfg));
    elseif nsd==3
      fJ(nf1,1)=fJ(nf1,1)+NwfT*(mu_m^(-1/2)*(ny.*Bzfg-nz.*Byfg));
      fJ(nf2,1)=fJ(nf2,1)+NwfT*(mu_m^(-1/2)*(nz.*Bxfg-nx.*Bzfg));
      fJ(nf3,1)=fJ(nf3,1)+NwfT*(mu_m^(-1/2)*(nx.*Byfg-ny.*Bxfg));
    end
    
    if nsd==2
      fG(nf1,1)=fG(nf1,1)+NwfT*(+VoigtB1*nx.*Bxfg);
      fG(nf2,1)=fG(nf2,1)+NwfT*(+VoigtB1*ny.*Byfg);
      fG(nf3,1)=fG(nf3,1)+NwfT*(+VoigtB3*ny.*Bxfg...
                                +VoigtB3*nx.*Byfg);
    elseif nsd==3
      fG(nf1,1)=fG(nf1,1)+NwfT*(+VoigtB1*nx.*Bxfg);
      fG(nf2,1)=fG(nf2,1)+NwfT*(+VoigtB1*ny.*Byfg);
      fG(nf3,1)=fG(nf3,1)+NwfT*(+VoigtB1*nz.*Bzfg);
      fG(nf4,1)=fG(nf4,1)+NwfT*(+VoigtB3*nx.*Byfg...
                                +VoigtB3*ny.*Bxfg);
      fG(nf5,1)=fG(nf5,1)+NwfT*(+VoigtB3*nx.*Bzfg...
                                +VoigtB3*nz.*Bxfg);
      fG(nf6,1)=fG(nf6,1)+NwfT*(+VoigtB3*ny.*Bzfg...
                                +VoigtB3*nz.*Byfg);
    end
    
    fb(nf1,1)=fb(nf1,1)-NwfT*(tauB*bxfg);
    fb(nf2,1)=fb(nf2,1)-NwfT*(tauB*byfg);
    if nsd==3
      fb(nf3,1)=fb(nf3,1)-NwfT*(tauB*bzfg);
    end
    
    if isCouplingTerms
      if nsd==2
        fb(nf1,1)=fb(nf1,1)-NwfT*((Bxfg.*Vyfg-Vxfg.*Byfg).*ny);
        fb(nf2,1)=fb(nf2,1)-NwfT*((Byfg.*Vxfg-Vyfg.*Bxfg).*nx);
      elseif nsd==3
        fb(nf1,1)=fb(nf1,1)-NwfT*((Bxfg.*Vyfg-Vxfg.*Byfg).*ny...
                                 +(Bxfg.*Vzfg-Vxfg.*Bzfg).*nz);
        fb(nf2,1)=fb(nf2,1)-NwfT*((Byfg.*Vxfg-Vyfg.*Bxfg).*nx...
                                  +(Byfg.*Vzfg-Vyfg.*Bzfg).*nz);
        fb(nf3,1)=fb(nf3,1)-NwfT*((Bzfg.*Vxfg-Vzfg.*Bxfg).*nx...
                                  +(Bzfg.*Vyfg-Vzfg.*Byfg).*ny);
      end
    end
    
    fb(nf1,1)=fb(nf1,1)+NwfT*(tauB*Bxfg);
    fb(nf2,1)=fb(nf2,1)+NwfT*(tauB*Byfg);
    if nsd==3
      fb(nf3,1)=fb(nf3,1)+NwfT*(tauB*Bzfg);
    end
    
    fq(nf1,1)=fq(nf1,1)+NwfT*(tauQ*qfg);
    
    fq(nf1,1)=fq(nf1,1)+NwfT*(Bxfg.*nx+Byfg.*ny-tauQ*Qfg);
    if nsd==3
      fq(nf1,1)=fq(nf1,1)+NwfT*(Bzfg.*nz);
    end
    
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
    
    if not(isDirichlet_p)
      fP(nefP1,1)=fP(nefP1,1)-NwfT*(tauP*pfg);
    end
    
    fP(nefP1,1)=fP(nefP1,1)+NwfT*(tauP*Pfg);
    
    if isDirichlet_p
      fP(nefP1,1)=fP(nefP1,1)-NwfT*(tauP*pDfg);
    end

    if not(isExterior) || isNeumann_s_x
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(+VoigtB1*nx.*Gxxfg...
                                      +VoigtB3*ny.*Gxyfg...
                                      +qfg.*nx);
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(+VoigtB1*nx.*Gxxfg...
                                      +VoigtB3*ny.*Gxyfg+VoigtB3*nz.*Gxzfg...
                                      +qfg.*nx);
      end
    end
    
    if not(isExterior) || isNeumann_s_y
      if nsd==2
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(+VoigtB1*ny.*Gyyfg...
                                      +VoigtB3*nx.*Gxyfg...
                                      +qfg.*ny);
      elseif nsd==3
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(+VoigtB1*ny.*Gyyfg...
                                      +VoigtB3*nx.*Gxyfg+VoigtB3*nz.*Gyzfg...
                                      +qfg.*ny);
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_s_z)
      fB(nefB3,1)=fB(nefB3,1)+NwfT*(+VoigtB1*nz.*Gzzfg...
                                    +VoigtB3*nx.*Gxzfg+VoigtB3*ny.*Gyzfg...
                                    +qfg.*nz);
    end
    
    if not(isDirichlet_b_x)
      fB(nefB1,1)=fB(nefB1,1)+NwfT*(+tauB*bxfg);
    end
    
    if not(isDirichlet_b_y)
      fB(nefB2,1)=fB(nefB2,1)+NwfT*(+tauB*byfg);
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      fB(nefB3,1)=fB(nefB3,1)+NwfT*(+tauB*bzfg);
    end
    
    if nsd==2
      fB(nefB1,1)=fB(nefB1,1)-NwfT*(tauB*Bxfg);
      fB(nefB2,1)=fB(nefB2,1)-NwfT*(tauB*Byfg);
    elseif nsd==3
      fB(nefB1,1)=fB(nefB1,1)-NwfT*(tauB*Bxfg);
      fB(nefB2,1)=fB(nefB2,1)-NwfT*(tauB*Byfg);
      fB(nefB3,1)=fB(nefB3,1)-NwfT*(tauB*Bzfg);
    end
    
    if isDirichlet_b_x
      fB(nefB1,1)=fB(nefB1,1)+NwfT*(tauB*bDxfg);
    end
    
    if isDirichlet_b_y
      fB(nefB2,1)=fB(nefB2,1)+NwfT*(tauB*bDyfg);
    end
    
    if nsd==3 && isDirichlet_b_z
      fB(nefB3,1)=fB(nefB3,1)+NwfT*(tauB*bDzfg);
    end
    
    if not(isDirichlet_q)
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*qfg);
    end
    
    fQ(nefQ1,1)=fQ(nefQ1,1)+NwfT*(tauQ*Qfg);
    
    if isDirichlet_q
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*qDfg);
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
    
    if isNeumann_s_x
      fB(nefB1,1)=fB(nefB1,1)+NwfT*(sNxfg);
    end
    
    if isNeumann_s_y
      fB(nefB2,1)=fB(nefB2,1)+NwfT*(sNyfg);
    end
    
    if nsd==3 && isNeumann_s_z
      fB(nefB3,1)=fB(nefB3,1)+NwfT*(sNzfg);
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
      KBG(nefB1,nf1)=KGB(nf1,nefB1)'*(not(isExterior) || isNeumann_s_x);
      KBG(nefB1,nf3)=KGB(nf3,nefB1)'*(not(isExterior) || isNeumann_s_x);
      KBG(nefB2,nf2)=KGB(nf2,nefB2)'*(not(isExterior) || isNeumann_s_y);
      KBG(nefB2,nf3)=KGB(nf3,nefB2)'*(not(isExterior) || isNeumann_s_y);
    elseif nsd==3
      KBG(nefB1,nf1)=KGB(nf1,nefB1)'*(not(isExterior) || isNeumann_s_x);
      KBG(nefB1,nf4)=KGB(nf4,nefB1)'*(not(isExterior) || isNeumann_s_x);
      KBG(nefB1,nf5)=KGB(nf5,nefB1)'*(not(isExterior) || isNeumann_s_x);
      KBG(nefB2,nf2)=KGB(nf2,nefB2)'*(not(isExterior) || isNeumann_s_y);
      KBG(nefB2,nf4)=KGB(nf4,nefB2)'*(not(isExterior) || isNeumann_s_y);
      KBG(nefB2,nf6)=KGB(nf6,nefB2)'*(not(isExterior) || isNeumann_s_y);
      KBG(nefB3,nf3)=KGB(nf3,nefB3)'*(not(isExterior) || isNeumann_s_z);
      KBG(nefB3,nf5)=KGB(nf5,nefB3)'*(not(isExterior) || isNeumann_s_z);
      KBG(nefB3,nf6)=KGB(nf6,nefB3)'*(not(isExterior) || isNeumann_s_z);
    end
    
    KBq(nefB1,nf1)=KqB(nf1,nefB1)'*(not(isExterior) || isNeumann_s_x);
    KBq(nefB2,nf1)=KqB(nf1,nefB2)'*(not(isExterior) || isNeumann_s_y);
    if nsd==3
      KBq(nefB3,nf1)=KqB(nf1,nefB3)'*(not(isExterior) || isNeumann_s_z);
    end
    
    KQq(nefQ1,nf1)=KqQ(nf1,nefQ1)'*not(isDirichlet_q);
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iL=1:msd*NumElementNodes;
iv=iL(end)+(1:nsd*NumElementNodes);
ip=iv(end)+(1:NumElementNodes);
iJ=ip(end)+(1:qsd*NumElementNodes);
iG=iJ(end)+(1:msd*NumElementNodes);
ib=iG(end)+(1:nsd*NumElementNodes);
iq=ib(end)+(1:NumElementNodes);
iV=reshape((0:NumElementFaces-1)*(nsd+1+nsd+1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iP=reshape((0:NumElementFaces-1)*(nsd+1+nsd+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+nsd*NumFaceNodes;
iB=reshape((0:NumElementFaces-1)*(nsd+1+nsd+1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[])+(nsd+1)*NumFaceNodes;
iQ=reshape((0:NumElementFaces-1)*(nsd+1+nsd+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+(nsd+1+nsd)*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+nsd+1+qsd+msd+nsd+1)*NumElementNodes,(msd+nsd+1+qsd+msd+nsd+1)*NumElementNodes);
LhsLG=zeros((msd+nsd+1+qsd+msd+nsd+1)*NumElementNodes,(nsd+1+nsd+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd+1+nsd+1)*NumElementFaces*NumFaceNodes,(msd+nsd+1+qsd+msd+nsd+1)*NumElementNodes);
LhsGG=zeros((nsd+1+nsd+1)*NumElementFaces*NumFaceNodes,(nsd+1+nsd+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+nsd+1+qsd+msd+nsd+1)*NumElementNodes,1);
RhsG=zeros((nsd+1+nsd+1)*NumElementFaces*NumFaceNodes,1);

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
LhsLL(iG,iG)=KGG;
LhsLL(iG,ib)=KGb;
LhsLL(ib,iv)=Kbv;
LhsLL(ib,iG)=KGb';
LhsLL(ib,ib)=Kbb;
LhsLL(ib,iq)=Kbq;
LhsLL(iq,ib)=Kbq';
LhsLL(iq,iq)=Kqq;

% Lhs local-global
LhsLG(iL,iV)=KLV;
LhsLG(iv,iV)=KvV;
LhsLG(iv,iP)=KvP;
LhsLG(ip,iV)=KpV;
LhsLG(ip,iP)=KpP;
LhsLG(iJ,iB)=KJB;
LhsLG(iG,iB)=KGB;
LhsLG(ib,iV)=KbV;
LhsLG(ib,iB)=KbB;
LhsLG(iq,iB)=KqB;
LhsLG(iq,iQ)=KqQ;

% Rhs local
RhsL(iL,1)=fL;
RhsL(iv,1)=fv;
RhsL(ip,1)=fp;
RhsL(iJ,1)=fJ;
RhsL(iG,1)=fG;
RhsL(ib,1)=fb;
RhsL(iq,1)=fq; 

% Lhs global-local
LhsGL(iV,iL)=KVL;
LhsGL(iV,iv)=KVv;
LhsGL(iV,ip)=KVp;
LhsGL(iP,ip)=KPp;
LhsGL(iB,iG)=KBG;
LhsGL(iB,ib)=KBb;
LhsGL(iB,iq)=KBq;
LhsGL(iQ,iq)=KQq;

% Lhs global-global
LhsGG(iV,iV)=KVV;
LhsGG(iP,iP)=KPP;
LhsGG(iB,iB)=KBB;
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
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.DynamicViscosity;
lambda=-2/3*Parameters.DynamicViscosity;
eta=Parameters.MagneticDiffusivity;
Xe=Nodes';

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
ve=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
Ge=reshape(SolutionLocal(:,msd+nsd+1+qsd+(1:msd)),[],1);
be=reshape(SolutionLocal(:,msd+nsd+1+qsd+msd+(1:nsd)),[],1);
Ue=SolutionGlobal;

% Initialize lhs
KVpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
KBpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
Ktp=zeros(nsd,nsd*NumElementNodesPost);
Krp=zeros(qsd,nsd*NumElementNodesPost);

% Initialize rhs
fVp=zeros(nsd*NumElementNodesPost,1);
fVt=zeros(nsd,1);
fVr=zeros(qsd,1);
fBp=zeros(nsd*NumElementNodesPost,1);
fBt=zeros(nsd,1);
fBr=zeros(qsd,1);

% Get reference data
FaceNodes=RefElement.FaceNodesElem;

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions(1,RefElement.PostLow,RefElement.Post,Xe,nsd);
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
vxe=ve(nle1);
vye=ve(nle2);
if nsd==3
  vze=ve(nle3);
end
if nsd==2
  Gxxe=Ge(nle1);
  Gyye=Ge(nle2);
  Gxye=Ge(nle3);
elseif nsd==3
  Gxxe=Ge(nle1);
  Gyye=Ge(nle2);
  Gzze=Ge(nle3);
  Gxye=Ge(nle4);
  Gxze=Ge(nle5);
  Gyze=Ge(nle6);
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
  Gxxeg=Nle*Gxxe;
  Gyyeg=Nle*Gyye;
  Gxyeg=Nle*Gxye;
elseif nsd==3
  Gxxeg=Nle*Gxxe;
  Gyyeg=Nle*Gyye;
  Gzzeg=Nle*Gzze;
  Gxyeg=Nle*Gxye;
  Gxzeg=Nle*Gxze;
  Gyzeg=Nle*Gyze;
end
bxeg=Nle*bxe;
byeg=Nle*bye;
if nsd==3
  bzeg=Nle*bze;
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
  VoigtV1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  VoigtV2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  VoigtV3=sqrt(mu);
elseif nsd==3
  VoigtV1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  VoigtV2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  VoigtV3=sqrt(mu);
end

% Compute material matrix
VoigtB1=sqrt(2*eta);
VoigtB3=sqrt(eta);

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

if nsd==2
  KBpp(ne1,ne1)=-VoigtB1*Kxxe-VoigtB3*Kyye;
  KBpp(ne1,ne2)=-VoigtB3*Kyxe;
  KBpp(ne2,ne1)=-VoigtB3*Kxye;
  KBpp(ne2,ne2)=-VoigtB1*Kyye-VoigtB3*Kxxe;
elseif nsd==3
  KBpp(ne1,ne1)=-VoigtB1*Kxxe-VoigtB3*Kyye-VoigtB3*Kzze;
  KBpp(ne1,ne2)=-VoigtB3*Kyxe;
  KBpp(ne1,ne3)=-VoigtB3*Kzxe;
  KBpp(ne2,ne1)=-VoigtB3*Kxye;
  KBpp(ne2,ne2)=-VoigtB1*Kyye-VoigtB3*Kxxe-VoigtB3*Kzze;
  KBpp(ne2,ne3)=-VoigtB3*Kzye;
  KBpp(ne3,ne1)=-VoigtB3*Kxze;
  KBpp(ne3,ne2)=-VoigtB3*Kyze;
  KBpp(ne3,ne3)=-VoigtB1*Kzze-VoigtB3*Kxxe-VoigtB3*Kyye;
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
  fBp(ne1,1)=NwexT*(Gxxeg)+NweyT*(Gxyeg);
  fBp(ne2,1)=NwexT*(Gxyeg)+NweyT*(Gyyeg);
elseif nsd==3
  fBp(ne1,1)=NwexT*(Gxxeg)+NweyT*(Gxyeg)+NwezT*(Gxzeg);
  fBp(ne2,1)=NwexT*(Gxyeg)+NweyT*(Gyyeg)+NwezT*(Gyzeg);
  fBp(ne3,1)=NwexT*(Gxzeg)+NweyT*(Gyzeg)+NwezT*(Gzzeg);
end

fBt(1,1)=Nw1eT*(bxeg);
fBt(2,1)=Nw1eT*(byeg);
if nsd==3
  fBt(3,1)=Nw1eT*(bzeg);
end

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions(0,RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  
  % Indices
  nlefU1=(iFace-1)*(nsd+1+nsd+1)*NumFaceNodes+(1:NumFaceNodes);
  nlefU2=nlefU1+NumFaceNodes;
  nlefU3=nlefU2+NumFaceNodes;
  nlefU4=nlefU3+NumFaceNodes;
  nlefU5=nlefU4+NumFaceNodes;
  nlefU6=nlefU5+NumFaceNodes;
  nlefU7=nlefU6+NumFaceNodes;
  
  % Flip face
  Node2Match1stNode1=Faces.Interior(2,iFace);
  FlipFace=max(Node2Match1stNode1);
  if FlipFace
    order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
    nlefU1=nlefU1(order);
    nlefU2=nlefU2(order);
    nlefU3=nlefU3(order);
    nlefU4=nlefU4(order);
    nlefU5=nlefU5(order);
    nlefU6=nlefU6(order);
    nlefU7=nlefU7(order);
  end
  
  % Compute variables at nodes
  Vxf=Ue(nlefU1);
  Vyf=Ue(nlefU2);
  if nsd==3
    Vzf=Ue(nlefU3);
  end
  if nsd==2
    Bxf=Ue(nlefU4);
    Byf=Ue(nlefU5);
  elseif nsd==3
    Bxf=Ue(nlefU5);
    Byf=Ue(nlefU6);
    Bzf=Ue(nlefU7);
  end
  
  % Compute variables at Gauss points
  Vxfg=Nlf*Vxf;
  Vyfg=Nlf*Vyf;
  if nsd==3
    Vzfg=Nlf*Vzf;
  end
  Bxfg=Nlf*Bxf;
  Byfg=Nlf*Byf;
  if nsd==3
    Bzfg=Nlf*Bzf;
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
  
  if nsd==2
    fBr(1)=fBr(1)+Nw1fT*(-Bxfg.*ny+Byfg.*nx);
  elseif nsd==3
    fBr(1)=fBr(1)+Nw1fT*(-Byfg.*nz+Bzfg.*ny);
    fBr(2)=fBr(2)+Nw1fT*(+Bxfg.*nz-Bzfg.*nx);
    fBr(3)=fBr(3)+Nw1fT*(-Bxfg.*ny+Byfg.*nx);
  end
end

% Indices
iVp=1:nsd*NumElementNodesPost;
iVt=iVp(end)+(1:nsd);
iVr=iVt(end)+(1:qsd);
iBp=iVr(end)+(1:nsd*NumElementNodesPost); iBp2=iVp(end)+(1:nsd*NumElementNodesPost);
iBt=iBp(end)+(1:nsd);
iBr=iBt(end)+(1:qsd);

% Initialization of lhs and rhs
LhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd+nsd*NumElementNodesPost+nsd+qsd,...
              nsd*NumElementNodesPost+nsd*NumElementNodesPost);
RhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd+nsd*NumElementNodesPost+nsd+qsd,1);

% Lhs for post-processing
LhsPost(iVp,iVp)=KVpp;
LhsPost(iVt,iVp)=Ktp;
LhsPost(iVr,iVp)=Krp;
LhsPost(iBp,iBp2)=KBpp;
LhsPost(iBt,iBp2)=Ktp;
LhsPost(iBr,iBp2)=Krp;

% Rhs for post-processing
RhsPost(iVp,1)=fVp;
RhsPost(iVt,1)=fVt;
RhsPost(iVr,1)=fVr;
RhsPost(iBp,1)=fBp;
RhsPost(iBt,1)=fBt;
RhsPost(iBr,1)=fBr;

end