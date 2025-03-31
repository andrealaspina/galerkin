classdef Electromagnetic_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim-1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) NumSpaceDim;

    % Discretization type
    DiscretizationType='HDG';

    % Time derivative order
    TimeDerOrder=1;
    
    % Domain
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
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,~)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          Parameters(iD).MagneticField(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).ElectricField(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              Parameters(iD).MagneticField(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).ElectricField(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
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
          Parameters,Time,RefElement,Sizes);
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
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,Time,RefElement,Sizes)
      RefElement.PostLowHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree  ,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+1,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHighHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+2,Parameters.Degree+2,Parameters.NodesDistribution);
      NodesElem=Elements.Nodes;
      SolutionLocalElem=Elements.SolutionLocal;
      SolutionOldElem=Elements.SolutionOld;
      LhsPost=cell(Sizes.NumElements,1);
      RhsPost=cell(Sizes.NumElements,1);
      parfor iElem=1:Sizes.NumElements
        [LhsPostElem,RhsPostElem]=...
          doPostProcessElement(NodesElem{iElem},...
          SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          Parameters,Time,RefElement,Sizes);
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
        Results(iD).MagneticField=[];
        Results(iD).ElectricField=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).MagneticField(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2);
      Results(iD).ElectricField(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+(1:Sizes(iD).NumSpaceDim));
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).ElectricFieldPost=Block(iD,iD).SolutionPost;
      end
    end
    
    %% Compute radiation
    function [Results]=computeRadiation(~,iD,Results,Elements,Mesh,Faces,RefElement,Sizes)
      Elements(iD).ResultsFFT=mat2cell(...
        [Results(iD).MagneticFieldFFT(Mesh(iD).Elements(:),:),...
         Results(iD).ElectricFieldFFT(Mesh(iD).Elements(:),:)],...
        ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp);
      NodesElem=Elements(iD).Nodes;
      ResultsFFTElem=Elements(iD).ResultsFFT;
      Radiation=0;
      for iFace=1:size(Faces(iD,iD).Radiation,1)
        iElem=Faces(iD,iD).Radiation(iFace,1);
        iElemFace=Faces(iD,iD).Radiation(iFace,2);
        [RadiationElem]=...
          computeRadiationElement(iElemFace,NodesElem{iElem},ResultsFFTElem{iElem},...
          RefElement,Sizes);
        Radiation=Radiation+RadiationElem;
      end
      Results(iD).Radiation=Radiation;
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.MagneticPermeability;
epsilon=Parameters.ElectricPermittivity;
sigma=Parameters.ElectricConductivity;
eD=Parameters.ElectricField;
j=Parameters.CurrentDensity;
g=Parameters.IncidentField;
Xe=Nodes';
tauE=Parameters.StabElectricField;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
He=reshape(SolutionLocal(:,1:qsd),[],1);
ee=reshape(SolutionLocal(:,qsd+(1:nsd)),[],1);
Ee=SolutionGlobal;
if isTimeDependent
  Holde=reshape(SolutionOld(:,1:qsd,:),[],BDFo);
  eolde=reshape(SolutionOld(:,qsd+(1:nsd),:),[],BDFo);
end

% Initialize lhs
KHH=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KHe=zeros(qsd*NumElementNodes,nsd*NumElementNodes);
KHE=zeros(qsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
Kee=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KeE=zeros(nsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KEE=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);

% Initialize rhs
fH=zeros(qsd*NumElementNodes,1);
fe=zeros(nsd*NumElementNodes,1);
fE=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Hze=He(ne1);
elseif nsd==3
  Hxe=He(ne1);
  Hye=He(ne2);
  Hze=He(ne3);
end
exe=ee(ne1);
eYe=ee(ne2);
if nsd==3
  eze=ee(ne3);
end
if isTimeDependent
  if nsd==2
    Holdze=Holde(ne1,:);
  elseif nsd==3
    Holdxe=Holde(ne1,:);
    Holdye=Holde(ne2,:);
    Holdze=Holde(ne3,:);
  end
  eoldxe=eolde(ne1,:);
  eoldye=eolde(ne2,:);
  if nsd==3
    eoldze=eolde(ne3,:);
  end
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
if nsd==2
  Hzeg=Ne*Hze;
elseif nsd==3
  Hxeg=Ne*Hxe;
  Hyeg=Ne*Hye;
  Hzeg=Ne*Hze;
end
exeg=Ne*exe;
eyeg=Ne*eYe;
if nsd==3
  ezeg=Ne*eze;
end
if isTimeDependent
  if nsd==2
    Holdzeg=Ne*Holdze;
  elseif nsd==3
    Holdxeg=Ne*Holdxe;
    Holdyeg=Ne*Holdye;
    Holdzeg=Ne*Holdze;
  end
  eoldxeg=Ne*eoldxe;
  eoldyeg=Ne*eoldye;
  if nsd==3
    eoldzeg=Ne*eoldze;
  end
end
jeg=j(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
jxeg=jeg(:,1);
jyeg=jeg(:,2);
if nsd==3
  jzeg=jeg(:,3);
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

% Compute lhs
if isTimeDependent
  KHH(ne1,ne1)=-mu*alpha(1)/dt*Me;
  if nsd==3
    KHH(ne2,ne2)=-mu*alpha(1)/dt*Me;
    KHH(ne3,ne3)=-mu*alpha(1)/dt*Me;
  end
end

if nsd==2
  KHe(ne1,ne1)=-Cye;
  KHe(ne1,ne2)=+Cxe;
elseif nsd==3
  KHe(ne2,ne1)=+Cze;
  KHe(ne3,ne1)=-Cye;
  KHe(ne1,ne2)=-Cze;
  KHe(ne3,ne2)=+Cxe;
  KHe(ne1,ne3)=+Cye;
  KHe(ne2,ne3)=-Cxe;
end

if isTimeDependent
  Kee(ne1,ne1)=epsilon*alpha(1)/dt*Me;
  Kee(ne2,ne2)=epsilon*alpha(1)/dt*Me;
  if nsd==3
    Kee(ne3,ne3)=epsilon*alpha(1)/dt*Me;
  end
end

Kee(ne1,ne1)=Kee(ne1,ne1)+sigma*Me;
Kee(ne2,ne2)=Kee(ne2,ne2)+sigma*Me;
if nsd==3
  Kee(ne3,ne3)=Kee(ne3,ne3)+sigma*Me;
end

% Compute rhs
if isTimeDependent
  if nsd==2
    fH(ne1,1)=+NweT*(mu/dt*Hzeg*alpha(1)...
                    +mu/dt*Holdzeg*alpha(2:BDFo+1,1));
  elseif nsd==3
    fH(ne1,1)=+NweT*(mu/dt*Hxeg*alpha(1)...
                    +mu/dt*Holdxeg*alpha(2:BDFo+1,1));
    fH(ne2,1)=+NweT*(mu/dt*Hyeg*alpha(1)...
                    +mu/dt*Holdyeg*alpha(2:BDFo+1,1));
    fH(ne3,1)=+NweT*(mu/dt*Hzeg*alpha(1)...
                    +mu/dt*Holdzeg*alpha(2:BDFo+1,1));
  end
end

if nsd==2
  fH(ne1,1)=fH(ne1,1)-NwexT*(eyeg)...
                     +NweyT*(exeg);
elseif nsd==3
  fH(ne1,1)=fH(ne1,1)-NweyT*(ezeg)...
                     +NwezT*(eyeg);
  fH(ne2,1)=fH(ne2,1)-NwezT*(exeg)...
                     +NwexT*(ezeg);
  fH(ne3,1)=fH(ne3,1)-NwexT*(eyeg)...
                     +NweyT*(exeg);
end

if isTimeDependent
  fe(ne1,1)=-NweT*(epsilon/dt*exeg*alpha(1)...
                  +epsilon/dt*eoldxeg*alpha(2:BDFo+1,1));
  fe(ne2,1)=-NweT*(epsilon/dt*eyeg*alpha(1)...
                  +epsilon/dt*eoldyeg*alpha(2:BDFo+1,1));
  if nsd==3
    fe(ne3,1)=-NweT*(epsilon/dt*ezeg*alpha(1)...
                    +epsilon/dt*eoldzeg*alpha(2:BDFo+1,1));
  end
end

fe(ne1,1)=fe(ne1,1)-NweT*(sigma*exeg);
fe(ne2,1)=fe(ne2,1)-NweT*(sigma*eyeg);
if nsd==3
  fe(ne3,1)=fe(ne3,1)-NweT*(sigma*ezeg);
end

if nsd==2
  fe(ne1,1)=fe(ne1,1)+NweT*(+Ney*Hze);
  fe(ne2,1)=fe(ne2,1)+NweT*(-Nex*Hze);
elseif nsd==3
  fe(ne1,1)=fe(ne1,1)+NweT*(-Nez*Hye+Ney*Hze);
  fe(ne2,1)=fe(ne2,1)+NweT*(-Nex*Hze+Nez*Hxe);
  fe(ne3,1)=fe(ne3,1)+NweT*(-Ney*Hxe+Nex*Hye);
end

fe(ne1,1)=fe(ne1,1)-NweT*(jxeg);
fe(ne2,1)=fe(ne2,1)-NweT*(jyeg);
if nsd==3
  fe(ne3,1)=fe(ne3,1)-NweT*(jzeg);
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
    isDirichlet=Faces.Dirichlet(iFace);
    isAbsorbing=Faces.Absorbing(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nef1=(iFace-1)*(nsd-1)*NumFaceNodes+(1:NumFaceNodes);
    nef2=nef1+NumFaceNodes;
    nefR=1:nsd;
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
      nef1=nef1(order);
      nef2=nef2(order);
      nefR=circshift(flip(nefR(1:nsd)),Node2Match1stNode1);
    end

    % Gram-Schmidt orthonormalization (straight faces)
    V=Xf(nefR,1:nsd)';
    W=(V(:,2:end)-V(:,1))./vecnorm(V(:,2:end)-V(:,1));
    R(:,1)=W(:,1);
    if nsd==3
      R(:,2)=(W(:,2)-(R(:,1)'*W(:,2))*R(:,1))/norm(W(:,2)-(R(:,1)'*W(:,2))*R(:,1));
    end
    
    % Get rotation operator
    if nsd==2
      Rxt1=R(1,1);
      Ryt1=R(2,1);
    elseif nsd==3
      Rxt1=R(1,1); Rxt2=R(1,2);
      Ryt1=R(2,1); Ryt2=R(2,2);
      Rzt1=R(3,1); Rzt2=R(3,2);
    end
    
    % Compute variables at nodes
    if nsd==2
      Hzf=Hze(nf1);
    elseif nsd==3
      Hxf=Hxe(nf1);
      Hyf=Hye(nf1);
      Hzf=Hze(nf1);
    end
    exf=exe(nf1);
    eyf=eYe(nf1);
    if nsd==3
      ezf=eze(nf1);
    end
    Et1f=Ee(nef1);
    if nsd==3
      Et2f=Ee(nef2);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    if nsd==2
      Hzfg=Nf*Hzf;
    elseif nsd==3
      Hxfg=Nf*Hxf;
      Hyfg=Nf*Hyf;
      Hzfg=Nf*Hzf;
    end
    exfg=Nf*exf;
    eyfg=Nf*eyf;
    if nsd==3
      ezfg=Nf*ezf;
    end
    Et1fg=Nf*Et1f;
    if nsd==3
      Et2fg=Nf*Et2f;
    end
    if nsd==2
      Etxfg=Rxt1*Et1fg;
      Etyfg=Ryt1*Et1fg;
    elseif nsd==3
      Etxfg=Rxt1*Et1fg+Rxt2*Et2fg;
      Etyfg=Ryt1*Et1fg+Ryt2*Et2fg;
      Etzfg=Rzt1*Et1fg+Rzt2*Et2fg;
    end
    if isDirichlet
      eDfg=eD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      if nsd==2
        eDt1fg=Rxt1*eDfg(:,1)+Ryt1*eDfg(:,2);
      elseif nsd==3
        eDt1fg=Rxt1*eDfg(:,1)+Ryt1*eDfg(:,2)+Rzt1*eDfg(:,3);
        eDt2fg=Rxt2*eDfg(:,1)+Ryt2*eDfg(:,2)+Rzt2*eDfg(:,3);
      end
      if nsd==2
        Etxfg=Rxt1*eDt1fg;
        Etyfg=Ryt1*eDt1fg;
      elseif nsd==3
        Etxfg=Rxt1*eDt1fg+Rxt2*eDt2fg;
        Etyfg=Ryt1*eDt1fg+Ryt2*eDt2fg;
        Etzfg=Rzt1*eDt1fg+Rzt2*eDt2fg;
      end
    elseif isAbsorbing
      gfg=g(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      gxfg=gfg(:,1);
      gyfg=gfg(:,2);
      if nsd==3
        gzfg=gfg(:,3);
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
    
    % Compute lhs
    if nsd==2
      Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*Mfnyny;
      Kee(nf2,nf1)=Kee(nf2,nf1)-tauE*Mfnxny;
      Kee(nf1,nf2)=Kee(nf1,nf2)-tauE*Mfnxny;
      Kee(nf2,nf2)=Kee(nf2,nf2)+tauE*Mfnxnx;
    elseif nsd==3
      Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*Mfnyny+tauE*Mfnznz;
      Kee(nf2,nf1)=Kee(nf2,nf1)-tauE*Mfnxny;
      Kee(nf3,nf1)=Kee(nf3,nf1)-tauE*Mfnxnz;
      Kee(nf1,nf2)=Kee(nf1,nf2)-tauE*Mfnxny;
      Kee(nf2,nf2)=Kee(nf2,nf2)+tauE*Mfnxnx+tauE*Mfnznz;
      Kee(nf3,nf2)=Kee(nf3,nf2)-tauE*Mfnynz;
      Kee(nf1,nf3)=Kee(nf1,nf3)-tauE*Mfnxnz;
      Kee(nf2,nf3)=Kee(nf2,nf3)-tauE*Mfnynz;
      Kee(nf3,nf3)=Kee(nf3,nf3)+tauE*Mfnxnx+tauE*Mfnyny;
    end
    
    if not(isDirichlet)                               
      if nsd==2
        KHE(nf1,nef1)=KHE(nf1,nef1)+Rxt1*Mfny...
                                   -Ryt1*Mfnx;
      elseif nsd==3
        KHE(nf1,nef1)=KHE(nf1,nef1)-Rzt1*Mfny...
                                   +Ryt1*Mfnz;
        KHE(nf2,nef1)=KHE(nf2,nef1)-Rxt1*Mfnz...
                                   +Rzt1*Mfnx;
        KHE(nf3,nef1)=KHE(nf3,nef1)-Ryt1*Mfnx...
                                   +Rxt1*Mfny;
        KHE(nf1,nef2)=KHE(nf1,nef2)-Rzt2*Mfny...
                                   +Ryt2*Mfnz;
        KHE(nf2,nef2)=KHE(nf2,nef2)-Rxt2*Mfnz...
                                   +Rzt2*Mfnx;
        KHE(nf3,nef2)=KHE(nf3,nef2)-Ryt2*Mfnx...
                                   +Rxt2*Mfny;
      end
      
      if nsd==2
        KeE(nf1,nef1)=KeE(nf1,nef1)-Rxt1*tauE*Mf;
        KeE(nf2,nef1)=KeE(nf2,nef1)-Ryt1*tauE*Mf;
      elseif nsd==3
        KeE(nf1,nef1)=KeE(nf1,nef1)-Rxt1*tauE*Mf;
        KeE(nf2,nef1)=KeE(nf2,nef1)-Ryt1*tauE*Mf;
        KeE(nf3,nef1)=KeE(nf3,nef1)-Rzt1*tauE*Mf;
        KeE(nf1,nef2)=KeE(nf1,nef2)-Rxt2*tauE*Mf;
        KeE(nf2,nef2)=KeE(nf2,nef2)-Ryt2*tauE*Mf;
        KeE(nf3,nef2)=KeE(nf3,nef2)-Rzt2*tauE*Mf;
      end
      
      if nsd==2
        KEE(nef1,nef1)=KEE(nef1,nef1)+Rxt1*Rxt1*tauE*Mf...
                                     +Ryt1*Ryt1*tauE*Mf;
      elseif nsd==3
        KEE(nef1,nef1)=KEE(nef1,nef1)+Rxt1*Rxt1*tauE*Mf...
                                     +Ryt1*Ryt1*tauE*Mf...
                                     +Rzt1*Rzt1*tauE*Mf;
        KEE(nef2,nef1)=KEE(nef2,nef1)+Rxt1*Rxt2*tauE*Mf...
                                     +Ryt1*Ryt2*tauE*Mf...
                                     +Rzt1*Rzt2*tauE*Mf;
        KEE(nef1,nef2)=KEE(nef1,nef2)+Rxt1*Rxt2*tauE*Mf...
                                     +Ryt1*Ryt2*tauE*Mf...
                                     +Rzt1*Rzt2*tauE*Mf;
        KEE(nef2,nef2)=KEE(nef2,nef2)+Rxt2*Rxt2*tauE*Mf...
                                     +Ryt2*Ryt2*tauE*Mf...
                                     +Rzt2*Rzt2*tauE*Mf;
      end
    end
    
    if isAbsorbing
      if nsd==2
        KEE(nef1,nef1)=KEE(nef1,nef1)+Rxt1*Rxt1*sqrt(epsilon/mu)*Mf...
                                     +Ryt1*Ryt1*sqrt(epsilon/mu)*Mf;
      elseif nsd==3
        KEE(nef1,nef1)=KEE(nef1,nef1)+Rxt1*Rxt1*sqrt(epsilon/mu)*Mf...
                                     +Ryt1*Ryt1*sqrt(epsilon/mu)*Mf...
                                     +Rzt1*Rzt1*sqrt(epsilon/mu)*Mf;
        KEE(nef2,nef1)=KEE(nef2,nef1)+Rxt1*Rxt2*sqrt(epsilon/mu)*Mf...
                                     +Ryt1*Ryt2*sqrt(epsilon/mu)*Mf...
                                     +Rzt1*Rzt2*sqrt(epsilon/mu)*Mf;
        KEE(nef1,nef2)=KEE(nef1,nef2)+Rxt1*Rxt2*sqrt(epsilon/mu)*Mf...
                                     +Ryt1*Ryt2*sqrt(epsilon/mu)*Mf...
                                     +Rzt1*Rzt2*sqrt(epsilon/mu)*Mf;
        KEE(nef2,nef2)=KEE(nef2,nef2)+Rxt2*Rxt2*sqrt(epsilon/mu)*Mf...
                                     +Ryt2*Ryt2*sqrt(epsilon/mu)*Mf...
                                     +Rzt2*Rzt2*sqrt(epsilon/mu)*Mf;
      end
    end
    
    % Compute rhs
    if nsd==2
      fH(nf1,1)=fH(nf1,1)+NwfT*(nx.*Etyfg-ny.*Etxfg);
    elseif nsd==3
      fH(nf1,1)=fH(nf1,1)+NwfT*(ny.*Etzfg-nz.*Etyfg);
      fH(nf2,1)=fH(nf2,1)+NwfT*(nz.*Etxfg-nx.*Etzfg);
      fH(nf3,1)=fH(nf3,1)+NwfT*(nx.*Etyfg-ny.*Etxfg);
    end
    
    if nsd==2
      fe(nf1,1)=fe(nf1,1)+NwfT*(tauE*(+ny.*(nx.*eyfg-ny.*exfg)));
      fe(nf2,1)=fe(nf2,1)+NwfT*(tauE*(-nx.*(nx.*eyfg-ny.*exfg)));
    elseif nsd==3
      fe(nf1,1)=fe(nf1,1)+NwfT*(tauE*(+ny.*(nx.*eyfg-ny.*exfg)...
                                      -nz.*(nz.*exfg-nx.*ezfg)));
      fe(nf2,1)=fe(nf2,1)+NwfT*(tauE*(+nz.*(ny.*ezfg-nz.*eyfg)...
                                      -nx.*(nx.*eyfg-ny.*exfg)));
      fe(nf3,1)=fe(nf3,1)+NwfT*(tauE*(+nx.*(nz.*exfg-nx.*ezfg)...
                                      -ny.*(ny.*ezfg-nz.*eyfg)));
    end
    
    fe(nf1,1)=fe(nf1,1)+NwfT*(tauE*Etxfg);
    fe(nf2,1)=fe(nf2,1)+NwfT*(tauE*Etyfg);
    if nsd==3
      fe(nf3,1)=fe(nf3,1)+NwfT*(tauE*Etzfg);
    end
    
    if not(isDirichlet)
      if nsd==2
        fE(nef1,1)=fE(nef1,1)+NwfT*(Rxt1*(-ny.*Hzfg)+...
                                    Ryt1*(+nx.*Hzfg));
      elseif nsd==3
        fE(nef1,1)=fE(nef1,1)+NwfT*(Rxt1*(nz.*Hyfg-ny.*Hzfg)+...
                                    Ryt1*(nx.*Hzfg-nz.*Hxfg)+...
                                    Rzt1*(ny.*Hxfg-nx.*Hyfg));
        fE(nef2,1)=fE(nef2,1)+NwfT*(Rxt2*(nz.*Hyfg-ny.*Hzfg)+...
                                    Ryt2*(nx.*Hzfg-nz.*Hxfg)+...
                                    Rzt2*(ny.*Hxfg-nx.*Hyfg));
      end
      
      if nsd==2
        fE(nef1,1)=fE(nef1,1)+NwfT*(Rxt1*(tauE*exfg)+...
                                    Ryt1*(tauE*eyfg));
      elseif nsd==3
        fE(nef1,1)=fE(nef1,1)+NwfT*(Rxt1*(tauE*exfg)+...
                                    Ryt1*(tauE*eyfg)+...
                                    Rzt1*(tauE*ezfg));
        fE(nef2,1)=fE(nef2,1)+NwfT*(Rxt2*(tauE*exfg)+...
                                    Ryt2*(tauE*eyfg)+...
                                    Rzt2*(tauE*ezfg));
      end
      
      if nsd==2
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*(tauE*Etxfg)+...
                                    Ryt1*(tauE*Etyfg));
      elseif nsd==3
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*(tauE*Etxfg)+...
                                    Ryt1*(tauE*Etyfg)+...
                                    Rzt1*(tauE*Etzfg));
        fE(nef2,1)=fE(nef2,1)-NwfT*(Rxt2*(tauE*Etxfg)+...
                                    Ryt2*(tauE*Etyfg)+...
                                    Rzt2*(tauE*Etzfg));
      end
    end
    
    if isAbsorbing
      if nsd==2
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*(sqrt(epsilon/mu)*Etxfg)+...
                                    Ryt1*(sqrt(epsilon/mu)*Etyfg));
      elseif nsd==3
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*(sqrt(epsilon/mu)*Etxfg)+...
                                    Ryt1*(sqrt(epsilon/mu)*Etyfg)+...
                                    Rzt1*(sqrt(epsilon/mu)*Etzfg));
        fE(nef2,1)=fE(nef2,1)-NwfT*(Rxt2*(sqrt(epsilon/mu)*Etxfg)+...
                                    Ryt2*(sqrt(epsilon/mu)*Etyfg)+...
                                    Rzt2*(sqrt(epsilon/mu)*Etzfg));
      end
      
      if nsd==2
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*(gxfg)+...
                                    Ryt1*(gyfg));
      elseif nsd==3
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*(gxfg)+...
                                    Ryt1*(gyfg)+...
                                    Rzt1*(gzfg));
        fE(nef2,1)=fE(nef2,1)-NwfT*(Rxt2*(gxfg)+...
                                    Ryt2*(gyfg)+...
                                    Rzt2*(gzfg));
      end
    end
    
    % Remove undetermination
    if isDirichlet
      KEE(nef1,nef1)=eye(NumFaceNodes);
      if nsd==3
        KEE(nef2,nef2)=eye(NumFaceNodes);
      end
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iH=1:qsd*NumElementNodes;
ie=iH(end)+(1:nsd*NumElementNodes);
iE=1:(nsd-1)*NumElementFaces*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((qsd+nsd)*NumElementNodes,(qsd+nsd)*NumElementNodes);
LhsLG=zeros((qsd+nsd)*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(qsd+nsd)*NumElementNodes);
LhsGG=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((qsd+nsd)*NumElementNodes,1);
RhsG=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iH,iH)=KHH;
LhsLL(iH,ie)=KHe;
LhsLL(ie,iH)=KHe';
LhsLL(ie,ie)=Kee;

% Lhs local-global
LhsLG(iH,iE)=KHE;
LhsLG(ie,iE)=KeE;

% Rhs local
RhsL(iH,1)=fH;
RhsL(ie,1)=fe;

% Lhs global-local
LhsGL(iE,iH)=KHE';
LhsGL(iE,ie)=KeE';

% Lhs global-global
LhsGG(iE,iE)=KEE;

% Rhs global
RhsG(iE,1)=fE;

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
  Nodes,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementNodesPostHigh=size(RefElement.PostHighHigh.NodesCoordElem,1);
mu=Parameters.MagneticPermeability;
Xe=Nodes';
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
He=reshape(SolutionLocal(:,1:qsd),[],1);
ee=reshape(SolutionLocal(:,qsd+(1:nsd)),[],1);
if isTimeDependent
  Holde=reshape(SolutionOld(:,1:qsd,:),[],BDFo);
end

% Initialize lhs
Kpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
Khp=zeros(NumElementNodesPostHigh,nsd*NumElementNodesPost);

% Initialize rhs
fp=zeros(nsd*NumElementNodesPost,1);
fh=zeros(NumElementNodesPostHigh,1);

% Get reference data
Nhe=RefElement.PostHigh.ShapeFunctionsElem;

% Compute weights at Gauss points
[~,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions(1,RefElement.PostLow,RefElement.Post,Xe,nsd);
[~,Nhhex,Nhhey,Nhhez,wheg,Nlhe]=mapShapeFunctions(1,RefElement.PostLowHigh,...
                                                    RefElement.PostHighHigh,Xe,nsd);

% Indices
ne1=1:NumElementNodesPost;
ne2=ne1+NumElementNodesPost;
ne3=ne2+NumElementNodesPost;
nle1=1:NumElementNodes;
nle2=nle1+NumElementNodes;
nle3=nle2+NumElementNodes;
nhe1=1:NumElementNodesPostHigh;

% Compute variables at nodes
if nsd==2
  Hze=He(nle1);
elseif nsd==3
  Hxe=He(nle1);
  Hye=He(nle2);
  Hze=He(nle3);
end
exe=ee(nle1);
eYe=ee(nle2);
if nsd==3
  eze=ee(nle3);
end
if isTimeDependent
  if nsd==2
    Holdze=Holde(nle1,:);
  elseif nsd==3
    Holdxe=Holde(nle1,:);
    Holdye=Holde(nle2,:);
    Holdze=Holde(nle3,:);
  end
end

% Compute variables at Gauss points
if nsd==2
  Hzeg=Nle*Hze;
elseif nsd==3
  Hxeg=Nle*Hxe;
  Hyeg=Nle*Hye;
  Hzeg=Nle*Hze;
end
exeg=Nlhe*exe;
eyeg=Nlhe*eYe;
if nsd==3
  ezeg=Nlhe*eze;
end
if isTimeDependent
  if nsd==2
    Holdzeg=Nle*Holdze;
  elseif nsd==3
    Holdxeg=Nle*Holdxe;
    Holdyeg=Nle*Holdye;
    Holdzeg=Nle*Holdze;
  end
end

% Compute basic matrices
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

% Compute lhs
if nsd==2
  Kpp(ne1,ne1)=-Kyye;
  Kpp(ne1,ne2)=+Kyxe;
  Kpp(ne2,ne1)=+Kxye;
  Kpp(ne2,ne2)=-Kxxe;
elseif nsd==3
  Kpp(ne1,ne1)=-Kyye-Kzze;
  Kpp(ne1,ne2)=+Kyxe;
  Kpp(ne1,ne3)=+Kzxe;
  Kpp(ne2,ne1)=+Kxye;
  Kpp(ne2,ne2)=-Kxxe-Kzze;
  Kpp(ne2,ne3)=+Kzye;
  Kpp(ne3,ne1)=+Kxze;
  Kpp(ne3,ne2)=+Kyze;
  Kpp(ne3,ne3)=-Kxxe-Kyye;
end

Khp(nhe1,ne1)=NwhhexT*(Nhe);
Khp(nhe1,ne2)=NwhheyT*(Nhe);
if nsd==3
  Khp(nhe1,ne3)=NwhhezT*(Nhe);
end

% Compute rhs
if nsd==2
  fp(ne1,1)=-NweyT*(mu/dt*Hzeg*alpha(1)+mu/dt*Holdzeg*alpha(2:BDFo+1,1));
  fp(ne2,1)=+NwexT*(mu/dt*Hzeg*alpha(1)+mu/dt*Holdzeg*alpha(2:BDFo+1,1));
elseif nsd==3
  fp(ne1,1)=+NwezT*(mu/dt*Hyeg*alpha(1)+mu/dt*Holdyeg*alpha(2:BDFo+1,1))...
            -NweyT*(mu/dt*Hzeg*alpha(1)+mu/dt*Holdzeg*alpha(2:BDFo+1,1));
  fp(ne2,1)=+NwexT*(mu/dt*Hzeg*alpha(1)+mu/dt*Holdzeg*alpha(2:BDFo+1,1))...
            -NwezT*(mu/dt*Hxeg*alpha(1)+mu/dt*Holdxeg*alpha(2:BDFo+1,1));
  fp(ne3,1)=+NweyT*(mu/dt*Hxeg*alpha(1)+mu/dt*Holdxeg*alpha(2:BDFo+1,1))...
            -NwexT*(mu/dt*Hyeg*alpha(1)+mu/dt*Holdyeg*alpha(2:BDFo+1,1));
end

fh(nhe1,1)=NwhhexT*(exeg)...
          +NwhheyT*(eyeg);
if nsd==3
  fh(nhe1,1)=fh(nhe1,1)+NwhhezT*(ezeg);
end

% Indices
ip=1:nsd*NumElementNodesPost;
ih=ip(end)+(1:NumElementNodesPostHigh);

% Initialization of lhs and rhs
LhsPost=zeros(nsd*NumElementNodesPost+NumElementNodesPostHigh,nsd*NumElementNodesPost);
RhsPost=zeros(nsd*NumElementNodesPost+NumElementNodesPostHigh,1);

% Lhs for post-processing
LhsPost(ip,ip)=Kpp;
LhsPost(ih,ip)=Khp;

% Rhs for post-processing
RhsPost(ip,1)=fp;
RhsPost(ih,1)=fh;

end

%% Compute radiation element
function [Radiation]=computeRadiationElement(...
  iFace,Nodes,ResultsFFT,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
Xe=Nodes';

% Get solution
Hwe=reshape(ResultsFFT(:,1:qsd),[],1);
ewe=reshape(ResultsFFT(:,qsd+(1:nsd)),[],1);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Hwze=Hwe(ne1);
elseif nsd==3
  Hwxe=Hwe(ne1);
  Hwye=Hwe(ne2);
  Hwze=Hwe(ne3);
end
ewxe=ewe(ne1);
ewye=ewe(ne2);
if nsd==3
  ewze=ewe(ne3);
end

% Compute weights at Gauss points
FaceNodes=RefElement.FaceNodesElem;
Xf=Xe(FaceNodes(iFace,:),:);
[Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf,nsd);

% Indices
nf1=FaceNodes(iFace,:);

% Compute variables at nodes
if nsd==2
  Hwzf=Hwze(nf1);
elseif nsd==3
  Hwxf=Hwxe(nf1);
  Hwyf=Hwye(nf1);
  Hwzf=Hwze(nf1);
end
ewxf=ewxe(nf1);
ewyf=ewye(nf1);
if nsd==3
  ewzf=ewze(nf1);
end

% Compute variables at Gauss points
if nsd==2
  Hwzfg=Nf*Hwzf;
elseif nsd==3
  Hwxfg=Nf*Hwxf;
  Hwyfg=Nf*Hwyf;
  Hwzfg=Nf*Hwzf;
end
ewxfg=Nf*ewxf;
ewyfg=Nf*ewyf;
if nsd==3
  ewzfg=Nf*ewzf;
end

% Poynting vector at Gauss points
if nsd==2
  Sxfg=real(+ewyfg.*conj(Hwzfg));
  Syfg=real(-ewxfg.*conj(Hwzfg));
elseif nsd==3
  Sxfg=real(ewyfg.*conj(Hwzfg)-ewzfg.*conj(Hwyfg));
  Syfg=real(ewzfg.*conj(Hwxfg)-ewxfg.*conj(Hwzfg));
  Szfg=real(ewxfg.*conj(Hwyfg)-ewyfg.*conj(Hwxfg));
end

% Compute radiation
if nsd==2
  Radiation=sum((Sxfg.*nx+Syfg.*ny).*wfg);
elseif nsd==3
  Radiation=sum((Sxfg.*nx+Syfg.*ny+Szfg.*nz).*wfg);
end

end