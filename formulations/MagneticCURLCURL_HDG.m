classdef MagneticCURLCURL_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim-1+1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim+1;
    
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
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
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
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,~,RefElement,Sizes)
      RefElement.PostLowHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree  ,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+1,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHighHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+2,Parameters.Degree+2,Parameters.NodesDistribution);
      NodesElem=Elements.Nodes;
      SolutionLocalElem=Elements.SolutionLocal;
      LhsPost=cell(Sizes.NumElements,1);
      RhsPost=cell(Sizes.NumElements,1);
      parfor iElem=1:Sizes.NumElements
        [LhsPostElem,RhsPostElem]=...
          doPostProcessElement(NodesElem{iElem},...
          SolutionLocalElem{iElem},...
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
        Results(iD).ScaledMagneticCurl=[];
        Results(iD).MagneticInduction=[];
        Results(iD).LagrangeMultiplier=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledMagneticCurl(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2);
      Results(iD).MagneticInduction(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+(1:Sizes(iD).NumSpaceDim));
      Results(iD).LagrangeMultiplier(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim+1);
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).MagneticInductionPost=Block(iD,iD).SolutionPost;
      end
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isConvectiveTerm=strcmp(Parameters.ConvectiveTerm,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
eta=Parameters.MagneticDiffusivity;
bD=Parameters.MagneticInduction;
qD=Parameters.LagrangeMultiplier;
eNt=Parameters.TangentialElectricField;
bNn=Parameters.NormalMagneticInduction;
g=Parameters.Source;
if isConvectiveTerm
  v=Parameters.Velocity;
end
Xe=Nodes';
tauB=Parameters.StabMagneticInduction;
tauQ=Parameters.StabLagrangeMultiplier;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
Je=reshape(SolutionLocal(:,1:qsd),[],1);
be=reshape(SolutionLocal(:,qsd+(1:nsd)),[],1);
qe=reshape(SolutionLocal(:,qsd+nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  bolde=reshape(SolutionOld(:,qsd+(1:nsd),:),[],BDFo);
end

% Initialize lhs
KJJ=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KJb=zeros(qsd*NumElementNodes,nsd*NumElementNodes);
KJB=zeros(qsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KbJ=zeros(nsd*NumElementNodes,qsd*NumElementNodes);
Kbb=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kbq=zeros(nsd*NumElementNodes,NumElementNodes);
KbB=zeros(nsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KbQ=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Kqb=zeros(NumElementNodes,nsd*NumElementNodes);
Kqq=zeros(NumElementNodes,NumElementNodes);
KqQ=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KBJ=zeros((nsd-1)*NumElementFaces*NumFaceNodes,qsd*NumElementNodes);
KBb=zeros((nsd-1)*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KBB=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KBQ=zeros((nsd-1)*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KQb=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KQq=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KQQ=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fJ=zeros(qsd*NumElementNodes,1);
fb=zeros(nsd*NumElementNodes,1);
fq=zeros(NumElementNodes,1);
fB=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);
fQ=zeros(NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
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
  boldxe=bolde(ne1,:);
  boldye=bolde(ne2,:);
  if nsd==3
    boldze=bolde(ne3,:);
  end
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
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
  boldxeg=Ne*boldxe;
  boldyeg=Ne*boldye;
  if nsd==3
    boldzeg=Ne*boldze;
  end
end
geg=g(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
gxeg=geg(:,1);
gyeg=geg(:,2);
if nsd==3
  gzeg=geg(:,3);
end
if isConvectiveTerm
  veg=v(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
  vxeg=veg(:,1);
  vyeg=veg(:,2);
  if nsd==3
    vzeg=veg(:,3);
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

% Compute lhs
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

if isConvectiveTerm
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

if isConvectiveTerm
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf,nsd);
    
    % Check boundary
    isDirichlet=Faces.Dirichlet(iFace);
    isNatural=Faces.Natural(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nefU1=(iFace-1)*(nsd-1+1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
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
    Bt1f=Ue(nefU1);
    if nsd==3
      Bt2f=Ue(nefU2);
    end
    if nsd==2
      Qf=Ue(nefU2);
    elseif nsd==3
      Qf=Ue(nefU3);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
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
    if isDirichlet
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
    elseif isNatural
      eNtfg=eNt(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      eNtxfg=eNtfg(:,1);
      eNtyfg=eNtfg(:,2);
      if nsd==3
        eNtzfg=eNtfg(:,3);
      end
      bNnfg=bNn(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isConvectiveTerm
      vfg=v(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      vxfg=vfg(:,1);
      vyfg=vfg(:,2);
      if nsd==3
        vzfg=vfg(:,3);
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
    
    if isConvectiveTerm
      if nsd==2
        Kbb(nf1,nf1)=Kbb(nf1,nf1)+NwfT*((+vyfg.*ny).*Nf);
        Kbb(nf2,nf1)=Kbb(nf2,nf1)+NwfT*((-vyfg.*nx).*Nf);
        Kbb(nf1,nf2)=Kbb(nf1,nf2)+NwfT*((-vxfg.*ny).*Nf);
        Kbb(nf2,nf2)=Kbb(nf2,nf2)+NwfT*((+vxfg.*nx).*Nf);
      elseif nsd==3
        Kbb(nf1,nf1)=Kbb(nf1,nf1)+NwfT*((+vyfg.*ny+vzfg.*nz).*Nf);
        Kbb(nf2,nf1)=Kbb(nf2,nf1)+NwfT*((-vyfg.*nx).*Nf);
        Kbb(nf3,nf1)=Kbb(nf3,nf1)+NwfT*((-vzfg.*nx).*Nf);
        Kbb(nf1,nf2)=Kbb(nf1,nf2)+NwfT*((-vxfg.*ny).*Nf);
        Kbb(nf2,nf2)=Kbb(nf2,nf2)+NwfT*((+vxfg.*nx+vzfg.*nz).*Nf);
        Kbb(nf3,nf2)=Kbb(nf3,nf2)+NwfT*((-vzfg.*ny).*Nf);
        Kbb(nf1,nf3)=Kbb(nf1,nf3)+NwfT*((-vxfg.*nz).*Nf);
        Kbb(nf2,nf3)=Kbb(nf2,nf3)+NwfT*((-vyfg.*nz).*Nf);
        Kbb(nf3,nf3)=Kbb(nf3,nf3)+NwfT*((+vxfg.*nx+vyfg.*ny).*Nf);
      end
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
    
    if isConvectiveTerm && not(isDirichlet)
      if nsd==2
        KBb(nefB1,nf1)=KBb(nefB1,nf1)-NwfT*((+Rxt1.*vyfg.*ny...
                                             -Ryt1.*vyfg.*nx).*Nf);
        KBb(nefB1,nf2)=KBb(nefB1,nf2)-NwfT*((-Rxt1.*vxfg.*ny...
                                             +Ryt1.*vxfg.*nx).*Nf);
      elseif nsd==3
        KBb(nefB1,nf1)=KBb(nefB1,nf1)-NwfT*((+Rxt1.*(vyfg.*ny+vzfg.*nz)...
                                             -Ryt1.*(vyfg.*nx)...
                                             -Rzt1.*(vzfg.*nx)).*Nf);
        KBb(nefB2,nf1)=KBb(nefB2,nf1)-NwfT*((+Rxt2.*(vyfg.*ny+vzfg.*nz)...
                                             -Ryt2.*(vyfg.*nx)...
                                             -Rzt2.*(vzfg.*nx)).*Nf);
        KBb(nefB1,nf2)=KBb(nefB1,nf2)-NwfT*((-Rxt1.*(vxfg.*ny)...
                                             +Ryt1.*(vxfg.*nx+vzfg.*nz)...
                                             -Rzt1.*(vzfg.*ny)).*Nf);
        KBb(nefB2,nf2)=KBb(nefB2,nf2)-NwfT*((-Rxt2.*(vxfg.*ny)...
                                             +Ryt2.*(vxfg.*nx+vzfg.*nz)...
                                             -Rzt2.*(vzfg.*ny)).*Nf);
        KBb(nefB1,nf3)=KBb(nefB1,nf3)-NwfT*((-Rxt1.*(vxfg.*nz)...
                                             -Ryt1.*(vyfg.*nz)...
                                             +Rzt1.*(vxfg.*nx+vyfg.*ny)).*Nf);
        KBb(nefB2,nf3)=KBb(nefB2,nf3)-NwfT*((-Rxt2.*(vxfg.*nz)...
                                             -Ryt2.*(vyfg.*nz)...
                                             +Rzt2.*(vxfg.*nx+vyfg.*ny)).*Nf);
      end
    end
    
    if not(isDirichlet)
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
    
    if not(isDirichlet)
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
    
    if not(isDirichlet)
      KQb(nefQ1,nf1)=KQb(nefQ1,nf1)+Mfnx;
      KQb(nefQ1,nf2)=KQb(nefQ1,nf2)+Mfny;
      if nsd==3
        KQb(nefQ1,nf3)=KQb(nefQ1,nf3)+Mfnz;
      end
    end
    
    KQQ(nefQ1,nefQ1)=KQQ(nefQ1,nefQ1)-tauQ*Mf;
    
    % Compute rhs
    if nsd==2
      fJ(nf1,1)=fJ(nf1,1)-NwfT*(sqrt(eta)*(nx.*Btyfg-ny.*Btxfg));
    elseif nsd==3
      fJ(nf1,1)=fJ(nf1,1)-NwfT*(sqrt(eta)*(ny.*Btzfg-nz.*Btyfg));
      fJ(nf2,1)=fJ(nf2,1)-NwfT*(sqrt(eta)*(nz.*Btxfg-nx.*Btzfg));
      fJ(nf3,1)=fJ(nf3,1)-NwfT*(sqrt(eta)*(nx.*Btyfg-ny.*Btxfg));
    end
    
    if isConvectiveTerm
      if nsd==2
        fb(nf1,1)=fb(nf1,1)-NwfT*((bxfg.*vyfg-vxfg.*byfg).*ny);
        fb(nf2,1)=fb(nf2,1)-NwfT*((byfg.*vxfg-vyfg.*bxfg).*nx);
      elseif nsd==3
        fb(nf1,1)=fb(nf1,1)-NwfT*((bxfg.*vyfg-vxfg.*byfg).*ny...
                                 +(bxfg.*vzfg-vxfg.*bzfg).*nz);
        fb(nf2,1)=fb(nf2,1)-NwfT*((byfg.*vxfg-vyfg.*bxfg).*nx...
                                 +(byfg.*vzfg-vyfg.*bzfg).*nz);
        fb(nf3,1)=fb(nf3,1)-NwfT*((bzfg.*vxfg-vzfg.*bxfg).*nx...
                                 +(bzfg.*vyfg-vzfg.*byfg).*ny);
      end
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
    
    if isConvectiveTerm && not(isDirichlet)
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*(bxfg.*vyfg-vxfg.*byfg).*ny...
                                     +Ryt1.*(byfg.*vxfg-vyfg.*bxfg).*nx);
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(Rxt1.*((bxfg.*vyfg-vxfg.*byfg).*ny+...
                                             (bxfg.*vzfg-vxfg.*bzfg).*nz)...
                                     +Ryt1.*((byfg.*vxfg-vyfg.*bxfg).*nx+...
                                             (byfg.*vzfg-vyfg.*bzfg).*nz)...
                                     +Rzt1.*((bzfg.*vxfg-vzfg.*bxfg).*nx+...
                                             (bzfg.*vyfg-vzfg.*byfg).*ny));
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(Rxt2.*((bxfg.*vyfg-vxfg.*byfg).*ny+...
                                             (bxfg.*vzfg-vxfg.*bzfg).*nz)...
                                     +Ryt2.*((byfg.*vxfg-vyfg.*bxfg).*nx+...
                                             (byfg.*vzfg-vyfg.*bzfg).*nz)...
                                     +Rzt2.*((bzfg.*vxfg-vzfg.*bxfg).*nx+...
                                             (bzfg.*vyfg-vzfg.*byfg).*ny));
      end
    end
    
    if not(isDirichlet)
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
    
    if not(isDirichlet)
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
    
    if not(isDirichlet)
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
    
    if isDirichlet
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
    
    if isNatural
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
    
    if not(isDirichlet)
      if nsd==2
        fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(bxfg.*nx+byfg.*ny);
      elseif nsd==3
        fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(bxfg.*nx+byfg.*ny+bzfg.*nz);
      end
    end
    
    if not(isDirichlet)
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*qfg);
    end
    
    fQ(nefQ1,1)=fQ(nefQ1,1)+NwfT*(tauQ*Qfg);
    
    if isDirichlet
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*qDfg);
    end
    
    if isNatural
      fQ(nefQ1,1)=fQ(nefQ1,1)+NwfT*(bNnfg);
    end
    
    % Clear blocks
    if nsd==2
      KBJ(nefB1,nf1)=KJB(nf1,nefB1)'*not(isDirichlet);
    elseif nsd==3
      KBJ(nefB1,nf1)=KJB(nf1,nefB1)'*not(isDirichlet);
      KBJ(nefB1,nf2)=KJB(nf2,nefB1)'*not(isDirichlet);
      KBJ(nefB1,nf3)=KJB(nf3,nefB1)'*not(isDirichlet);
      KBJ(nefB2,nf1)=KJB(nf1,nefB2)'*not(isDirichlet);
      KBJ(nefB2,nf2)=KJB(nf2,nefB2)'*not(isDirichlet);
      KBJ(nefB2,nf3)=KJB(nf3,nefB2)'*not(isDirichlet);
    end
    
    KQq(nefQ1,nf1)=KqQ(nf1,nefQ1)'*not(isDirichlet);
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iJ=1:qsd*NumElementNodes;
ib=iJ(end)+(1:nsd*NumElementNodes);
iq=ib(end)+(1:NumElementNodes);
iB=reshape((0:NumElementFaces-1)*(nsd-1+1)*NumFaceNodes+repmat((1:(nsd-1)*NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iQ=reshape((0:NumElementFaces-1)*(nsd-1+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+(nsd-1)*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((qsd+nsd+1)*NumElementNodes,(qsd+nsd+1)*NumElementNodes);
LhsLG=zeros((qsd+nsd+1)*NumElementNodes,(nsd-1+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd-1+1)*NumElementFaces*NumFaceNodes,(qsd+nsd+1)*NumElementNodes);
LhsGG=zeros((nsd-1+1)*NumElementFaces*NumFaceNodes,(nsd-1+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((qsd+nsd+1)*NumElementNodes,1);
RhsG=zeros((nsd-1+1)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iJ,iJ)=KJJ;
LhsLL(iJ,ib)=KJb;
LhsLL(ib,iJ)=KbJ;
LhsLL(ib,ib)=Kbb;
LhsLL(ib,iq)=Kbq;
LhsLL(iq,ib)=Kqb;
LhsLL(iq,iq)=Kqq;

% Lhs local-global
LhsLG(iJ,iB)=KJB;
LhsLG(ib,iB)=KbB;
LhsLG(ib,iQ)=KbQ;
LhsLG(iq,iQ)=KqQ;

% Rhs local
RhsL(iJ,1)=fJ;
RhsL(ib,1)=fb;
RhsL(iq,1)=fq;

% Lhs global-local
LhsGL(iB,iJ)=KBJ;
LhsGL(iB,ib)=KBb;
LhsGL(iQ,ib)=KQb;
LhsGL(iQ,iq)=KQq;

% Lhs global-global
LhsGG(iB,iB)=KBB;
LhsGG(iB,iQ)=KBQ;
LhsGG(iQ,iQ)=KQQ;

% Rhs global
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
  Nodes,SolutionLocal,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementNodesPostHigh=size(RefElement.PostHighHigh.NodesCoordElem,1);
eta=Parameters.MagneticDiffusivity;
Xe=Nodes';

% Get solution
Je=reshape(SolutionLocal(:,1:qsd),[],1);
be=reshape(SolutionLocal(:,qsd+(1:nsd)),[],1);

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
  Kpp(ne1,ne1)=+sqrt(eta)*Kyye;
  Kpp(ne1,ne2)=-sqrt(eta)*Kyxe;
  Kpp(ne2,ne1)=-sqrt(eta)*Kxye;
  Kpp(ne2,ne2)=+sqrt(eta)*Kxxe;
elseif nsd==3
  Kpp(ne1,ne1)=+sqrt(eta)*Kyye+sqrt(eta)*Kzze;
  Kpp(ne1,ne2)=-sqrt(eta)*Kyxe;
  Kpp(ne1,ne3)=-sqrt(eta)*Kzxe;
  Kpp(ne2,ne1)=-sqrt(eta)*Kxye;
  Kpp(ne2,ne2)=+sqrt(eta)*Kxxe+sqrt(eta)*Kzze;
  Kpp(ne2,ne3)=-sqrt(eta)*Kzye;
  Kpp(ne3,ne1)=-sqrt(eta)*Kxze;
  Kpp(ne3,ne2)=-sqrt(eta)*Kyze;
  Kpp(ne3,ne3)=+sqrt(eta)*Kxxe+sqrt(eta)*Kyye;
end

Khp(nhe1,ne1)=NwhhexT*(Nhe);
Khp(nhe1,ne2)=NwhheyT*(Nhe);
if nsd==3
  Khp(nhe1,ne3)=NwhhezT*(Nhe);
end

% Compute rhs
if nsd==2
  fp(ne1,1)=-NweyT*(Jzeg);
  fp(ne2,1)=+NwexT*(Jzeg);
elseif nsd==3
  fp(ne1,1)=+NwezT*(Jyeg)...
            -NweyT*(Jzeg);
  fp(ne2,1)=+NwexT*(Jzeg)...
            -NwezT*(Jxeg);
  fp(ne3,1)=+NweyT*(Jxeg)...
            -NwexT*(Jyeg);
end

fh(nhe1,1)=NwhhexT*(bxeg)...
          +NwhheyT*(byeg);
if nsd==3
  fh(nhe1,1)=fh(nhe1,1)+NwhhezT*(bzeg);
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