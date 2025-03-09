classdef ElectromagneticPML_HDG_fast < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim-1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim+...
                                NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim;
    
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
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,~)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          Parameters(iD).MagneticField(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).ElectricField(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).MagneticFieldAux(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).ElectricFieldAux(...
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
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).MagneticFieldAux(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).ElectricFieldAux(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD,Block,Elements,~,Parameters,~,~,Time,RefElement,Sizes)
      % Extract element data
      NodesElem=Elements(iD).Nodes;
      FacesElem=Elements(iD).Faces;
      
      % Global lhs and local matrices
      if Time.TimeStep<=Time.BDFOrderEff
        LhsCoef=zeros(Sizes(iD).NumElementLhsCoef,Sizes(iD).NumElements);
        for iElem=1:Sizes(iD).NumElements
          [LhsGlobalElem,MatLocalElem,LhsLLinvElem,LhsGLElem]=...
            buildBlockElementLhs(NodesElem{iElem},FacesElem(iElem),...
            Parameters,Time,RefElement.Value,Sizes);
          LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
          Elements(iD).MatLocal{iElem}=MatLocalElem;
          Elements(iD).LhsLLinv{iElem}=LhsLLinvElem;
          Elements(iD).LhsGL{iElem}=LhsGLElem;
        end
        Block(iD,iD).LhsGlobal=fsparse(Block(iD,iD).LhsRowIndices,...
                                       Block(iD,iD).LhsColIndices,LhsCoef(:));
      end
      
      % Extract element data
      SolutionOldElem=Elements(iD).SolutionOld;
      
      % Element rhs
      RhsL=cell(Sizes(iD).NumElements,1);
      RhsG=cell(Sizes(iD).NumElements,1);
      parfor iElem=1:Sizes(iD).NumElements
        [RhsLElem,RhsGElem]=...
          buildBlockElementRhs(NodesElem{iElem},FacesElem(iElem),...
          SolutionOldElem{iElem},...
          Parameters,Time,RefElement.Value,Sizes); %#ok
        RhsL{iElem}=RhsLElem;
        RhsG{iElem}=RhsGElem;
      end
      
      % Global rhs and local vectors
      RhsCoef=zeros(Sizes(iD).NumElementRhsCoef,Sizes(iD).NumElements);
      for iElem=1:Sizes(iD).NumElements
        Elements(iD).VecLocal{iElem}=Elements(iD).LhsLLinv{iElem}*RhsL{iElem};
        RhsGlobalElem=RhsG{iElem}-Elements(iD).LhsGL{iElem}*Elements(iD).VecLocal{iElem};
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
      end
      Block(iD,iD).RhsGlobal=fsparse(Block(iD,iD).RhsRowIndices,1,RhsCoef(:));
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
        Results(iD).MagneticFieldAux=[];
        Results(iD).ElectricFieldAux=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).MagneticField(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2);
      Results(iD).ElectricField(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+(1:Sizes(iD).NumSpaceDim));
      Results(iD).MagneticFieldAux(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim+...
        (1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2));
      Results(iD).ElectricFieldAux(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim+...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+(1:Sizes(iD).NumSpaceDim));
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).ElectricFieldPost=Block(iD,iD).SolutionPost;
      end
    end
    
    %% Compute radiation
    function [Results]=computeRadiation(~,iD,Results,Elements,Mesh,Faces,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      ResultsFFTElem=mat2cell(...
        [Results(iD).MagneticFieldFFT(Mesh(iD).Elements(:),:),...
         Results(iD).ElectricFieldFFT(Mesh(iD).Elements(:),:)],...
        ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes,Sizes(iD).NumLocalComp/2);
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
    
    %% Compute divergence
    function [Results]=computeDivergence(~,iD,Results,Elements,Parameters,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      SolutionLocalElem=Elements(iD).SolutionLocal;
      MagneticFieldDiv=zeros(Sizes(iD).NumElementNodes,Sizes(iD).NumElements);
      ElectricFieldDiv=zeros(Sizes(iD).NumElementNodes,Sizes(iD).NumElements);
      parfor iElem=1:Sizes(iD).NumElements
        [MagneticFieldDivElem,ElectricFielDivdElem]=...
          computeDivergenceElement(NodesElem{iElem},...
          SolutionLocalElem{iElem},...
          Parameters,RefElement,Sizes);
        MagneticFieldDiv(:,iElem)=reshape(MagneticFieldDivElem',[],1);
        ElectricFieldDiv(:,iElem)=reshape(ElectricFielDivdElem',[],1);
      end
      Results(iD).MagneticFieldDiv=MagneticFieldDiv(:);
      Results(iD).ElectricFieldDiv=ElectricFieldDiv(:);
    end
    
  end
  
end

%% Build block element (lhs)
function [LhsGlobal,MatLocal,LhsLLinv,LhsGL]=buildBlockElementLhs(...
  Nodes,Faces,Parameters,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.MagneticPermeability;
epsilon=Parameters.ElectricPermittivity;
sigma=Parameters.ElectricConductivity;
Sigma=Parameters.DampingFunction;
Xe=Nodes';
tauE=Parameters.StabElectricField;
dt=Time.TimeStepSize;
alpha=Time.BDF1stDerEff;

% Initialize lhs
KHH=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KHe=zeros(qsd*NumElementNodes,nsd*NumElementNodes);
KHM=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KHE=zeros(qsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
Kee=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KeJ=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KeE=zeros(nsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KMH=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KMM=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KJe=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KJJ=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KEE=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe(1:nsd+1,:),nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
Sigmae=Sigma(Xe(:,1),Xe(:,2),Xe(:,3));
sigmaxe=Sigmae(:,1);
sigmaye=Sigmae(:,2);
if nsd==3
  sigmaze=Sigmae(:,3);
end

% Compute variables at Gauss points
sigmaxeg=Ne*sigmaxe;
sigmayeg=Ne*sigmaye;
if nsd==3
  sigmazeg=Ne*sigmaze;
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
KHH(ne1,ne1)=-mu*alpha(1)/dt*Me;
if nsd==3
  KHH(ne2,ne2)=-mu*alpha(1)/dt*Me;
  KHH(ne3,ne3)=-mu*alpha(1)/dt*Me;
end

if nsd==2
  KHe(ne1,ne1)=-Cye;
  KHe(ne1,ne2)=+Cxe;
else
  KHe(ne2,ne1)=+Cze;
  KHe(ne3,ne1)=-Cye;
  KHe(ne1,ne2)=-Cze;
  KHe(ne3,ne2)=+Cxe;
  KHe(ne1,ne3)=+Cye;
  KHe(ne2,ne3)=-Cxe;
end

if nsd==2
  KHH(ne1,ne1)=KHH(ne1,ne1)-mu*NweT*(((sigmaxeg+sigmayeg)).*Ne);
else
  KHH(ne1,ne1)=KHH(ne1,ne1)-mu*NweT*(((sigmayeg+sigmazeg-sigmaxeg)).*Ne);
  KHH(ne2,ne2)=KHH(ne2,ne2)-mu*NweT*(((sigmazeg+sigmaxeg-sigmayeg)).*Ne);
  KHH(ne3,ne3)=KHH(ne3,ne3)-mu*NweT*(((sigmaxeg+sigmayeg-sigmazeg)).*Ne);
end

KHM(ne1,ne1)=-Me;
if nsd==3
  KHM(ne2,ne2)=-Me;
  KHM(ne3,ne3)=-Me;
end

Kee(ne1,ne1)=epsilon*alpha(1)/dt*Me;
Kee(ne2,ne2)=epsilon*alpha(1)/dt*Me;
if nsd==3
  Kee(ne3,ne3)=epsilon*alpha(1)/dt*Me;
end

Kee(ne1,ne1)=Kee(ne1,ne1)+sigma*Me;
Kee(ne2,ne2)=Kee(ne2,ne2)+sigma*Me;
if nsd==3
  Kee(ne3,ne3)=Kee(ne3,ne3)+sigma*Me;
end

if nsd==2
  Kee(ne1,ne1)=Kee(ne1,ne1)+epsilon*NweT*(((sigmayeg-sigmaxeg)).*Ne);
  Kee(ne2,ne2)=Kee(ne2,ne2)+epsilon*NweT*(((sigmaxeg-sigmayeg)).*Ne);
else
  Kee(ne1,ne1)=Kee(ne1,ne1)+epsilon*NweT*(((sigmayeg+sigmazeg-sigmaxeg)).*Ne);
  Kee(ne2,ne2)=Kee(ne2,ne2)+epsilon*NweT*(((sigmazeg+sigmaxeg-sigmayeg)).*Ne);
  Kee(ne3,ne3)=Kee(ne3,ne3)+epsilon*NweT*(((sigmaxeg+sigmayeg-sigmazeg)).*Ne);
end

KeJ(ne1,ne1)=Me;
KeJ(ne2,ne2)=Me;
if nsd==3
  KeJ(ne3,ne3)=Me;
end

KMM(ne1,ne1)=alpha(1)/dt*Me;
if nsd==3
  KMM(ne2,ne2)=alpha(1)/dt*Me;
  KMM(ne3,ne3)=alpha(1)/dt*Me;
end

if nsd==3
  KMM(ne1,ne1)=KMM(ne1,ne1)+NweT*((sigmaxeg).*Ne);
  KMM(ne2,ne2)=KMM(ne2,ne2)+NweT*((sigmayeg).*Ne);
  KMM(ne3,ne3)=KMM(ne3,ne3)+NweT*((sigmazeg).*Ne);
end

if nsd==2
  KMH(ne1,ne1)=-mu*NweT*(((-sigmaxeg).*(-sigmayeg)).*Ne);
else
  KMH(ne1,ne1)=-mu*NweT*(((sigmaxeg-sigmayeg).*(sigmaxeg-sigmazeg)).*Ne);
  KMH(ne2,ne2)=-mu*NweT*(((sigmayeg-sigmaxeg).*(sigmayeg-sigmazeg)).*Ne);
  KMH(ne3,ne3)=-mu*NweT*(((sigmazeg-sigmaxeg).*(sigmazeg-sigmayeg)).*Ne);
end

KJJ(ne1,ne1)=-alpha(1)/dt*Me;
KJJ(ne2,ne2)=-alpha(1)/dt*Me;
if nsd==3
  KJJ(ne3,ne3)=-alpha(1)/dt*Me;
end

KJJ(ne1,ne1)=KJJ(ne1,ne1)-NweT*((sigmaxeg).*Ne);
KJJ(ne2,ne2)=KJJ(ne2,ne2)-NweT*((sigmayeg).*Ne);
if nsd==3
  KJJ(ne3,ne3)=KJJ(ne3,ne3)-NweT*((sigmazeg).*Ne);
end

if nsd==2
  KJe(ne1,ne1)=+epsilon*NweT*(((sigmaxeg-sigmayeg).*(sigmaxeg)).*Ne);
  KJe(ne2,ne2)=+epsilon*NweT*(((sigmayeg-sigmaxeg).*(sigmayeg)).*Ne);  
else
  KJe(ne1,ne1)=+epsilon*NweT*(((sigmaxeg-sigmayeg).*(sigmaxeg-sigmazeg)).*Ne);
  KJe(ne2,ne2)=+epsilon*NweT*(((sigmayeg-sigmaxeg).*(sigmayeg-sigmazeg)).*Ne);
  KJe(ne3,ne3)=+epsilon*NweT*(((sigmazeg-sigmaxeg).*(sigmazeg-sigmayeg)).*Ne);
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf(1:nsd,:),nsd);
    
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
    
    % Flip face (periodic (slave))
    if matchField(Faces,'Periodic') && not(isempty(Faces.Periodic))
      Node2Match1stNode1=Faces.Periodic(2,iFace);
      FlipFace=max(Node2Match1stNode1);
      if FlipFace
        order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
        nef1=nef1(order);
        nef2=nef2(order);
        nefR=circshift(flip(nefR(1:nsd)),Node2Match1stNode1);
      end
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
    else
      Rxt1=R(1,1); Rxt2=R(1,2);
      Ryt1=R(2,1); Ryt2=R(2,2);
      Rzt1=R(3,1); Rzt2=R(3,2);
    end
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    Mf=NwfT*Nf;
    Mf=(Mf+Mf')/2;
    Mfnx=Mf*nx;
    Mfny=Mf*ny;
    if nsd==3
      Mfnz=Mf*nz;
    end
    if nsd==2
      Mfnxnx=Mf*nx*nx;
      Mfnxny=Mf*nx*ny;
      Mfnyny=Mf*ny*ny;
    else
      Mfnxnx=Mf*nx*nx;
      Mfnxny=Mf*nx*ny;
      Mfnxnz=Mf*nx*nz;
      Mfnyny=Mf*ny*ny;
      Mfnynz=Mf*ny*nz;
      Mfnznz=Mf*nz*nz;
    end
    
    % Compute lhs
    if nsd==2
      Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*Mfnyny;
      Kee(nf2,nf1)=Kee(nf2,nf1)-tauE*Mfnxny;
      Kee(nf1,nf2)=Kee(nf1,nf2)-tauE*Mfnxny;
      Kee(nf2,nf2)=Kee(nf2,nf2)+tauE*Mfnxnx;
    else
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
      else
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
      else
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
      else
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
      else
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
    
    % Remove undetermination
    if isDirichlet
      KEE(nef1,nef1)=eye(NumFaceNodes);
      if nsd==3
        KEE(nef2,nef2)=eye(NumFaceNodes);
      end
    end
  end
end

% Compute elemental contributions to lhs

% Indices
iH=1:qsd*NumElementNodes;
ie=iH(end)+(1:nsd*NumElementNodes);
iM=ie(end)+(1:qsd*NumElementNodes);
iJ=iM(end)+(1:nsd*NumElementNodes);
iE=1:(nsd-1)*NumElementFaces*NumFaceNodes;

% Initialization of lhs
LhsLL=zeros((qsd+nsd+qsd+nsd)*NumElementNodes,(qsd+nsd+qsd+nsd)*NumElementNodes);
LhsLG=zeros((qsd+nsd+qsd+nsd)*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(qsd+nsd+qsd+nsd)*NumElementNodes);
LhsGG=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);

% Lhs local-local
LhsLL(iH,iH)=KHH;
LhsLL(iH,ie)=KHe;
LhsLL(iH,iM)=KHM;
LhsLL(ie,iH)=KHe';
LhsLL(ie,ie)=Kee;
LhsLL(ie,iJ)=KeJ;
LhsLL(iM,iH)=KMH;
LhsLL(iM,iM)=KMM;
LhsLL(iJ,ie)=KJe;
LhsLL(iJ,iJ)=KJJ;

% Add missing term for axisymmetric formulation
if nsd==2 && matchField(Parameters,'Axisymmetric','yes')
  KeH_axisym=zeros(nsd*NumElementNodes,qsd*NumElementNodes);
  KeH_axisym(ne1,ne1)=-NweT*((1./(Ne*Xe(:,2))).*Ne);
  LhsLL(ie,iH)=LhsLL(ie,iH)+KeH_axisym;
end

% Lhs local-global
LhsLG(iH,iE)=KHE;
LhsLG(ie,iE)=KeE;

% Lhs global-local
LhsGL(iE,iH)=KHE';
LhsGL(iE,ie)=KeE';

% Lhs global-global
LhsGG(iE,iE)=KEE;

% Inverse of lhs local-local
LhsLLinv=pinv(LhsLL);

% Matrix for local problem
MatLocal=LhsLLinv*LhsLG;

% Lhs for global problem
LhsGlobal=LhsGG-LhsGL*MatLocal;

end

%% Build block element (rhs)
function [RhsL,RhsG]=buildBlockElementRhs(...
  Nodes,Faces,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.MagneticPermeability;
epsilon=Parameters.ElectricPermittivity;
eD=Parameters.ElectricField;
j=Parameters.CurrentDensity;
g=Parameters.IncidentField;
Xe=Nodes';
tauE=Parameters.StabElectricField;
t=Time.Time;
dt=Time.TimeStepSize;
BDFo=Time.BDFOrderEff;
alpha=Time.BDF1stDerEff(2:BDFo+1,1);

% Get solution
if nsd==2
  Holde=reshape(SolutionOld(:,1,:),[],BDFo);
  eolde=reshape(SolutionOld(:,[2,3],:),[],BDFo);
  Molde=reshape(SolutionOld(:,4,:),[],BDFo);
  Jolde=reshape(SolutionOld(:,[5,6],:),[],BDFo);
else
  Holde=reshape(SolutionOld(:,[1,2,3],:),[],BDFo);
  eolde=reshape(SolutionOld(:,[4,5,6],:),[],BDFo);
  Molde=reshape(SolutionOld(:,[7,8,9],:),[],BDFo);
  Jolde=reshape(SolutionOld(:,[10,11,12],:),[],BDFo);
end

% Initialize rhs
fH=zeros(qsd*NumElementNodes,1);
fe=zeros(nsd*NumElementNodes,1);
fM=zeros(qsd*NumElementNodes,1);
fJ=zeros(nsd*NumElementNodes,1);
fE=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,~,~,~,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe(1:nsd+1,:),nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Holdze=Holde(ne1,:);
else
  Holdxe=Holde(ne1,:);
  Holdye=Holde(ne2,:);
  Holdze=Holde(ne3,:);
end
eoldxe=eolde(ne1,:);
eoldye=eolde(ne2,:);
if nsd==3
  eoldze=eolde(ne3,:);
end
if nsd==2
  Moldze=Molde(ne1,:);
else
  Moldxe=Molde(ne1,:);
  Moldye=Molde(ne2,:);
  Moldze=Molde(ne3,:);
end
Joldxe=Jolde(ne1,:);
Joldye=Jolde(ne2,:);
if nsd==3
  Joldze=Jolde(ne3,:);
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
if nsd==2
  Holdzeg=Ne*Holdze;
else
  Holdxeg=Ne*Holdxe;
  Holdyeg=Ne*Holdye;
  Holdzeg=Ne*Holdze;
end
eoldxeg=Ne*eoldxe;
eoldyeg=Ne*eoldye;
if nsd==3
  eoldzeg=Ne*eoldze;
end
if nsd==2
  Moldzeg=Ne*Moldze;
else
  Moldxeg=Ne*Moldxe;
  Moldyeg=Ne*Moldye;
  Moldzeg=Ne*Moldze;
end
Joldxeg=Ne*Joldxe;
Joldyeg=Ne*Joldye;
if nsd==3
  Joldzeg=Ne*Joldze;
end
jeg=j(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
jxeg=jeg(:,1);
jyeg=jeg(:,2);
if nsd==3
  jzeg=jeg(:,3);
end

% Compute basic matrices
NweT=(weg.*Ne)';

% Compute rhs
if nsd==2
  fH(ne1,1)=+mu/dt*NweT*(Holdzeg*alpha);
else
  fH(ne1,1)=+mu/dt*NweT*(Holdxeg*alpha);
  fH(ne2,1)=+mu/dt*NweT*(Holdyeg*alpha);
  fH(ne3,1)=+mu/dt*NweT*(Holdzeg*alpha);
end

fe(ne1,1)=-epsilon/dt*NweT*(eoldxeg*alpha);
fe(ne2,1)=-epsilon/dt*NweT*(eoldyeg*alpha);
if nsd==3
  fe(ne3,1)=-epsilon/dt*NweT*(eoldzeg*alpha);
end

fe(ne1,1)=fe(ne1,1)-NweT*(jxeg);
fe(ne2,1)=fe(ne2,1)-NweT*(jyeg);
if nsd==3
  fe(ne3,1)=fe(ne3,1)-NweT*(jzeg);
end

if nsd==2
  fM(ne1,1)=-1/dt*NweT*(Moldzeg*alpha);
else
  fM(ne1,1)=-1/dt*NweT*(Moldxeg*alpha);
  fM(ne2,1)=-1/dt*NweT*(Moldyeg*alpha);
  fM(ne3,1)=-1/dt*NweT*(Moldzeg*alpha);
end

fJ(ne1,1)=+1/dt*NweT*(Joldxeg*alpha);
fJ(ne2,1)=+1/dt*NweT*(Joldyeg*alpha);
if nsd==3
  fJ(ne3,1)=+1/dt*NweT*(Joldzeg*alpha);
end

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=Faces.Exterior(iFace);
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement.FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf(1:nsd,:),nsd);
    
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
    
    % Flip face (periodic (slave))
    if matchField(Faces,'Periodic') && not(isempty(Faces.Periodic))
      Node2Match1stNode1=Faces.Periodic(2,iFace);
      FlipFace=max(Node2Match1stNode1);
      if FlipFace
        order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
        nef1=nef1(order);
        nef2=nef2(order);
        nefR=circshift(flip(nefR(1:nsd)),Node2Match1stNode1);
      end
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
    else
      Rxt1=R(1,1); Rxt2=R(1,2);
      Ryt1=R(2,1); Ryt2=R(2,2);
      Rzt1=R(3,1); Rzt2=R(3,2);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    if isDirichlet
      eDfg=eD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      if nsd==2
        eDt1fg=Rxt1*eDfg(:,1)+Ryt1*eDfg(:,2);
      else
        eDt1fg=Rxt1*eDfg(:,1)+Ryt1*eDfg(:,2)+Rzt1*eDfg(:,3);
        eDt2fg=Rxt2*eDfg(:,1)+Ryt2*eDfg(:,2)+Rzt2*eDfg(:,3);
      end
      if nsd==2
        eDtxfg=Rxt1*eDt1fg;
        eDtyfg=Ryt1*eDt1fg;
      else
        eDtxfg=Rxt1*eDt1fg+Rxt2*eDt2fg;
        eDtyfg=Ryt1*eDt1fg+Ryt2*eDt2fg;
        eDtzfg=Rzt1*eDt1fg+Rzt2*eDt2fg;
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
    
    % Compute rhs
    if isDirichlet      
      if nsd==2
        fH(nf1,1)=fH(nf1,1)+NwfT*(nx.*eDtyfg-ny.*eDtxfg);
      else
        fH(nf1,1)=fH(nf1,1)+NwfT*(ny.*eDtzfg-nz.*eDtyfg);
        fH(nf2,1)=fH(nf2,1)+NwfT*(nz.*eDtxfg-nx.*eDtzfg);
        fH(nf3,1)=fH(nf3,1)+NwfT*(nx.*eDtyfg-ny.*eDtxfg);
      end
      
      fe(nf1,1)=fe(nf1,1)+tauE*NwfT*(eDtxfg);
      fe(nf2,1)=fe(nf2,1)+tauE*NwfT*(eDtyfg);
      if nsd==3
        fe(nf3,1)=fe(nf3,1)+tauE*NwfT*(eDtzfg);
      end
    elseif isAbsorbing      
      if nsd==2
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*gxfg+...
                                    Ryt1*gyfg);
      else
        fE(nef1,1)=fE(nef1,1)-NwfT*(Rxt1*gxfg+...
                                    Ryt1*gyfg+...
                                    Rzt1*gzfg);
        fE(nef2,1)=fE(nef2,1)-NwfT*(Rxt2*gxfg+...
                                    Ryt2*gyfg+...
                                    Rzt2*gzfg);
      end
    end
  end
end

% Compute elemental contributions to rhs

% Indices
iH=1:qsd*NumElementNodes;
ie=iH(end)+(1:nsd*NumElementNodes);
iM=ie(end)+(1:qsd*NumElementNodes);
iJ=iM(end)+(1:nsd*NumElementNodes);
iE=1:(nsd-1)*NumElementFaces*NumFaceNodes;

% Initialization of rhs
RhsL=zeros((qsd+nsd+qsd+nsd)*NumElementNodes,1);
RhsG=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);

% Rhs local
RhsL(iH,1)=fH;
RhsL(ie,1)=fe;
RhsL(iM,1)=fM;
RhsL(iJ,1)=fJ;

% Rhs global
RhsG(iE,1)=fE;

end

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  Nodes,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementNodesPostHigh=size(RefElement.PostHighHigh.NodesCoordElem,1);
mu=Parameters.MagneticPermeability;
Sigma=Parameters.DampingFunction;
Xe=Nodes';
dt=Time.TimeStepSize;
BDFo=Time.BDFOrderEff;
alpha=Time.BDF1stDerEff;

% Get solution
He=reshape(SolutionLocal(:,1:qsd),[],1);
ee=reshape(SolutionLocal(:,qsd+(1:nsd)),[],1);
Me=reshape(SolutionLocal(:,qsd+nsd+(1:qsd)),[],1);
Holde=reshape(SolutionOld(:,1:qsd,:),[],BDFo);

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
else
  Hxe=He(nle1);
  Hye=He(nle2);
  Hze=He(nle3);
end
exe=ee(nle1);
eYe=ee(nle2);
if nsd==3
  eze=ee(nle3);
end
if nsd==2
  Mze=Me(nle1);
else
  Mxe=Me(nle1);
  Mye=Me(nle2);
  Mze=Me(nle3);
end
if nsd==2
  Holdze=Holde(nle1,:);
else
  Holdxe=Holde(nle1,:);
  Holdye=Holde(nle2,:);
  Holdze=Holde(nle3,:);
end
Sigmae=Sigma(Xe(nle1,1),Xe(nle1,2),Xe(nle1,3));
sigmaxe=Sigmae(:,1);
sigmaye=Sigmae(:,2);
if nsd==3
  sigmaze=Sigmae(:,3);
end

% Compute variables at Gauss points
if nsd==2
  Hzeg=Nle*Hze;
else
  Hxeg=Nle*Hxe;
  Hyeg=Nle*Hye;
  Hzeg=Nle*Hze;
end
exeg=Nlhe*exe;
eyeg=Nlhe*eYe;
if nsd==3
  ezeg=Nlhe*eze;
end
if nsd==2
  Mzeg=Nle*Mze;
else
  Mxeg=Nle*Mxe;
  Myeg=Nle*Mye;
  Mzeg=Nle*Mze;
end
if nsd==2
  Holdzeg=Nle*Holdze;
else
  Holdxeg=Nle*Holdxe;
  Holdyeg=Nle*Holdye;
  Holdzeg=Nle*Holdze;
end
sigmaxeg=Nle*sigmaxe;
sigmayeg=Nle*sigmaye;
if nsd==3
  sigmazeg=Nle*sigmaze;
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
else
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
else
  fp(ne1,1)=+NwezT*(mu/dt*Hyeg*alpha(1)+mu/dt*Holdyeg*alpha(2:BDFo+1,1))...
            -NweyT*(mu/dt*Hzeg*alpha(1)+mu/dt*Holdzeg*alpha(2:BDFo+1,1));
  fp(ne2,1)=+NwexT*(mu/dt*Hzeg*alpha(1)+mu/dt*Holdzeg*alpha(2:BDFo+1,1))...
            -NwezT*(mu/dt*Hxeg*alpha(1)+mu/dt*Holdxeg*alpha(2:BDFo+1,1));
  fp(ne3,1)=+NweyT*(mu/dt*Hxeg*alpha(1)+mu/dt*Holdxeg*alpha(2:BDFo+1,1))...
            -NwexT*(mu/dt*Hyeg*alpha(1)+mu/dt*Holdyeg*alpha(2:BDFo+1,1));
end

if nsd==2
  fp(ne1,1)=fp(ne1,1)-NweyT*(mu*(sigmaxeg+sigmayeg).*Hzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(mu*(sigmaxeg+sigmayeg).*Hzeg);
else
  fp(ne1,1)=fp(ne1,1)+NwezT*(mu*(sigmazeg+sigmaxeg-sigmayeg).*Hyeg)...
                     -NweyT*(mu*(sigmaxeg+sigmayeg-sigmazeg).*Hzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(mu*(sigmaxeg+sigmayeg-sigmazeg).*Hzeg)...
                     -NwezT*(mu*(sigmayeg+sigmazeg-sigmaxeg).*Hxeg);
  fp(ne3,1)=fp(ne3,1)+NweyT*(mu*(sigmayeg+sigmazeg-sigmaxeg).*Hxeg)...
                     -NwexT*(mu*(sigmazeg+sigmaxeg-sigmayeg).*Hyeg);
end

if nsd==2
  fp(ne1,1)=fp(ne1,1)-NweyT*(Mzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(Mzeg);
else
  fp(ne1,1)=fp(ne1,1)+NwezT*(Myeg)...
                     -NweyT*(Mzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(Mzeg)...
                     -NwezT*(Mxeg);
  fp(ne3,1)=fp(ne3,1)+NweyT*(Mxeg)...
                     -NwexT*(Myeg);
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
else
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
[Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf(1:nsd,:),nsd);

% Indices
nf1=FaceNodes(iFace,:);

% Compute variables at nodes
if nsd==2
  Hwzf=Hwze(nf1);
else
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
else
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
else
  Sxfg=real(ewyfg.*conj(Hwzfg)-ewzfg.*conj(Hwyfg));
  Syfg=real(ewzfg.*conj(Hwxfg)-ewxfg.*conj(Hwzfg));
  Szfg=real(ewxfg.*conj(Hwyfg)-ewyfg.*conj(Hwxfg));
end

% Compute radiation
if nsd==2
  Radiation=sum((Sxfg.*nx+Syfg.*ny).*wfg);
else
  Radiation=sum((Sxfg.*nx+Syfg.*ny+Szfg.*nz).*wfg);
end

end

%% Compute divergence element
function [MagneticFieldDiv,ElectricFieldDiv]=computeDivergenceElement(...
  Nodes,SolutionLocal,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
Xe=Nodes';

% Get solution
He=reshape(SolutionLocal(:,1:qsd),[],1);
ee=reshape(SolutionLocal(:,qsd+(1:nsd)),[],1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,~,~,pinvNe]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Hze=He(ne1);
else
  Hxe=He(ne1);
  Hye=He(ne2);
  Hze=He(ne3);
end
exe=ee(ne1);
eYe=ee(ne2);
if nsd==3
  eze=ee(ne3);
end

% Compute divergence of magnetic field
if nsd==2
  MagneticFieldDiv=zeros(NumElementNodes,1);
else
  MagneticFieldDiv=pinvNe*(Nex*Hxe+Ney*Hye+Nez*Hze);
end

% Compute divergence of electric field
if nsd==2
  ElectricFieldDiv=pinvNe*(Nex*exe+Ney*eYe);
else
  ElectricFieldDiv=pinvNe*(Nex*exe+Ney*eYe+Nez*eze);
end

% Add missing term for axisymmetric formulation
if nsd==2 && matchField(Parameters,'Axisymmetric','yes')
  ElectricFieldDiv=ElectricFieldDiv+pinvNe*(1./(Ne*Xe(:,2)).*(Ne*eYe));
end

end