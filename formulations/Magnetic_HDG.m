classdef Magnetic_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim+1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2+NumSpaceDim+1;
    
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
    
    %% Store results
    function [Results]=storeResults(~,iD,iST,Results,Block,~,Parameters,~,Time,~,Sizes)
      if iST==1
        Results(iD).Time=[];
        Results(iD).ScaledMagneticGradient=[];
        Results(iD).MagneticInduction=[];
        Results(iD).LagrangeMultiplier=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledMagneticGradient(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        1:Sizes(iD).NumVoigtComp);
      Results(iD).MagneticInduction(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+(1:Sizes(iD).NumSpaceDim));
      Results(iD).LagrangeMultiplier(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1);
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
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
eta=Parameters.MagneticDiffusivity;
zeta=0*Parameters.MagneticDiffusivity;
bD=Parameters.MagneticInduction;
qD=Parameters.LagrangeMultiplier;
sN=Parameters.PseudoTraction;
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
Ge=reshape(SolutionLocal(:,1:msd),[],1);
be=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
qe=reshape(SolutionLocal(:,msd+nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  bolde=reshape(SolutionOld(:,msd+(1:nsd),:),[],BDFo);
end

% Initialize lhs
KGG=zeros(msd*NumElementNodes,msd*NumElementNodes);
KGb=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KGB=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kbb=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kbq=zeros(nsd*NumElementNodes,NumElementNodes);
KbB=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kqq=zeros(NumElementNodes,NumElementNodes);
KqB=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KqQ=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KBG=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KBb=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KBq=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementNodes);
KBB=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KQQ=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fG=zeros(msd*NumElementNodes,1);
fb=zeros(nsd*NumElementNodes,1);
fq=zeros(NumElementNodes,1);
fB=zeros(nsd*NumElementFaces*NumFaceNodes,1);
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
  boldxe=bolde(ne1,:);
  boldye=bolde(ne2,:);
  if nsd==3
    boldze=bolde(ne3,:);
  end
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
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
CxeT=Cxe';
CyeT=Cye';
if nsd==3
  CzeT=Cze';
end

% Compute material matrix
if nsd==2
  Voigt1=(sqrt(2*(eta+zeta))+sqrt(2*eta))/2;
  Voigt2=(sqrt(2*(eta+zeta))-sqrt(2*eta))/2;
  Voigt3=sqrt(eta);
elseif nsd==3
  Voigt1=(sqrt(2*eta+3*zeta)+2*sqrt(2*eta))/3;
  Voigt2=(sqrt(2*eta+3*zeta)-1*sqrt(2*eta))/3;
  Voigt3=sqrt(eta);
end

% Compute lhs
KGG(ne1,ne1)=-Me;
KGG(ne2,ne2)=-Me;
KGG(ne3,ne3)=-Me;
if nsd==3
  KGG(ne4,ne4)=-Me;
  KGG(ne5,ne5)=-Me;
  KGG(ne6,ne6)=-Me;
end

if nsd==2
  KGb(ne1,ne1)=Voigt1*Cxe;
  KGb(ne2,ne1)=Voigt2*Cxe;
  KGb(ne3,ne1)=Voigt3*Cye;
  KGb(ne1,ne2)=Voigt2*Cye;
  KGb(ne2,ne2)=Voigt1*Cye;
  KGb(ne3,ne2)=Voigt3*Cxe;
elseif nsd==3
  KGb(ne1,ne1)=Voigt1*Cxe;
  KGb(ne2,ne1)=Voigt2*Cxe;
  KGb(ne3,ne1)=Voigt2*Cxe;
  KGb(ne4,ne1)=Voigt3*Cye;
  KGb(ne5,ne1)=Voigt3*Cze;
  KGb(ne1,ne2)=Voigt2*Cye;
  KGb(ne2,ne2)=Voigt1*Cye;
  KGb(ne3,ne2)=Voigt2*Cye;
  KGb(ne4,ne2)=Voigt3*Cxe;
  KGb(ne6,ne2)=Voigt3*Cze;
  KGb(ne1,ne3)=Voigt2*Cze;
  KGb(ne2,ne3)=Voigt2*Cze;
  KGb(ne3,ne3)=Voigt1*Cze;
  KGb(ne5,ne3)=Voigt3*Cxe;
  KGb(ne6,ne3)=Voigt3*Cye;
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

Kbq(ne1,ne1)=Kbq(ne1,ne1)+CxeT;
Kbq(ne2,ne1)=Kbq(ne2,ne1)+CyeT;
if nsd==3
  Kbq(ne3,ne1)=Kbq(ne3,ne1)+CzeT;
end

% Compute rhs
if nsd==2
  fG(ne1,1)=+NweT*(Gxxeg)...
            -NwexT*(Voigt1*bxeg)...
            -NweyT*(Voigt2*byeg);
  fG(ne2,1)=+NweT*(Gyyeg)...
            -NwexT*(Voigt2*bxeg)...
            -NweyT*(Voigt1*byeg);
  fG(ne3,1)=+NweT*(Gxyeg)...
            -NweyT*(Voigt3*bxeg)...
            -NwexT*(Voigt3*byeg);
elseif nsd==3
  fG(ne1,1)=+NweT*(Gxxeg)...
            -NwexT*(Voigt1*bxeg)...
            -NweyT*(Voigt2*byeg)...
            -NwezT*(Voigt2*bzeg);
  fG(ne2,1)=+NweT*(Gyyeg)...
            -NwexT*(Voigt2*bxeg)...
            -NweyT*(Voigt1*byeg)...
            -NwezT*(Voigt2*bzeg);
  fG(ne3,1)=+NweT*(Gzzeg)...
            -NwexT*(Voigt2*bxeg)...
            -NweyT*(Voigt2*byeg)...
            -NwezT*(Voigt1*bzeg);
  fG(ne4,1)=+NweT*(Gxyeg)...
            -NwexT*(Voigt3*byeg)...
            -NweyT*(Voigt3*bxeg);
  fG(ne5,1)=+NweT*(Gxzeg)...
            -NwexT*(Voigt3*bzeg)...
            -NwezT*(Voigt3*bxeg);
  fG(ne6,1)=+NweT*(Gyzeg)...
            -NweyT*(Voigt3*bzeg)...
            -NwezT*(Voigt3*byeg);
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
  fb(ne1,1)=fb(ne1,1)-NweT*(+Voigt1*(Nex*Gxxe)+Voigt2*(Nex*Gyye)...
                            +Voigt3*(Ney*Gxye)...
                            +(Nex*qe));
  fb(ne2,1)=fb(ne2,1)-NweT*(+Voigt2*(Ney*Gxxe)+Voigt1*(Ney*Gyye)...
                            +Voigt3*(Nex*Gxye)...
                            +(Ney*qe));
elseif nsd==3
  fb(ne1,1)=fb(ne1,1)-NweT*(+Voigt1*(Nex*Gxxe)+Voigt2*(Nex*Gyye)+Voigt2*(Nex*Gzze)...
                            +Voigt3*(Ney*Gxye)+Voigt3*(Nez*Gxze)...
                            +(Nex*qe));
  fb(ne2,1)=fb(ne2,1)-NweT*(+Voigt2*(Ney*Gxxe)+Voigt1*(Ney*Gyye)+Voigt2*(Ney*Gzze)...
                            +Voigt3*(Nex*Gxye)+Voigt3*(Nez*Gyze)...
                            +(Ney*qe));
  fb(ne3,1)=fb(ne3,1)-NweT*(+Voigt2*(Nez*Gxxe)+Voigt2*(Nez*Gyye)+Voigt1*(Nez*Gzze)...
                            +Voigt3*(Nex*Gxze)+Voigt3*(Ney*Gyze)...
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions('Face',RefElement,RefElement,Xf,nsd);
    
    % Check boundary
    isExterior=Faces.Exterior(iFace);
    isDirichlet_b_x=Faces.Dirichlet_b_x(iFace);
    isDirichlet_b_y=Faces.Dirichlet_b_y(iFace);
    if nsd==3; isDirichlet_b_z=Faces.Dirichlet_b_z(iFace); end
    isDirichlet_q=Faces.Dirichlet_q(iFace);
    isNeumann_s_x=Faces.Neumann_s_x(iFace);
    isNeumann_s_y=Faces.Neumann_s_y(iFace);
    if nsd==3; isNeumann_s_z=Faces.Neumann_s_z(iFace); end
    isDirichlet_b=isDirichlet_b_x || isDirichlet_b_y || (nsd==3 && isDirichlet_b_z);
    isNeumann_s=isNeumann_s_x || isNeumann_s_y || (nsd==3 && isNeumann_s_z);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nf4=nf3+NumElementNodes;
    nf5=nf4+NumElementNodes;
    nf6=nf5+NumElementNodes;
    nefU1=(iFace-1)*(nsd+1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
    nefU4=nefU3+NumFaceNodes;
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
      nefB1=nefB1(order);
      nefB2=nefB2(order);
      nefB3=nefB3(order);
      nefQ1=nefQ1(order);
    end
    
    % Compute variables at nodes
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
    Bxf=Ue(nefU1);
    Byf=Ue(nefU2);
    if nsd==3
      Bzf=Ue(nefU3);
    end
    if nsd==2
      Qf=Ue(nefU3);
    elseif nsd==3
      Qf=Ue(nefU4);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
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
    Bxfg=Nf*Bxf;
    Byfg=Nf*Byf;
    if nsd==3
      Bzfg=Nf*Bzf;
    end
    Qfg=Nf*Qf;
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
    
    % Compute common terms
    if isDirichlet_b_x
      Bxfg=bDxfg;
    end
    if isDirichlet_b_y
      Byfg=bDyfg;
    end
    if nsd==3 && isDirichlet_b_z
      Bzfg=bDzfg;
    end
    if isDirichlet_q
      Qfg=qDfg;
    end
    
    % Compute lhs
    Kbb(nf1,nf1)=Kbb(nf1,nf1)+tauB*Mf;
    Kbb(nf2,nf2)=Kbb(nf2,nf2)+tauB*Mf;
    if nsd==3
      Kbb(nf3,nf3)=Kbb(nf3,nf3)+tauB*Mf;
    end
    
    Kqq(nf1,nf1)=Kqq(nf1,nf1)-tauQ*Mf;
    
    if not(isDirichlet_b_x)
      if nsd==2
        KGB(nf1,nefB1)=KGB(nf1,nefB1)-Voigt1*Mfnx;
        KGB(nf2,nefB1)=KGB(nf2,nefB1)-Voigt2*Mfnx;
        KGB(nf3,nefB1)=KGB(nf3,nefB1)-Voigt3*Mfny;
      elseif nsd==3
        KGB(nf1,nefB1)=KGB(nf1,nefB1)-Voigt1*Mfnx;
        KGB(nf2,nefB1)=KGB(nf2,nefB1)-Voigt2*Mfnx;
        KGB(nf3,nefB1)=KGB(nf3,nefB1)-Voigt2*Mfnx;
        KGB(nf4,nefB1)=KGB(nf4,nefB1)-Voigt3*Mfny;
        KGB(nf5,nefB1)=KGB(nf5,nefB1)-Voigt3*Mfnz;
      end
    end
    
    if not(isDirichlet_b_y)
      if nsd==2
        KGB(nf1,nefB2)=KGB(nf1,nefB2)-Voigt2*Mfny;
        KGB(nf2,nefB2)=KGB(nf2,nefB2)-Voigt1*Mfny;
        KGB(nf3,nefB2)=KGB(nf3,nefB2)-Voigt3*Mfnx;
      elseif nsd==3
        KGB(nf1,nefB2)=KGB(nf1,nefB2)-Voigt2*Mfny;
        KGB(nf2,nefB2)=KGB(nf2,nefB2)-Voigt1*Mfny;
        KGB(nf3,nefB2)=KGB(nf3,nefB2)-Voigt2*Mfny;
        KGB(nf4,nefB2)=KGB(nf4,nefB2)-Voigt3*Mfnx;
        KGB(nf6,nefB2)=KGB(nf6,nefB2)-Voigt3*Mfnz;
      end
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      KGB(nf1,nefB3)=KGB(nf1,nefB3)-Voigt2*Mfnz;
      KGB(nf2,nefB3)=KGB(nf2,nefB3)-Voigt2*Mfnz;
      KGB(nf3,nefB3)=KGB(nf3,nefB3)-Voigt1*Mfnz;
      KGB(nf5,nefB3)=KGB(nf5,nefB3)-Voigt3*Mfnx;
      KGB(nf6,nefB3)=KGB(nf6,nefB3)-Voigt3*Mfny;
    end
    
    if not(isDirichlet_b_x) && isConvectiveTerm
      if nsd==2
        KbB(nf1,nefB1)=KbB(nf1,nefB1)+NwfT*((+vyfg.*ny).*Nf);
        KbB(nf2,nefB1)=KbB(nf2,nefB1)+NwfT*((-vyfg.*nx).*Nf);
      elseif nsd==3
        KbB(nf1,nefB1)=KbB(nf1,nefB1)+NwfT*((+vyfg.*ny+vzfg.*nz).*Nf);
        KbB(nf2,nefB1)=KbB(nf2,nefB1)+NwfT*((-vyfg.*nx).*Nf);
        KbB(nf3,nefB1)=KbB(nf3,nefB1)+NwfT*((-vzfg.*nx).*Nf);
      end
    end
    
    if not(isDirichlet_b_y) && isConvectiveTerm
      if nsd==2
        KbB(nf1,nefB2)=KbB(nf1,nefB2)+NwfT*((-vxfg.*ny).*Nf);
        KbB(nf2,nefB2)=KbB(nf2,nefB2)+NwfT*((+vxfg.*nx).*Nf);
      elseif nsd==3
        KbB(nf1,nefB2)=KbB(nf1,nefB2)+NwfT*((-vxfg.*ny).*Nf);
        KbB(nf2,nefB2)=KbB(nf2,nefB2)+NwfT*((+vxfg.*nx+vzfg.*nz).*Nf);
        KbB(nf3,nefB2)=KbB(nf3,nefB2)+NwfT*((-vzfg.*ny).*Nf);
      end
    end
    
    if nsd==3 && not(isDirichlet_b_z) && isConvectiveTerm
      KbB(nf1,nefB3)=KbB(nf1,nefB3)+NwfT*((-vxfg.*nz).*Nf);
      KbB(nf2,nefB3)=KbB(nf2,nefB3)+NwfT*((-vyfg.*nz).*Nf);
      KbB(nf3,nefB3)=KbB(nf3,nefB3)+NwfT*((+vxfg.*nx+vyfg.*ny).*Nf);
    end
    
    if not(isDirichlet_b_x)
      KbB(nf1,nefB1)=KbB(nf1,nefB1)-tauB*Mf;
    end
    
    if not(isDirichlet_b_y)
      KbB(nf2,nefB2)=KbB(nf2,nefB2)-tauB*Mf;
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      KbB(nf3,nefB3)=KbB(nf3,nefB3)-tauB*Mf;
    end
    
    if not(isDirichlet_b_x)
      KqB(nf1,nefB1)=KqB(nf1,nefB1)-Mfnx;
    end
    
    if not(isDirichlet_b_y)
      KqB(nf1,nefB2)=KqB(nf1,nefB2)-Mfny;
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      KqB(nf1,nefB3)=KqB(nf1,nefB3)-Mfnz;
    end
    
    if not(isDirichlet_q)
      KqQ(nf1,nefQ1)=KqQ(nf1,nefQ1)+tauQ*Mf;
    end
    
    if not(isDirichlet_b_x)
      KBb(nefB1,nf1)=KBb(nefB1,nf1)-tauB*Mf;
    end
    
    if not(isDirichlet_b_y)
      KBb(nefB2,nf2)=KBb(nefB2,nf2)-tauB*Mf;
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      KBb(nefB3,nf3)=KBb(nefB3,nf3)-tauB*Mf;
    end
    
    if not(isDirichlet_b_x)
      KBB(nefB1,nefB1)=KBB(nefB1,nefB1)+tauB*Mf;
    end
    
    if not(isDirichlet_b_y)
      KBB(nefB2,nefB2)=KBB(nefB2,nefB2)+tauB*Mf;
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      KBB(nefB3,nefB3)=KBB(nefB3,nefB3)+tauB*Mf;
    end
    
    if not(isDirichlet_q)
      KQQ(nefQ1,nefQ1)=KQQ(nefQ1,nefQ1)-tauQ*Mf;
    end
    
    % Compute rhs
    fb(nf1,1)=fb(nf1,1)-NwfT*(tauB*bxfg);
    fb(nf2,1)=fb(nf2,1)-NwfT*(tauB*byfg);
    if nsd==3
      fb(nf3,1)=fb(nf3,1)-NwfT*(tauB*bzfg);
    end
    
    fq(nf1,1)=fq(nf1,1)+NwfT*(tauQ*qfg);
    
    if nsd==2
      fG(nf1,1)=fG(nf1,1)+NwfT*(+Voigt1*nx.*Bxfg...
                                +Voigt2*ny.*Byfg);
      fG(nf2,1)=fG(nf2,1)+NwfT*(+Voigt2*nx.*Bxfg...
                                +Voigt1*ny.*Byfg);
      fG(nf3,1)=fG(nf3,1)+NwfT*(+Voigt3*ny.*Bxfg...
                                +Voigt3*nx.*Byfg);
    elseif nsd==3
      fG(nf1,1)=fG(nf1,1)+NwfT*(+Voigt1*nx.*Bxfg...
                                +Voigt2*ny.*Byfg...
                                +Voigt2*nz.*Bzfg);
      fG(nf2,1)=fG(nf2,1)+NwfT*(+Voigt2*nx.*Bxfg...
                                +Voigt1*ny.*Byfg...
                                +Voigt2*nz.*Bzfg);
      fG(nf3,1)=fG(nf3,1)+NwfT*(+Voigt2*nx.*Bxfg...
                                +Voigt2*ny.*Byfg...
                                +Voigt1*nz.*Bzfg);
      fG(nf4,1)=fG(nf4,1)+NwfT*(+Voigt3*nx.*Byfg...
                                +Voigt3*ny.*Bxfg);
      fG(nf5,1)=fG(nf5,1)+NwfT*(+Voigt3*nx.*Bzfg...
                                +Voigt3*nz.*Bxfg);
      fG(nf6,1)=fG(nf6,1)+NwfT*(+Voigt3*ny.*Bzfg...
                                +Voigt3*nz.*Byfg);
    end
    
    if isConvectiveTerm
      if nsd==2
        fb(nf1,1)=fb(nf1,1)-NwfT*((Bxfg.*vyfg-vxfg.*Byfg).*ny);
        fb(nf2,1)=fb(nf2,1)-NwfT*((Byfg.*vxfg-vyfg.*Bxfg).*nx);
      elseif nsd==3
        fb(nf1,1)=fb(nf1,1)-NwfT*((Bxfg.*vyfg-vxfg.*Byfg).*ny...
                                 +(Bxfg.*vzfg-vxfg.*Bzfg).*nz);
        fb(nf2,1)=fb(nf2,1)-NwfT*((Byfg.*vxfg-vyfg.*Bxfg).*nx...
                                 +(Byfg.*vzfg-vyfg.*Bzfg).*nz);
        fb(nf3,1)=fb(nf3,1)-NwfT*((Bzfg.*vxfg-vzfg.*Bxfg).*nx...
                                 +(Bzfg.*vyfg-vzfg.*Byfg).*ny);
      end
    end
    
    fb(nf1,1)=fb(nf1,1)+NwfT*(tauB*Bxfg);
    fb(nf2,1)=fb(nf2,1)+NwfT*(tauB*Byfg);
    if nsd==3
      fb(nf3,1)=fb(nf3,1)+NwfT*(tauB*Bzfg);
    end
    
    fq(nf1,1)=fq(nf1,1)+NwfT*(Bxfg.*nx+Byfg.*ny-tauQ*Qfg);
    if nsd==3
      fq(nf1,1)=fq(nf1,1)+NwfT*(Bzfg.*nz);
    end

    if not(isExterior) || isNeumann_s_x
      if nsd==2
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(+Voigt1*nx.*Gxxfg+Voigt2*nx.*Gyyfg...
                                      +Voigt3*ny.*Gxyfg...
                                      +qfg.*nx);
      elseif nsd==3
        fB(nefB1,1)=fB(nefB1,1)+NwfT*(+Voigt1*nx.*Gxxfg+Voigt2*nx.*Gyyfg+Voigt2*nx.*Gzzfg...
                                      +Voigt3*ny.*Gxyfg+Voigt3*nz.*Gxzfg...
                                      +qfg.*nx);
      end
    end
    
    if not(isExterior) || isNeumann_s_y
      if nsd==2
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(+Voigt2*ny.*Gxxfg+Voigt1*ny.*Gyyfg...
                                      +Voigt3*nx.*Gxyfg...
                                      +qfg.*ny);
      elseif nsd==3
        fB(nefB2,1)=fB(nefB2,1)+NwfT*(+Voigt2*ny.*Gxxfg+Voigt1*ny.*Gyyfg+Voigt2*ny.*Gzzfg...
                                      +Voigt3*nx.*Gxyfg+Voigt3*nz.*Gyzfg...
                                      +qfg.*ny);
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_s_z)
      fB(nefB3,1)=fB(nefB3,1)+NwfT*(+Voigt2*nz.*Gxxfg+Voigt2*nz.*Gyyfg+Voigt1*nz.*Gzzfg...
                                    +Voigt3*nx.*Gxzfg+Voigt3*ny.*Gyzfg...
                                    +qfg.*nz);
    end
    
    if not(isDirichlet_b_x)
      fB(nefB1,1)=fB(nefB1,1)+NwfT*(+tauB*(bxfg-Bxfg));
    end
    
    if not(isDirichlet_b_y)
      fB(nefB2,1)=fB(nefB2,1)+NwfT*(+tauB*(byfg-Byfg));
    end
    
    if nsd==3 && not(isDirichlet_b_z)
      fB(nefB3,1)=fB(nefB3,1)+NwfT*(+tauB*(bzfg-Bzfg));
    end
      
    if not(isDirichlet_q)
      fQ(nefQ1,1)=fQ(nefQ1,1)-NwfT*(tauQ*(qfg-Qfg));
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
    
    % Remove undetermination
    if isDirichlet_b_x
      KBB(nefB1,nefB1)=eye(NumFaceNodes);
    end
    if isDirichlet_b_y
      KBB(nefB2,nefB2)=eye(NumFaceNodes);
    end
    if nsd==3 && isDirichlet_b_z
      KBB(nefB3,nefB3)=eye(NumFaceNodes);
    end
    if isDirichlet_q
      KQQ(nefQ1,nefQ1)=eye(NumFaceNodes);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iG=1:msd*NumElementNodes;
ib=iG(end)+(1:nsd*NumElementNodes);
iq=ib(end)+(1:NumElementNodes);
iB=reshape((0:NumElementFaces-1)*(nsd+1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iQ=reshape((0:NumElementFaces-1)*(nsd+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+nsd*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+nsd+1)*NumElementNodes,(msd+nsd+1)*NumElementNodes);
LhsLG=zeros((msd+nsd+1)*NumElementNodes,(nsd+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd+1)*NumElementFaces*NumFaceNodes,(msd+nsd+1)*NumElementNodes);
LhsGG=zeros((nsd+1)*NumElementFaces*NumFaceNodes,(nsd+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+nsd+1)*NumElementNodes,1);
RhsG=zeros((nsd+1)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iG,iG)=KGG;
LhsLL(iG,ib)=KGb;
LhsLL(ib,iG)=KGb';
LhsLL(ib,ib)=Kbb;
LhsLL(ib,iq)=Kbq;
LhsLL(iq,ib)=Kbq';
LhsLL(iq,iq)=Kqq;

% Lhs local-global
LhsLG(iG,iB)=KGB;
LhsLG(ib,iB)=KbB;
LhsLG(iq,iB)=KqB;
LhsLG(iq,iQ)=KqQ;

% Rhs local
RhsL(iG,1)=fG;
RhsL(ib,1)=fb;
RhsL(iq,1)=fq;

% Lhs global-local
LhsGL(iB,iG)=KBG;
LhsGL(iB,ib)=KBb;
LhsGL(iB,iq)=KBq;
LhsGL(iQ,iq)=KqQ';

% Lhs global-global
LhsGG(iB,iB)=KBB;
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
  Nodes,Faces,SolutionGlobal,SolutionLocal,Parameters,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
eta=Parameters.MagneticDiffusivity;
zeta=0*Parameters.MagneticDiffusivity;
bD=Parameters.MagneticInduction;
Xe=Nodes';
t=Time.Time;

% Get solution
Ge=reshape(SolutionLocal(:,1:msd),[],1);
be=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
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

% Compute material matrix
if nsd==2
  Voigt1=(sqrt(2*(eta+zeta))+sqrt(2*eta))/2;
  Voigt2=(sqrt(2*(eta+zeta))-sqrt(2*eta))/2;
  Voigt3=sqrt(eta);
elseif nsd==3
  Voigt1=(sqrt(2*eta+3*zeta)+2*sqrt(2*eta))/3;
  Voigt2=(sqrt(2*eta+3*zeta)-1*sqrt(2*eta))/3;
  Voigt3=sqrt(eta);
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
  fp(ne1,1)=NwexT*(Gxxeg)+NweyT*(Gxyeg);
  fp(ne2,1)=NwexT*(Gxyeg)+NweyT*(Gyyeg);
elseif nsd==3
  fp(ne1,1)=NwexT*(Gxxeg)+NweyT*(Gxyeg)+NwezT*(Gxzeg);
  fp(ne2,1)=NwexT*(Gxyeg)+NweyT*(Gyyeg)+NwezT*(Gyzeg);
  fp(ne3,1)=NwexT*(Gxzeg)+NweyT*(Gyzeg)+NwezT*(Gzzeg);
end

ft(1,1)=Nw1eT*(bxeg);
ft(2,1)=Nw1eT*(byeg);
if nsd==3
  ft(3,1)=Nw1eT*(bzeg);
end

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions('Face',RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  Xfg=Nlf*Xf;
  
  % Check boundary
  isDirichlet_b_x=Faces.Dirichlet_b_x(iFace);
  isDirichlet_b_y=Faces.Dirichlet_b_y(iFace);
  if nsd==3; isDirichlet_b_z=Faces.Dirichlet_b_z(iFace); end
  isDirichlet_b=isDirichlet_b_x || isDirichlet_b_y || (nsd==3 && isDirichlet_b_z);
  
  % Indices
  nlefU1=(iFace-1)*(nsd+1)*NumFaceNodes+(1:NumFaceNodes);
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
  Bxf=Ue(nlefU1);
  Byf=Ue(nlefU2);
  if nsd==3
    Bzf=Ue(nlefU3);
  end
  
  % Compute variables at Gauss points
  Bxfg=Nlf*Bxf;
  Byfg=Nlf*Byf;
  if nsd==3
    Bzfg=Nlf*Bzf;
  end
  if isDirichlet_b
    bDfg=bD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    bDxfg=bDfg(:,1);
    bDyfg=bDfg(:,2);
    if nsd==3
      bDzfg=bDfg(:,3);
    end
  end
  
  % Compute common terms
  if isDirichlet_b_x
    Bxfg=bDxfg;
  end
  if isDirichlet_b_y
    Byfg=bDyfg;
  end
  if nsd==3 && isDirichlet_b_z
    Bzfg=bDzfg;
  end
  
  % Compute basic matrices
  Nw1fT=(wfg.*N1f)';
  
  % Compute rhs
  if nsd==2
    fr(1)=fr(1)+Nw1fT*(-Bxfg.*ny+Byfg.*nx);
  elseif nsd==3
    fr(1)=fr(1)+Nw1fT*(-Byfg.*nz+Bzfg.*ny);
    fr(2)=fr(2)+Nw1fT*(+Bxfg.*nz-Bzfg.*nx);
    fr(3)=fr(3)+Nw1fT*(-Bxfg.*ny+Byfg.*nx);
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