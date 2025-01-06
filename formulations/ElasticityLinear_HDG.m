classdef ElasticityLinear_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) NumSpaceDim;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2+NumSpaceDim;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) NumSpaceDim;

    % Discretization type
    DiscretizationType='HDG';

    % Time derivative order
    TimeDerOrder=2;
    
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
          Parameters(iD).Displacement(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        Block(iD,iD).SolutionOld(:,:,2)=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
          Parameters(iD).Displacement(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',-Time.TimeStepSize)];% FIX THIS!!!
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
            zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
            Parameters(iD).Displacement(...
            Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,~,Parameters,~,Faces,Time,...
        RefElement,Sizes)
      NodesElem=Elements(iD1).Nodes;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionLocalElem=Elements(iD1).SolutionLocal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        NodesElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        SolutionGlobalElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        iD2=setdiff(1:2,iD1);
        iEF=sub2ind([Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces],...
          Faces(iD1,iD2).Interface(:,1),Faces(iD1,iD2).Interface(:,2));
        NodesElemCoupled(iEF)=Elements(iD2).Nodes(Faces(iD1,iD2).Interface(:,3));
        SolutionGlobalElemCoupled(iEF)=Elements(iD2).SolutionGlobal(Faces(iD1,iD2).Interface(:,3));
      else
        NodesElemCoupled=double.empty(Sizes.NumElements,0);
        SolutionGlobalElemCoupled=double.empty(Sizes.NumElements,0);
      end
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      MatLocal=cell(Sizes(iD1).NumElements,1);
      VecLocal=cell(Sizes(iD1).NumElements,1);
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem,MatLocalElem,VecLocalElem]=...
          buildBlockElement(iD1,NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          NodesElemCoupled(iElem,:),SolutionGlobalElemCoupled(iElem,:),...
          Parameters,Time,RefElement,Sizes);
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
    function [Elements]=doPostProcess(~,Elements,~,Parameters,Faces,Time,RefElement,Sizes)
      NodesElem=Elements.Nodes;
      SolutionGlobalElem=Elements.SolutionGlobal;
      SolutionLocalElem=Elements.SolutionLocal;
      LhsPost=cell(Sizes.NumElements,1);
      RhsPost=cell(Sizes.NumElements,1);
      parfor iElem=1:Sizes.NumElements
        [LhsPostElem,RhsPostElem]=...
          doPostProcessElement(iElem,NodesElem{iElem},...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},...
          Parameters,Faces,Time,RefElement,Sizes);
        LhsPost{iElem}=LhsPostElem;
        RhsPost{iElem}=RhsPostElem;
      end
      Elements.LhsPost=LhsPost;
      Elements.RhsPost=RhsPost;
    end
    
    %% Do coupling
    function [Block]=doCoupling(~,iD1,iD2,Block,~,~,Parameters,Mesh,Faces,~,RefElement,Sizes)
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        LhsCoupCoef=zeros(Sizes(iD1).NumElementLhsCoupCoef(iD2),Sizes(iD1).NumFacesInterface(iD2));
        for iFaceInterface=1:Sizes(iD1).NumFacesInterface(iD2)
          [LhsCoupElem]=...
            doCouplingElement(iFaceInterface,iD1,iD2,Block,Parameters,Mesh,Faces,RefElement,Sizes);
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
    function [Results]=storeResults(~,iD,iST,Results,Block,~,Parameters,~,Time,~,Sizes)
      if iST==1
        Results(iD).Time=[];
        Results(iD).ScaledStrain=[];
        Results(iD).Displacement=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledStrain(:,:,iST)=Block(iD,iD).SolutionLocal(:,1:Sizes(iD).NumVoigtComp);
      Results(iD).Displacement(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        (1:Sizes(iD).NumSpaceDim));
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).DisplacementPost=Block(iD,iD).SolutionPost;
      end
    end
    
    %% Data for Paraview
    function [PointData,CellData]=dataForParaview(~,Results,Parameters,Mesh,Sizes,isPostProcess)
      
      if not(isPostProcess)
        
        % Write scaled strain
        Q=Results.ScaledStrain(:,:,end);
        PointData=sprintf('\nTENSORS Scaled_strain float\n');
        if Sizes.NumSpaceDim==2
          PointData=[PointData,sprintf(['%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f\n'],[Q';zeros(6,size(Q,1))])];
        elseif Sizes.NumSpaceDim==3
          PointData=[PointData,sprintf(['%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f\n'],[Q';zeros(3,size(Q,1))])];
        end
        
        % Write displacement
        u=Results.Displacement(:,:,end);
        PointData=[PointData,sprintf('\nVECTORS Displacement float\n')];
        if Sizes.NumSpaceDim==2
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',[u';zeros(1,size(u,1))])];
        elseif Sizes.NumSpaceDim==3
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',u')];
        end
        
      elseif isPostProcess
        
        % Write postprocessed displacement
        up=Results.DisplacementPost(:,:,end);
        PointData=sprintf('\nVECTORS Displacement_post float\n');
        if Sizes.NumSpaceDim==2
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',[up';zeros(1,size(up,1))])];
        elseif Sizes.NumSpaceDim==3
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',up')];
        end
        
      end
      
      % Write density
      rho=Parameters.Density;
      rho_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        rho_DG(iElem,:)=rho;
      end
      CellData=[sprintf('\nSCALARS Density float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',rho_DG')];
      
      % Write Young's modulus
      E=Parameters.YoungsModulus;
      E_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        Ce=Mesh.Elements(:,iElem);
        Xe=Mesh.Nodes(:,Ce);
        Xm=mean(Xe,2)';
        E_DG(iElem,:)=E(Xm(:,1),Xm(:,2),Xm(:,3));
      end
      CellData=[CellData,sprintf('\nSCALARS Youngs_modulus float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',E_DG')];
      
      % Write Poisson's ratio
      nu=Parameters.PoissonsRatio;
      nu_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        Ce=Mesh.Elements(:,iElem);
        Xe=Mesh.Nodes(:,Ce);
        Xm=mean(Xe,2)';
        nu_DG(iElem,:)=nu(Xm(:,1),Xm(:,2),Xm(:,3));
      end
      CellData=[CellData,sprintf('\nSCALARS Poissons_ratio float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',nu_DG')];
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  iD1,Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,NodesCoupled,SolutionGlobalCoupled,...
  Parameters,Time,RefElement,Sizes)

% Get general parameters
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes(iD1).NumSpaceDim;
msd=nsd*(nsd+1)/2;
NumElementNodes=Sizes(iD1).NumElementNodes;
NumElementFaces=Sizes(iD1).NumElementFaces;
NumFaceNodes=Sizes(iD1).NumFaceNodes;
rho=Parameters(iD1).Density;
uD=Parameters(iD1).Displacement;
tN=Parameters(iD1).Traction;
f=Parameters(iD1).Force;
Xe=Nodes';
Xem=sum(Xe(1:NumElementFaces,:),1)/NumElementFaces;
tauU=Parameters(iD1).StabDisplacement;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  beta=Time.BDF2ndDerEff;
end

% Get solution
Qe=reshape(SolutionLocal(:,1:msd),[],1);
ue=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  uolde=reshape(SolutionOld(:,msd+(1:nsd),:),[],BDFo+1);
end

% Initialize lhs
KQQ=zeros(msd*NumElementNodes,msd*NumElementNodes);
KQu=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KQU=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kuu=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KuU=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KUU=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);

% Initialize rhs
fQ=zeros(msd*NumElementNodes,1);
fu=zeros(nsd*NumElementNodes,1);
fU=zeros(nsd*NumElementFaces*NumFaceNodes,1);

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
  Qxxe=Qe(ne1);
  Qyye=Qe(ne2);
  Qxye=Qe(ne3);
elseif nsd==3
  Qxxe=Qe(ne1);
  Qyye=Qe(ne2);
  Qzze=Qe(ne3);
  Qxye=Qe(ne4);
  Qxze=Qe(ne5);
  Qyze=Qe(ne6);
end
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

% Compute variables at Gauss points
Xeg=Ne*Xe;
if nsd==2
  Qxxeg=Ne*Qxxe;
  Qyyeg=Ne*Qyye;
  Qxyeg=Ne*Qxye;
elseif nsd==3
  Qxxeg=Ne*Qxxe;
  Qyyeg=Ne*Qyye;
  Qzzeg=Ne*Qzze;
  Qxyeg=Ne*Qxye;
  Qxzeg=Ne*Qxze;
  Qyzeg=Ne*Qyze;
end
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
Cxe=NwexT*Ne;
Cye=NweyT*Ne;
if nsd==3
  Cze=NwezT*Ne;
end

% Compute stress and its linearization
[~,~,D]=computeStress('no','yes','no',Parameters(iD1),Xem,[],Sizes(iD1));
D_12=D^(1/2);
Voigt1=D_12(1,1);
Voigt2=D_12(2,1);
Voigt3=D_12(msd,msd);

% Compute lhs
KQQ(ne1,ne1)=-Me;
KQQ(ne2,ne2)=-Me;
KQQ(ne3,ne3)=-Me;
if nsd==3
  KQQ(ne4,ne4)=-Me;
  KQQ(ne5,ne5)=-Me;
  KQQ(ne6,ne6)=-Me;
end

if nsd==2
  KQu(ne1,ne1)=Voigt1*Cxe;
  KQu(ne2,ne1)=Voigt2*Cxe;
  KQu(ne3,ne1)=Voigt3*Cye;
  KQu(ne1,ne2)=Voigt2*Cye;
  KQu(ne2,ne2)=Voigt1*Cye;
  KQu(ne3,ne2)=Voigt3*Cxe;
elseif nsd==3
  KQu(ne1,ne1)=Voigt1*Cxe;
  KQu(ne2,ne1)=Voigt2*Cxe;
  KQu(ne3,ne1)=Voigt2*Cxe;
  KQu(ne4,ne1)=Voigt3*Cye;
  KQu(ne5,ne1)=Voigt3*Cze;
  KQu(ne1,ne2)=Voigt2*Cye;
  KQu(ne2,ne2)=Voigt1*Cye;
  KQu(ne3,ne2)=Voigt2*Cye;
  KQu(ne4,ne2)=Voigt3*Cxe;
  KQu(ne6,ne2)=Voigt3*Cze;
  KQu(ne1,ne3)=Voigt2*Cze;
  KQu(ne2,ne3)=Voigt2*Cze;
  KQu(ne3,ne3)=Voigt1*Cze;
  KQu(ne5,ne3)=Voigt3*Cxe;
  KQu(ne6,ne3)=Voigt3*Cye;
end

if isTimeDependent
  Kuu(ne1,ne1)=beta(1)/dt^2*rho*Me;
  Kuu(ne2,ne2)=beta(1)/dt^2*rho*Me;
  if nsd==3
    Kuu(ne3,ne3)=beta(1)/dt^2*rho*Me;
  end
end

% Compute rhs
if nsd==2
  fQ(ne1,1)=+NweT*(Qxxeg)...
            -NwexT*(Voigt1*uxeg)...
            -NweyT*(Voigt2*uyeg);
  fQ(ne2,1)=+NweT*(Qyyeg)...
            -NwexT*(Voigt2*uxeg)...
            -NweyT*(Voigt1*uyeg);
  fQ(ne3,1)=+NweT*(Qxyeg)...
            -NweyT*(Voigt3*uxeg)...
            -NwexT*(Voigt3*uyeg);
elseif nsd==3
  fQ(ne1,1)=+NweT*(Qxxeg)...
            -NwexT*(Voigt1*uxeg)...
            -NweyT*(Voigt2*uyeg)...
            -NwezT*(Voigt2*uzeg);
  fQ(ne2,1)=+NweT*(Qyyeg)...
            -NwexT*(Voigt2*uxeg)...
            -NweyT*(Voigt1*uyeg)...
            -NwezT*(Voigt2*uzeg);
  fQ(ne3,1)=+NweT*(Qzzeg)...
            -NwexT*(Voigt2*uxeg)...
            -NweyT*(Voigt2*uyeg)...
            -NwezT*(Voigt1*uzeg);
  fQ(ne4,1)=+NweT*(Qxyeg)...
            -NwexT*(Voigt3*uyeg)...
            -NweyT*(Voigt3*uxeg);
  fQ(ne5,1)=+NweT*(Qxzeg)...
            -NwexT*(Voigt3*uzeg)...
            -NwezT*(Voigt3*uxeg);
  fQ(ne6,1)=+NweT*(Qyzeg)...
            -NweyT*(Voigt3*uzeg)...
            -NwezT*(Voigt3*uyeg);
end

if isTimeDependent
  fu(ne1,1)=-NweT*(rho/dt^2*uxeg*beta(1)...
                  +rho/dt^2*uoldxeg*beta(2:BDFo+2,1));
  fu(ne2,1)=-NweT*(rho/dt^2*uyeg*beta(1)...
                  +rho/dt^2*uoldyeg*beta(2:BDFo+2,1));
  if nsd==3
    fu(ne3,1)=-NweT*(rho/dt^2*uzeg*beta(1)...
                    +rho/dt^2*uoldzeg*beta(2:BDFo+2,1));
  end
end

if nsd==2
  fu(ne1,1)=fu(ne1,1)-NweT*(+Voigt1*(Nex*Qxxe)...
                            +Voigt2*(Nex*Qyye)...
                            +Voigt3*(Ney*Qxye));
  fu(ne2,1)=fu(ne2,1)-NweT*(+Voigt2*(Ney*Qxxe)...
                            +Voigt1*(Ney*Qyye)...
                            +Voigt3*(Nex*Qxye));
elseif nsd==3
  fu(ne1,1)=fu(ne1,1)-NweT*(+Voigt1*(Nex*Qxxe)...
                            +Voigt2*(Nex*Qyye)...
                            +Voigt2*(Nex*Qzze)...
                            +Voigt3*(Ney*Qxye)...
                            +Voigt3*(Nez*Qxze));
  fu(ne2,1)=fu(ne2,1)-NweT*(+Voigt2*(Ney*Qxxe)...
                            +Voigt1*(Ney*Qyye)...
                            +Voigt2*(Ney*Qzze)...
                            +Voigt3*(Nex*Qxye)...
                            +Voigt3*(Nez*Qyze));
  fu(ne3,1)=fu(ne3,1)-NweT*(+Voigt2*(Nez*Qxxe)...
                            +Voigt2*(Nez*Qyye)...
                            +Voigt1*(Nez*Qzze)...
                            +Voigt3*(Nex*Qxze)...
                            +Voigt3*(Ney*Qyze));
end

fu(ne1,1)=fu(ne1,1)+NweT*(fxeg);
fu(ne2,1)=fu(ne2,1)+NweT*(fyeg);
if nsd==3
  fu(ne3,1)=fu(ne3,1)+NweT*(fzeg);
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
    nf4=nf3+NumElementNodes;
    nf5=nf4+NumElementNodes;
    nf6=nf5+NumElementNodes;
    nef1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nef2=nef1+NumFaceNodes;
    nef3=nef2+NumFaceNodes;
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters(iD1).Degree,Node2Match1stNode1);
      nef1=nef1(order);
      nef2=nef2(order);
      nef3=nef3(order);
    end
    
    % Compute variables at nodes
    if nsd==2
      Qxxf=Qxxe(nf1);
      Qyyf=Qyye(nf1);
      Qxyf=Qxye(nf1);
    elseif nsd==3
      Qxxf=Qxxe(nf1);
      Qyyf=Qyye(nf1);
      Qzzf=Qzze(nf1);
      Qxyf=Qxye(nf1);
      Qxzf=Qxze(nf1);
      Qyzf=Qyze(nf1);
    end
    uxf=uxe(nf1);
    uyf=uye(nf1);
    if nsd==3
      uzf=uze(nf1);
    end
    Uxf=Ue(nef1);
    Uyf=Ue(nef2);
    if nsd==3
      Uzf=Ue(nef3);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    if nsd==2
      Qxxfg=Nf*Qxxf;
      Qyyfg=Nf*Qyyf;
      Qxyfg=Nf*Qxyf;
    elseif nsd==3
      Qxxfg=Nf*Qxxf;
      Qyyfg=Nf*Qyyf;
      Qzzfg=Nf*Qzzf;
      Qxyfg=Nf*Qxyf;
      Qxzfg=Nf*Qxzf;
      Qyzfg=Nf*Qyzf;
    end
    uxfg=Nf*uxf;
    uyfg=Nf*uyf;
    if nsd==3
      uzfg=Nf*uzf;
    end
    Uxfg=Nf*Uxf;
    Uyfg=Nf*Uyf;
    if nsd==3
      Uzfg=Nf*Uzf;
    end
    if isDirichlet
      uDfg=uD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      uDxfg=uDfg(:,1);
      uDyfg=uDfg(:,2);
      if nsd==3
        uDzfg=uDfg(:,3);
      end
    elseif isNeumann
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
    if isInterface
      % Get general parameters
      iD2=setdiff(1:2,iD1);
      iFace2=Faces.Interface(2,iFace);
      NumElementNodes2=Sizes(iD2).NumElementNodes;
      NumElementFaces2=Sizes(iD2).NumElementFaces;
      X2e=NodesCoupled{iFace}';
      X2em=sum(X2e(1:NumElementFaces2,:),1)/NumElementFaces2;
      gamma=Parameters(iD2).NitschePenalty;
      
      % Get solution
      u2e=reshape(SolutionGlobalCoupled{iFace},[],1);
      
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
          computeStress('yes','no','no',Parameters(iD2),X2em,...
          [F2xxfg,F2xyfg,F2yxfg,F2yyfg],Sizes(iD2));
        s2xxfg=s2fg(:,1); s2xyfg=s2fg(:,2);
        s2yxfg=s2fg(:,3); s2yyfg=s2fg(:,4);
      elseif nsd==3
        [s2fg]=...
          computeStress('yes','no','no',Parameters(iD2),X2em,...
          [F2xxfg,F2xyfg,F2xzfg,F2yxfg,F2yyfg,F2yzfg,F2zxfg,F2zyfg,F2zzfg],Sizes(iD2));
        s2xxfg=s2fg(:,1); s2xyfg=s2fg(:,2); s2xzfg=s2fg(:,3);
        s2yxfg=s2fg(:,4); s2yyfg=s2fg(:,5); s2yzfg=s2fg(:,6);
        s2zxfg=s2fg(:,7); s2zyfg=s2fg(:,8); s2zzfg=s2fg(:,9);
      end
      
      % Compute basic matrices
      Nw12fT=(w12fg.*N12f)';
    end
    % ----------------------------------------------------------------------------------------------
    
    % Compute lhs
    Kuu(nf1,nf1)=Kuu(nf1,nf1)+tauU*Mf;
    Kuu(nf2,nf2)=Kuu(nf2,nf2)+tauU*Mf;
    if nsd==3
      Kuu(nf3,nf3)=Kuu(nf3,nf3)+tauU*Mf;
    end
    
    if not(isDirichlet)
      if nsd==2
        KQU(nf1,nef1)=KQU(nf1,nef1)-Voigt1*Mfnx;
        KQU(nf2,nef1)=KQU(nf2,nef1)-Voigt2*Mfnx;
        KQU(nf3,nef1)=KQU(nf3,nef1)-Voigt3*Mfny;
        KQU(nf1,nef2)=KQU(nf1,nef2)-Voigt2*Mfny;
        KQU(nf2,nef2)=KQU(nf2,nef2)-Voigt1*Mfny;
        KQU(nf3,nef2)=KQU(nf3,nef2)-Voigt3*Mfnx;
      elseif nsd==3
        KQU(nf1,nef1)=KQU(nf1,nef1)-Voigt1*Mfnx;
        KQU(nf2,nef1)=KQU(nf2,nef1)-Voigt2*Mfnx;
        KQU(nf3,nef1)=KQU(nf3,nef1)-Voigt2*Mfnx;
        KQU(nf4,nef1)=KQU(nf4,nef1)-Voigt3*Mfny;
        KQU(nf5,nef1)=KQU(nf5,nef1)-Voigt3*Mfnz;
        KQU(nf1,nef2)=KQU(nf1,nef2)-Voigt2*Mfny;
        KQU(nf2,nef2)=KQU(nf2,nef2)-Voigt1*Mfny;
        KQU(nf3,nef2)=KQU(nf3,nef2)-Voigt2*Mfny;
        KQU(nf4,nef2)=KQU(nf4,nef2)-Voigt3*Mfnx;
        KQU(nf6,nef2)=KQU(nf6,nef2)-Voigt3*Mfnz;
        KQU(nf1,nef3)=KQU(nf1,nef3)-Voigt2*Mfnz;
        KQU(nf2,nef3)=KQU(nf2,nef3)-Voigt2*Mfnz;
        KQU(nf3,nef3)=KQU(nf3,nef3)-Voigt1*Mfnz;
        KQU(nf5,nef3)=KQU(nf5,nef3)-Voigt3*Mfnx;
        KQU(nf6,nef3)=KQU(nf6,nef3)-Voigt3*Mfny;
      end
      
      KuU(nf1,nef1)=KuU(nf1,nef1)-tauU*Mf;
      KuU(nf2,nef2)=KuU(nf2,nef2)-tauU*Mf;
      if nsd==3
        KuU(nf3,nef3)=KuU(nf3,nef3)-tauU*Mf;
      end
      
      KUU(nef1,nef1)=KUU(nef1,nef1)+tauU*Mf;
      KUU(nef2,nef2)=KUU(nef2,nef2)+tauU*Mf;
      if nsd==3
        KUU(nef3,nef3)=KUU(nef3,nef3)+tauU*Mf;
      end
      
      if isInterface
        KUU(nef1,nef1)=KUU(nef1,nef1)+gamma/h*Mf;
        KUU(nef2,nef2)=KUU(nef2,nef2)+gamma/h*Mf;
        if nsd==3
          KUU(nef3,nef3)=KUU(nef3,nef3)+gamma/h*Mf;
        end
      end
    end
    
    % Compute rhs
    fu(nf1,1)=fu(nf1,1)-NwfT*(tauU*uxfg);
    fu(nf2,1)=fu(nf2,1)-NwfT*(tauU*uyfg);
    if nsd==3
      fu(nf3,1)=fu(nf3,1)-NwfT*(tauU*uzfg);
    end
    
    if isDirichlet
      if nsd==2
        fQ(nf1,1)=fQ(nf1,1)+NwfT*(+Voigt1*nx.*uDxfg+Voigt2*ny.*uDyfg);
        fQ(nf2,1)=fQ(nf2,1)+NwfT*(+Voigt2*nx.*uDxfg+Voigt1*ny.*uDyfg);
        fQ(nf3,1)=fQ(nf3,1)+NwfT*(+Voigt3*ny.*uDxfg+Voigt3*nx.*uDyfg);
      elseif nsd==3
        fQ(nf1,1)=fQ(nf1,1)+NwfT*(+Voigt1*nx.*uDxfg+Voigt2*ny.*uDyfg+Voigt2*nz.*uDzfg);
        fQ(nf2,1)=fQ(nf2,1)+NwfT*(+Voigt2*nx.*uDxfg+Voigt1*ny.*uDyfg+Voigt2*nz.*uDzfg);
        fQ(nf3,1)=fQ(nf3,1)+NwfT*(+Voigt2*nx.*uDxfg+Voigt2*ny.*uDyfg+Voigt1*nz.*uDzfg);
        fQ(nf4,1)=fQ(nf4,1)+NwfT*(+Voigt3*nx.*uDyfg+Voigt3*ny.*uDxfg);
        fQ(nf5,1)=fQ(nf5,1)+NwfT*(+Voigt3*nx.*uDzfg+Voigt3*nz.*uDxfg);
        fQ(nf6,1)=fQ(nf6,1)+NwfT*(+Voigt3*ny.*uDzfg+Voigt3*nz.*uDyfg);
      end
      
      fu(nf1,1)=fu(nf1,1)+NwfT*(tauU*uDxfg);
      fu(nf2,1)=fu(nf2,1)+NwfT*(tauU*uDyfg);
      if nsd==3
        fu(nf3,1)=fu(nf3,1)+NwfT*(tauU*uDzfg);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        fQ(nf1,1)=fQ(nf1,1)+NwfT*(+Voigt1*nx.*Uxfg+Voigt2*ny.*Uyfg);
        fQ(nf2,1)=fQ(nf2,1)+NwfT*(+Voigt2*nx.*Uxfg+Voigt1*ny.*Uyfg);
        fQ(nf3,1)=fQ(nf3,1)+NwfT*(+Voigt3*ny.*Uxfg+Voigt3*nx.*Uyfg);
      elseif nsd==3
        fQ(nf1,1)=fQ(nf1,1)+NwfT*(+Voigt1*nx.*Uxfg+Voigt2*ny.*Uyfg+Voigt2*nz.*Uzfg);
        fQ(nf2,1)=fQ(nf2,1)+NwfT*(+Voigt2*nx.*Uxfg+Voigt1*ny.*Uyfg+Voigt2*nz.*Uzfg);
        fQ(nf3,1)=fQ(nf3,1)+NwfT*(+Voigt2*nx.*Uxfg+Voigt2*ny.*Uyfg+Voigt1*nz.*Uzfg);
        fQ(nf4,1)=fQ(nf4,1)+NwfT*(+Voigt3*nx.*Uyfg+Voigt3*ny.*Uxfg);
        fQ(nf5,1)=fQ(nf5,1)+NwfT*(+Voigt3*nx.*Uzfg+Voigt3*nz.*Uxfg);
        fQ(nf6,1)=fQ(nf6,1)+NwfT*(+Voigt3*ny.*Uzfg+Voigt3*nz.*Uyfg);
      end
      
      fu(nf1,1)=fu(nf1,1)+NwfT*(tauU*Uxfg);
      fu(nf2,1)=fu(nf2,1)+NwfT*(tauU*Uyfg);
      if nsd==3
        fu(nf3,1)=fu(nf3,1)+NwfT*(tauU*Uzfg);
      end
      
      if nsd==2
        fU(nef1,1)=fU(nef1,1)+NwfT*(+Voigt1*nx.*Qxxfg+Voigt2*nx.*Qyyfg...
                                    +Voigt3*ny.*Qxyfg...
                                    +tauU*(uxfg-Uxfg));
        fU(nef2,1)=fU(nef2,1)+NwfT*(+Voigt2*ny.*Qxxfg+Voigt1*ny.*Qyyfg...
                                    +Voigt3*nx.*Qxyfg...
                                    +tauU*(uyfg-Uyfg));
      elseif nsd==3
        fU(nef1,1)=fU(nef1,1)+NwfT*(+Voigt1*nx.*Qxxfg+Voigt2*nx.*Qyyfg+Voigt2*nx.*Qzzfg...
                                    +Voigt3*ny.*Qxyfg+Voigt3*nz.*Qxzfg...
                                    +tauU*(uxfg-Uxfg));
        fU(nef2,1)=fU(nef2,1)+NwfT*(+Voigt2*ny.*Qxxfg+Voigt1*ny.*Qyyfg+Voigt2*ny.*Qzzfg...
                                    +Voigt3*nx.*Qxyfg+Voigt3*nz.*Qyzfg...
                                    +tauU*(uyfg-Uyfg));
        fU(nef3,1)=fU(nef3,1)+NwfT*(+Voigt2*nz.*Qxxfg+Voigt2*nz.*Qyyfg+Voigt1*nz.*Qzzfg...
                                    +Voigt3*nx.*Qxzfg+Voigt3*ny.*Qyzfg...
                                    +tauU*(uzfg-Uzfg));
      end
    end
    
    if isInterface
      if nsd==2
        fU(nef1,1)=fU(nef1,1)-Nw12fT*(+s2xxfg.*(-n12x)+s2xyfg.*(-n12y)...
                                      -gamma/h*u2xfg);
        fU(nef2,1)=fU(nef2,1)-Nw12fT*(+s2yxfg.*(-n12x)+s2yyfg.*(-n12y)....
                                      -gamma/h*u2yfg);
      elseif nsd==3
        fU(nef1,1)=fU(nef1,1)-Nw12fT*(+s2xxfg.*(-n12x)+s2xyfg.*(-n12y)+s2xzfg.*(-n12z)...
                                      -gamma/h*u2xfg);
        fU(nef2,1)=fU(nef2,1)-Nw12fT*(+s2yxfg.*(-n12x)+s2yyfg.*(-n12y)+s2yzfg.*(-n12z)...
                                      -gamma/h*u2yfg);
        fU(nef3,1)=fU(nef3,1)-Nw12fT*(+s2zxfg.*(-n12x)+s2zyfg.*(-n12y)+s2zzfg.*(-n12z)...
                                      -gamma/h*u2zfg);
      end
      
      fU(nef1,1)=fU(nef1,1)-NwfT*(gamma/h*Uxfg);
      fU(nef2,1)=fU(nef2,1)-NwfT*(gamma/h*Uyfg);
      if nsd==3
        fU(nef3,1)=fU(nef3,1)-NwfT*(gamma/h*Uzfg);
      end
    end
    
    if isNeumann
      fU(nef1,1)=fU(nef1,1)+NwfT*(tNxfg);
      fU(nef2,1)=fU(nef2,1)+NwfT*(tNyfg);
      if nsd==3
        fU(nef3,1)=fU(nef3,1)+NwfT*(tNzfg);
      end
    end
    
    % Remove undetermination
    if isDirichlet
      KUU(nef1,nef1)=eye(NumFaceNodes);
      KUU(nef2,nef2)=eye(NumFaceNodes);
      if nsd==3
        KUU(nef3,nef3)=eye(NumFaceNodes);
      end
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iQ=1:msd*NumElementNodes;
iu=iQ(end)+(1:nsd*NumElementNodes);
iU=1:nsd*NumElementFaces*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+nsd)*NumElementNodes,(msd+nsd)*NumElementNodes);
LhsLG=zeros((msd+nsd)*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
LhsGL=zeros(nsd*NumElementFaces*NumFaceNodes,(msd+nsd)*NumElementNodes);
LhsGG=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+nsd)*NumElementNodes,1);
RhsG=zeros(nsd*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iQ,iQ)=KQQ;
LhsLL(iQ,iu)=KQu;
LhsLL(iu,iQ)=KQu';
LhsLL(iu,iu)=Kuu;

% Lhs local-global
LhsLG(iQ,iU)=KQU;
LhsLG(iu,iU)=KuU;

% Rhs local
RhsL(iQ,1)=fQ;
RhsL(iu,1)=fu;

% Lhs global-local
LhsGL(iU,iQ)=KQU';
LhsGL(iU,iu)=KuU';

% Lhs global-global
LhsGG(iU,iU)=KUU;

% Rhs global
RhsG(iU,1)=fU;

% Matrix and vector for local problem
MatVecLocal=LhsLL\[LhsLG,RhsL];

% Extract matrix for local problem
MatLocal=MatVecLocal(:,1:end-1);

% Extract vector for local problem
VecLocal=MatVecLocal(:,end);

% Lhs for global problem
LhsGlobal=LhsGG-LhsGL*MatLocal;
LhsGlobal=(LhsGlobal+LhsGlobal.')/2;

% Rhs for global problem
RhsGlobal=RhsG-LhsGL*VecLocal;

end

%% Do coupling element
function [LhsCoup]=doCouplingElement(...
  iFaceInterface,iD1,iD2,Block,Parameters,Mesh,Faces,RefElement,Sizes)       

% Get general parameters
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
X2e=Mesh(iD2).Nodes(:,C2e)';
X2em=sum(X2e(1:NumElementFaces2,:),1)/NumElementFaces2;
gamma=Parameters(iD2).NitschePenalty;

% Get solution
u2e=reshape(Block(iD2,iD2).SolutionGlobal(C2e,:),[],1);

% Initialize lhs
KU1u2=zeros(nsd*NumElementFaces1*NumFaceNodes1,nsd*NumElementNodes2);

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
n1ef1=(iFace1-1)*nsd*NumFaceNodes1+(1:NumFaceNodes1);
n1ef2=n1ef1+NumFaceNodes1;
n1ef3=n1ef2+NumFaceNodes1;
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
    computeStress('no','yes','no',Parameters(iD2),X2em,...
    [F2xxfg,F2xyfg,F2yxfg,F2yyfg],Sizes(iD2));
  ds2xxdF2xxfg=ds2dF2fg(:,1);  ds2xxdF2xyfg=ds2dF2fg(:,2);  ds2xxdF2yxfg=ds2dF2fg(:,3);  ds2xxdF2yyfg=ds2dF2fg(:,4);
                               ds2xydF2xyfg=ds2dF2fg(:,6);  ds2xydF2yxfg=ds2dF2fg(:,7);  ds2xydF2yyfg=ds2dF2fg(:,8);
                                                            ds2yxdF2yxfg=ds2dF2fg(:,11); ds2yxdF2yyfg=ds2dF2fg(:,12);
                                                                                         ds2yydF2yyfg=ds2dF2fg(:,16);
elseif nsd==3
  [~,ds2dF2fg]=...
    computeStress('no','yes','no',Parameters(iD2),X2em,...
    [F2xxfg,F2xyfg,F2xzfg,F2yxfg,F2yyfg,F2yzfg,F2zxfg,F2zyfg,F2zzfg],Sizes(iD2));
  ds2xxdF2xxfg=ds2dF2fg(:,1);  ds2xxdF2xyfg=ds2dF2fg(:,2);  ds2xxdF2xzfg=ds2dF2fg(:,3);  ds2xxdF2yxfg=ds2dF2fg(:,4);  ds2xxdF2yyfg=ds2dF2fg(:,5);  ds2xxdF2yzfg=ds2dF2fg(:,6);  ds2xxdF2zxfg=ds2dF2fg(:,7);  ds2xxdF2zyfg=ds2dF2fg(:,8);  ds2xxdF2zzfg=ds2dF2fg(:,9);
                               ds2xydF2xyfg=ds2dF2fg(:,11); ds2xydF2xzfg=ds2dF2fg(:,12); ds2xydF2yxfg=ds2dF2fg(:,13); ds2xydF2yyfg=ds2dF2fg(:,14); ds2xydF2yzfg=ds2dF2fg(:,15); ds2xydF2zxfg=ds2dF2fg(:,16); ds2xydF2zyfg=ds2dF2fg(:,17); ds2xydF2zzfg=ds2dF2fg(:,18);
                                                            ds2xzdF2xzfg=ds2dF2fg(:,21); ds2xzdF2yxfg=ds2dF2fg(:,22); ds2xzdF2yyfg=ds2dF2fg(:,23); ds2xzdF2yzfg=ds2dF2fg(:,24); ds2xzdF2zxfg=ds2dF2fg(:,25); ds2xzdF2zyfg=ds2dF2fg(:,26); ds2xzdF2zzfg=ds2dF2fg(:,27);
                                                                                         ds2yxdF2yxfg=ds2dF2fg(:,31); ds2yxdF2yyfg=ds2dF2fg(:,32); ds2yxdF2yzfg=ds2dF2fg(:,33); ds2yxdF2zxfg=ds2dF2fg(:,34); ds2yxdF2zyfg=ds2dF2fg(:,35); ds2yxdF2zzfg=ds2dF2fg(:,36);
                                                                                                                      ds2yydF2yyfg=ds2dF2fg(:,41); ds2yydF2yzfg=ds2dF2fg(:,42); ds2yydF2zxfg=ds2dF2fg(:,43); ds2yydF2zyfg=ds2dF2fg(:,44); ds2yydF2zzfg=ds2dF2fg(:,45);
                                                                                                                                                   ds2yzdF2yzfg=ds2dF2fg(:,51); ds2yzdF2zxfg=ds2dF2fg(:,52); ds2yzdF2zyfg=ds2dF2fg(:,53); ds2yzdF2zzfg=ds2dF2fg(:,54);
                                                                                                                                                                                ds2zxdF2zxfg=ds2dF2fg(:,61); ds2zxdF2zyfg=ds2dF2fg(:,62); ds2zxdF2zzfg=ds2dF2fg(:,63);
                                                                                                                                                                                                             ds2zydF2zyfg=ds2dF2fg(:,71); ds2zydF2zzfg=ds2dF2fg(:,72);
                                                                                                                                                                                                                                          ds2zzdF2zzfg=ds2dF2fg(:,81);
end

% Compute lhs
if nsd==2
  KU1u2(n1ef1,n2e1)=KU1u2(n1ef1,n2e1)...
            +Nw12fT*((+ds2xxdF2xxfg.*(-n12x)+ds2xxdF2xyfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2xxdF2xyfg.*(-n12x)+ds2xydF2xyfg.*(-n12y)).*N21yf);
  KU1u2(n1ef1,n2e2)=KU1u2(n1ef1,n2e2)...
            +Nw12fT*((+ds2xxdF2yxfg.*(-n12x)+ds2xydF2yxfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2xxdF2yyfg.*(-n12x)+ds2xydF2yyfg.*(-n12y)).*N21yf);
  KU1u2(n1ef2,n2e1)=KU1u2(n1ef2,n2e1)...
            +Nw12fT*((+ds2xxdF2yxfg.*(-n12x)+ds2xxdF2yyfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2xydF2yxfg.*(-n12x)+ds2xydF2yyfg.*(-n12y)).*N21yf);
  KU1u2(n1ef2,n2e2)=KU1u2(n1ef2,n2e2)...
            +Nw12fT*((+ds2yxdF2yxfg.*(-n12x)+ds2yxdF2yyfg.*(-n12y)).*N21xf)...
            +Nw12fT*((+ds2yxdF2yyfg.*(-n12x)+ds2yydF2yyfg.*(-n12y)).*N21yf);
elseif nsd==3
  KU1u2(n1ef1,n2e1)=KU1u2(n1ef1,n2e1)...
            +Nw12fT*((+ds2xxdF2xxfg.*(-n12x)+ds2xxdF2xyfg.*(-n12y)+ds2xxdF2xzfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xxdF2xyfg.*(-n12x)+ds2xydF2xyfg.*(-n12y)+ds2xydF2xzfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xxdF2xzfg.*(-n12x)+ds2xydF2xzfg.*(-n12y)+ds2xzdF2xzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef1,n2e2)=KU1u2(n1ef1,n2e2)...
            +Nw12fT*((+ds2xxdF2yxfg.*(-n12x)+ds2xydF2yxfg.*(-n12y)+ds2xzdF2yxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xxdF2yyfg.*(-n12x)+ds2xydF2yyfg.*(-n12y)+ds2xzdF2yyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xxdF2yzfg.*(-n12x)+ds2xydF2yzfg.*(-n12y)+ds2xzdF2yzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef1,n2e3)=KU1u2(n1ef1,n2e3)...
            +Nw12fT*((+ds2xxdF2zxfg.*(-n12x)+ds2xydF2zxfg.*(-n12y)+ds2xzdF2zxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xxdF2zyfg.*(-n12x)+ds2xydF2zyfg.*(-n12y)+ds2xzdF2zyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xxdF2zzfg.*(-n12x)+ds2xydF2zzfg.*(-n12y)+ds2xzdF2zzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef2,n2e1)=KU1u2(n1ef2,n2e1)...
            +Nw12fT*((+ds2xxdF2yxfg.*(-n12x)+ds2xxdF2yyfg.*(-n12y)+ds2xxdF2yzfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xydF2yxfg.*(-n12x)+ds2xydF2yyfg.*(-n12y)+ds2xydF2yzfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xzdF2yxfg.*(-n12x)+ds2xzdF2yyfg.*(-n12y)+ds2xzdF2yzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef2,n2e2)=KU1u2(n1ef2,n2e2)...
            +Nw12fT*((+ds2yxdF2yxfg.*(-n12x)+ds2yxdF2yyfg.*(-n12y)+ds2yxdF2yzfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2yxdF2yyfg.*(-n12x)+ds2yydF2yyfg.*(-n12y)+ds2yydF2yzfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2yxdF2yzfg.*(-n12x)+ds2yydF2yzfg.*(-n12y)+ds2yzdF2yzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef2,n2e3)=KU1u2(n1ef2,n2e3)...
            +Nw12fT*((+ds2yxdF2zxfg.*(-n12x)+ds2yydF2zxfg.*(-n12y)+ds2yzdF2zxfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2yxdF2zyfg.*(-n12x)+ds2yydF2zyfg.*(-n12y)+ds2yzdF2zyfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2yxdF2zzfg.*(-n12x)+ds2yydF2zzfg.*(-n12y)+ds2yzdF2zzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef3,n2e1)=KU1u2(n1ef3,n2e1)...
            +Nw12fT*((+ds2xxdF2zxfg.*(-n12x)+ds2xxdF2zyfg.*(-n12y)+ds2xxdF2zzfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2xydF2zxfg.*(-n12x)+ds2xydF2zyfg.*(-n12y)+ds2xydF2zzfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2xzdF2zxfg.*(-n12x)+ds2xzdF2zyfg.*(-n12y)+ds2xzdF2zzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef3,n2e2)=KU1u2(n1ef3,n2e2)...
            +Nw12fT*((+ds2yxdF2zxfg.*(-n12x)+ds2yxdF2zyfg.*(-n12y)+ds2yxdF2zzfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2yydF2zxfg.*(-n12x)+ds2yydF2zyfg.*(-n12y)+ds2yydF2zzfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2yzdF2zxfg.*(-n12x)+ds2yzdF2zyfg.*(-n12y)+ds2yzdF2zzfg.*(-n12z)).*N21zf);
  KU1u2(n1ef3,n2e3)=KU1u2(n1ef3,n2e3)...
            +Nw12fT*((+ds2zxdF2zxfg.*(-n12x)+ds2zxdF2zyfg.*(-n12y)+ds2zxdF2zzfg.*(-n12z)).*N21xf)...
            +Nw12fT*((+ds2zxdF2zyfg.*(-n12x)+ds2zydF2zyfg.*(-n12y)+ds2zydF2zzfg.*(-n12z)).*N21yf)...
            +Nw12fT*((+ds2zxdF2zzfg.*(-n12x)+ds2zydF2zzfg.*(-n12y)+ds2zzdF2zzfg.*(-n12z)).*N21zf);
end

KU1u2(n1ef1,n2f1)=KU1u2(n1ef1,n2f1)-gamma/h*M12f;
KU1u2(n1ef2,n2f2)=KU1u2(n1ef2,n2f2)-gamma/h*M12f;
if nsd==3
  KU1u2(n1ef3,n2f3)=KU1u2(n1ef3,n2f3)-gamma/h*M12f;
end

% Compute elemental contributions to lhs
LhsCoup=KU1u2;

end

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  iElem,Nodes,SolutionGlobal,SolutionLocal,Parameters,Faces,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
uD=Parameters.Displacement;
Xe=Nodes';
Xem=sum(Xe(1:NumElementFaces,:),1)/NumElementFaces;
t=Time.Time;

% Get solution
Qe=reshape(SolutionLocal(:,1:msd),[],1);
ue=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
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
  Qxxe=Qe(nle1);
  Qyye=Qe(nle2);
  Qxye=Qe(nle3);
elseif nsd==3
  Qxxe=Qe(nle1);
  Qyye=Qe(nle2);
  Qzze=Qe(nle3);
  Qxye=Qe(nle4);
  Qxze=Qe(nle5);
  Qyze=Qe(nle6);
end
uxe=ue(nle1);
uye=ue(nle2);
if nsd==3
    uze=ue(nle3);
end

% Compute variables at Gauss points
if nsd==2
  Qxxeg=Nle*Qxxe;
  Qyyeg=Nle*Qyye;
  Qxyeg=Nle*Qxye;
elseif nsd==3
  Qxxeg=Nle*Qxxe;
  Qyyeg=Nle*Qyye;
  Qzzeg=Nle*Qzze;
  Qxyeg=Nle*Qxye;
  Qxzeg=Nle*Qxze;
  Qyzeg=Nle*Qyze;
end
uxeg=Nle*uxe;
uyeg=Nle*uye;
if nsd==3
  uzeg=Nle*uze;
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

% Compute stress and its linearization
[~,~,D]=computeStress('no','yes','no',Parameters,Xem,[],Sizes);
D_12=D^(1/2);
Voigt1=D_12(1,1);
Voigt2=D_12(2,1);
Voigt3=D_12(msd,msd);

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
  fp(ne1,1)=NwexT*(Qxxeg)+NweyT*(Qxyeg);
  fp(ne2,1)=NwexT*(Qxyeg)+NweyT*(Qyyeg);
elseif nsd==3
  fp(ne1,1)=NwexT*(Qxxeg)+NweyT*(Qxyeg)+NwezT*(Qxzeg);
  fp(ne2,1)=NwexT*(Qxyeg)+NweyT*(Qyyeg)+NwezT*(Qyzeg);
  fp(ne3,1)=NwexT*(Qxzeg)+NweyT*(Qyzeg)+NwezT*(Qzzeg);
end

ft(1,1)=Nw1eT*(uxeg);
ft(2,1)=Nw1eT*(uyeg);
if nsd==3
  ft(3,1)=Nw1eT*(uzeg);
end
    
% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions('Face',RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  
  % Check boundary
  isDirichlet=max((Faces.Dirichlet(:,1)==iElem) &...
                  (Faces.Dirichlet(:,2)==iFace));
  
  % Indices
  nlef1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
  nlef2=nlef1+NumFaceNodes;
  nlef3=nlef2+NumFaceNodes;
  
  % Flip face
  [FlipFace,FaceRow]=max((Faces.Interior(:,3)==iElem) &...
                         (Faces.Interior(:,4)==iFace));
  if FlipFace
    Node2Match1stNode1=Faces.Interior(FaceRow,5);
    order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
    nlef1=nlef1(order);
    nlef2=nlef2(order);
    nlef3=nlef3(order);
  end
  
  % Compute variables at nodes
  Uxf=Ue(nlef1);
  Uyf=Ue(nlef2);
  if nsd==3
    Uzf=Ue(nlef3);
  end
  
  % Compute variables at Gauss points
  Xfg=Nlf*Xf;
  Uxfg=Nlf*Uxf;
  Uyfg=Nlf*Uyf;
  if nsd==3
    Uzfg=Nlf*Uzf;
  end
  if isDirichlet
    uDfg=uD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    uDxfg=uDfg(:,1);
    uDyfg=uDfg(:,2);
    if nsd==3
      uDzfg=uDfg(:,3);
    end
  end
  
  % Compute basic matrices
  Nw1fT=(wfg.*N1f)';
  
  % Compute rhs
  if isDirichlet
    if nsd==2
      fr(1)=fr(1)+Nw1fT*(-uDxfg.*ny+uDyfg.*nx);
    elseif nsd==3
      fr(1)=fr(1)+Nw1fT*(-uDyfg.*nz+uDzfg.*ny);
      fr(2)=fr(2)+Nw1fT*(+uDxfg.*nz-uDzfg.*nx);
      fr(3)=fr(3)+Nw1fT*(-uDxfg.*ny+uDyfg.*nx);
    end
  end
  
  if not(isDirichlet)
    if nsd==2
      fr(1)=fr(1)+Nw1fT*(-Uxfg.*ny+Uyfg.*nx);
    elseif nsd==3
      fr(1)=fr(1)+Nw1fT*(-Uyfg.*nz+Uzfg.*ny);
      fr(2)=fr(2)+Nw1fT*(+Uxfg.*nz-Uzfg.*nx);
      fr(3)=fr(3)+Nw1fT*(-Uxfg.*ny+Uyfg.*nx);
    end
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