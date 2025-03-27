classdef Thermal_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim+1;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) 1;

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
      Block(iD,iD).SolutionGlobal=MP*zeros(Sizes(iD).NumGlobalNodes*Sizes(iD).NumGlobalComp,1);
      Block(iD,iD).SolutionLocal=MP*zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumLocalComp,1);
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionOld=MP*zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumLocalComp,...
          Parameters(iD).TimeDerOrder);
      end
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
          Parameters(iD).Temperature(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
              Parameters(iD).Temperature(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,Simulation,Parameters,~,Faces,Time,...
        RefElement,Sizes)
      NodesElem=Elements(iD1).Nodes;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionLocalElem=Elements(iD1).SolutionLocal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      LhsCoef=MP*zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=MP*zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      MatLocal=cell(Sizes(iD1).NumElements,1);
      VecLocal=cell(Sizes(iD1).NumElements,1);

      % Extract coupling data (coinciding elements)
      SolutionLocalSameElemCoupled=cell(Sizes(iD1).NumElements,1);
      if strcmp(Simulation.Problem,'ProjectionCorrection')
        iD2=2;
        SolutionLocalSameElemCoupled=Elements(iD2).SolutionLocal;
      end
      
      % Extract coupling data
      NodesElemCoupled=MP*double.empty(Sizes.NumElements,0);
      SolutionGlobalElemCoupled=double.empty(Sizes.NumElements,0);
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        iD2=setdiff(1:2,iD1);
        NodesElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        SolutionGlobalElemCoupled=cell(Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces);
        iEF=sub2ind([Sizes(iD1).NumElements,Sizes(iD1).NumElementFaces],...
          Faces(iD1,iD2).Interface(:,1),Faces(iD1,iD2).Interface(:,2));
        NodesElemCoupled(iEF)=Elements(iD2).Nodes(Faces(iD1,iD2).Interface(:,3));
        SolutionGlobalElemCoupled(iEF)=Elements(iD2).SolutionGlobal(Faces(iD1,iD2).Interface(:,3));
      end
      
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem,MatLocalElem,VecLocalElem]=...
          buildBlockElement(iD1,NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          SolutionLocalSameElemCoupled{iElem},...
          NodesElemCoupled(iElem,:),SolutionGlobalElemCoupled(iElem,:),...
          Simulation,Parameters,Time,RefElement.Value,Sizes); %#ok
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
        MatLocal{iElem}=MatLocalElem;
        VecLocal{iElem}=VecLocalElem;
      end
      if Simulation.Digits==16
        Block(iD1,iD1).LhsGlobal=fsparse(Block(iD1,iD1).LhsRowIndices,...
                                         Block(iD1,iD1).LhsColIndices,LhsCoef(:));
        Block(iD1,iD1).RhsGlobal=fsparse(Block(iD1,iD1).RhsRowIndices,1,RhsCoef(:));
      else
        Block(iD1,iD1).LhsGlobal=sparse(Block(iD1,iD1).LhsRowIndices,...
                                        Block(iD1,iD1).LhsColIndices,LhsCoef(:));
        Block(iD1,iD1).RhsGlobal=sparse(Block(iD1,iD1).RhsRowIndices,1,RhsCoef(:));
      end
      Elements(iD1).MatLocal=MatLocal;
      Elements(iD1).VecLocal=VecLocal;
    end
    
    %% Do post-process
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,~,RefElement,Sizes)
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
    
    %% Do coupling
    function [Block]=doCoupling(~,iD1,iD2,Block,~,Simulation,Parameters,Mesh,Faces,~,RefElement,...
        Sizes)
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        LhsCoupCoef=MP*zeros(Sizes(iD1).NumElementLhsCoupCoef(iD2),...
                             Sizes(iD1).NumFacesInterface(iD2));
        for iFaceInterface=1:Sizes(iD1).NumFacesInterface(iD2)
          [LhsCoupElem]=...
            doCouplingElement(iFaceInterface,iD1,iD2,Parameters,Mesh,Faces,RefElement.Value,Sizes);
          LhsCoupCoef(:,iFaceInterface)=reshape(LhsCoupElem',[],1);
        end
        if Simulation.Digits==16
          Block(iD1,iD2).LhsGlobal=fsparse(Block(iD1,iD2).LhsRowIndices,...
                                           Block(iD1,iD2).LhsColIndices,[LhsCoupCoef(:);0]);
        else
          Block(iD1,iD2).LhsGlobal=sparse(Block(iD1,iD2).LhsRowIndices,...
                                          Block(iD1,iD2).LhsColIndices,[LhsCoupCoef(:);0]);
        end
      else
        if Simulation.Digits==16
          Block(iD1,iD2).LhsGlobal=fsparse(Block(iD1,iD2).LhsRowIndices,...
                                           Block(iD1,iD2).LhsColIndices,0);
        else
          Block(iD1,iD2).LhsGlobal=sparse(Block(iD1,iD2).LhsRowIndices,...
                                          Block(iD1,iD2).LhsColIndices,0);
        end
      end
    end
    
    %% Store results
    function [Results]=storeResults(~,iD,iST,Results,Block,~,Parameters,~,Time,~,Sizes)
      if iST==1
        Results(iD).Time=MP*[];
        Results(iD).ScaledTemperatureGradient=MP*[];
        Results(iD).Temperature=MP*[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledTemperatureGradient(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        1:Sizes(iD).NumSpaceDim);
      Results(iD).Temperature(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumSpaceDim+1);
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).TemperaturePost=Block(iD,iD).SolutionPost;
      end
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  iD1,Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,...
  SolutionLocalSameCoupled,...
  NodesCoupled,SolutionGlobalCoupled,...
  Simulation,Parameters,Time,RefElement,Sizes)

% Get general parameters
isTimeDependent=strcmp(Time.TimeDependent,'yes');
isProjectionCorrection=strcmp(Simulation.Problem,'ProjectionCorrection');
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes=Sizes(iD1).NumElementNodes;
NumElementFaces=Sizes(iD1).NumElementFaces;
NumFaceNodes=Sizes(iD1).NumFaceNodes;
rho=Parameters(iD1).Density;
cp=Parameters(iD1).SpecificHeatCapacity;
kappa=Parameters(iD1).ThermalConductivity;
H=Parameters(iD1).ConvectionCoefficient;
uinf=Parameters(iD1).AmbientTemperature;
uD=Parameters(iD1).Temperature;
fN=Parameters(iD1).ThermalFlux;
s=Parameters(iD1).HeatSource;
MP=Parameters(iD1).MP;
Xe=Nodes';
tauU=Parameters(iD1).StabTemperature;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end
if isProjectionCorrection
  isAxisymmetric=strcmp(Parameters(iD1).Axisymmetric,'yes');
  epsilon=Parameters.ElectricPermittivity;
  q=Parameters.ParticleCharge;
  qi=Parameters.IonCharge;
  m=Parameters.ParticleMass;
  mi=Parameters.IonMass;
  ri=Parameters.IonDensity;
end

% Get solution
qe=reshape(SolutionLocal(:,1:nsd),[],1);
ue=reshape(SolutionLocal(:,nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  uolde=reshape(SolutionOld(:,nsd+1,:),[],BDFo);
end
if isProjectionCorrection
  re=reshape(SolutionLocalSameCoupled(:,1),[],1);
  ee=reshape(SolutionLocalSameCoupled(:,1+(1:nsd)),[],1);
end

% Initialize lhs
Kqq=MP*zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kqu=MP*zeros(nsd*NumElementNodes,NumElementNodes);
KqU=MP*zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Kuu=MP*zeros(NumElementNodes,NumElementNodes);
KuU=MP*zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KUU=MP*zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fq=MP*zeros(nsd*NumElementNodes,1);
fu=MP*zeros(NumElementNodes,1);
fU=MP*zeros(NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement(iD1,iD1),RefElement(iD1,iD1),Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
qxe=qe(ne1);
qye=qe(ne2);
if nsd==3
  qze=qe(ne3);
end
if isProjectionCorrection
  exe=ee(ne1);
  eYe=ee(ne2);
  if nsd==3
    eze=ee(ne3);
  end
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
qxeg=Ne*qxe;
qyeg=Ne*qye;
if nsd==3
  qzeg=Ne*qze;
end
ueg=Ne*ue;
if isTimeDependent
  uoldeg=Ne*uolde;
end
seg=s(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
if isProjectionCorrection
  reg=Ne*re;
  if nsd==2
    diveeg=Nex*exe+Ney*eYe;
  else
    diveeg=Nex*exe+Ney*eYe+Nez*eze;
  end
  if isAxisymmetric
    diveeg=diveeg+(Ne*eYe)./Xeg(:,2);
  end
  if isa(epsilon,'function_handle')
    epsiloneg=epsilon(Xeg(:,1),Xeg(:,2),Xeg(:,3));
  else
    epsiloneg=epsilon;
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
Kqq(ne1,ne1)=-Me;
Kqq(ne2,ne2)=-Me;
if nsd==3
  Kqq(ne3,ne3)=-Me;
end

Kqu(ne1,ne1)=+sqrt(kappa)*Cxe;
Kqu(ne2,ne1)=+sqrt(kappa)*Cye;
if nsd==3
  Kqu(ne3,ne1)=+sqrt(kappa)*Cze;
end

if isTimeDependent
  Kuu(ne1,ne1)=rho*cp*alpha(1)/dt*Me;
end

% Compute rhs
fq(ne1,1)=+NweT*(qxeg)...
          -sqrt(kappa)*NwexT*(ueg);
fq(ne2,1)=+NweT*(qyeg)...
          -sqrt(kappa)*NweyT*(ueg);
if nsd==3
  fq(ne3,1)=+NweT*(qzeg)...
            -sqrt(kappa)*NwezT*(ueg);
end

if isTimeDependent
  fu(ne1,1)=-NweT*(rho*cp/dt*ueg*alpha(1)...
                  +rho*cp/dt*uoldeg*alpha(2:BDFo+1,1));
end

fu(ne1,1)=fu(ne1,1)-sqrt(kappa)*NweT*(Nex*qxe+Ney*qye);
if nsd==3
  fu(ne1,1)=fu(ne1,1)-sqrt(kappa)*NweT*(Nez*qze);
end

fu(ne1,1)=fu(ne1,1)+NweT*(seg);

if isProjectionCorrection
  fu(ne1,1)=fu(ne1,1)+NweT*(diveeg-(q/m*reg+qi/mi*ri)./epsiloneg);
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement(iD1,iD1),RefElement(iD1,iD1),Xf,nsd);
    
    % Check boundary
    Boundary=Faces.Boundary(iFace);
    isDirichlet=Faces.Dirichlet(iFace);
    isNeumann=Faces.Neumann(iFace);
    isRobin=Faces.Robin(iFace);
    if matchField(Faces,'Interface')
      isInterface=Faces.Interface(1,iFace);
    else
      isInterface=false;
    end
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nef1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters(iD1).Degree,Node2Match1stNode1);
      nef1=nef1(order);
    end
    
    % Flip face (periodic (slave))
    if matchField(Faces,'Periodic') && not(isempty(Faces.Periodic))
      Node2Match1stNode1=Faces.Periodic(2,iFace);
      FlipFace=max(Node2Match1stNode1);
      if FlipFace
        order=flipFace(nsd,Parameters(iD1).Degree,Node2Match1stNode1);
        nef1=nef1(order);
      end
    end
    
    % Compute variables at nodes
    qxf=qxe(nf1);
    qyf=qye(nf1);
    if nsd==3
      qzf=qze(nf1);
    end
    uf=ue(nf1);
    Uf=Ue(nef1);
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    qxfg=Nf*qxf;
    qyfg=Nf*qyf;
    if nsd==3
      qzfg=Nf*qzf;
    end
    ufg=Nf*uf;
    Ufg=Nf*Uf;
    if isDirichlet
      uDfg=uD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    elseif isNeumann
      fNfg=fN(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    elseif isRobin
      Hfg=H(Xfg(:,1),Xfg(:,2),Xfg(:,3),Boundary);
      uinffg=uinf(Xfg(:,1),Xfg(:,2),Xfg(:,3),Boundary);
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
      X2e=NodesCoupled{iFace}';
      gamma=Parameters(iD2).NitschePenalty;
      
      % Get solution
      u2e=SolutionGlobalCoupled{iFace};
      
      % Compute weights at Gauss points
      [~,N2ex,N2ey,N2ez,~,~,pinvN2e]=mapShapeFunctions(1,RefElement(iD2,iD2),...
                                                         RefElement(iD2,iD2),X2e,nsd);
      
      % Compute variables at nodes
      Du2xe=pinvN2e*(N2ex*u2e);
      Du2ye=pinvN2e*(N2ey*u2e);
      if nsd==3
        Du2ze=pinvN2e*(N2ez*u2e);
      end
      
      % Compute weights at Gauss points
      FaceNodes2=RefElement(iD2,iD2).FaceNodesElem;
      X2f=X2e(FaceNodes2(iFace2,:),:);
      [N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions(0,RefElement(iD1,iD2),...
                                                      RefElement(iD1,iD2),Xf,nsd);
      N21f=RefElement(iD2,iD1).ShapeFunctionsFace;
      [~,~,~,~,w2fg]=mapShapeFunctions(0,RefElement(iD2,iD2),RefElement(iD2,iD2),X2f,nsd);
      
      % Compute characteristic element size
      h=sum(w2fg);
      
      % Indices
      n2f1=FaceNodes2(iFace2,:);
      
      % Flip face
      Node2Match1stNode1=Faces.Interface(3,iFace);
      order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
      n2f1=n2f1(order);
      
      % Compute variables at nodes
      u2f=u2e(n2f1);
      Du2xf=Du2xe(n2f1);
      Du2yf=Du2ye(n2f1);
      if nsd==3
        Du2zf=Du2ze(n2f1);
      end
      
      % Compute variables at Gauss points
      u2fg=N21f*u2f;
      Du2xfg=N21f*Du2xf;
      Du2yfg=N21f*Du2yf;
      if nsd==3
        Du2zfg=N21f*Du2zf;
      end
      
      % Compute basic matrices
      Nw12fT=(w12fg.*N12f)';
    end
    % ----------------------------------------------------------------------------------------------
    
    % Compute lhs
    Kuu(nf1,nf1)=Kuu(nf1,nf1)+tauU*Mf;
    
    if not(isDirichlet)
      KqU(nf1,nef1)=KqU(nf1,nef1)-sqrt(kappa)*Mfnx';
      KqU(nf2,nef1)=KqU(nf2,nef1)-sqrt(kappa)*Mfny';
      if nsd==3
        KqU(nf3,nef1)=KqU(nf3,nef1)-sqrt(kappa)*Mfnz';
      end
      
      KuU(nf1,nef1)=KuU(nf1,nef1)-tauU*Mf;
      
      KUU(nef1,nef1)=KUU(nef1,nef1)+tauU*Mf;
    end
    
    if isInterface
      KUU(nef1,nef1)=KUU(nef1,nef1)+gamma/h*Mf;
    end

    if isRobin
      KUU(nef1,nef1)=KUU(nef1,nef1)+NwfT*(Hfg.*Nf);
    end
    
    % Compute rhs
    fu(nf1,1)=fu(nf1,1)-NwfT*(tauU*ufg);
    
    if isDirichlet
      fq(nf1,1)=fq(nf1,1)+sqrt(kappa)*NwfT*(nx.*uDfg);
      fq(nf2,1)=fq(nf2,1)+sqrt(kappa)*NwfT*(ny.*uDfg);
      if nsd==3
        fq(nf3,1)=fq(nf3,1)+sqrt(kappa)*NwfT*(nz.*uDfg);
      end
      
      fu(nf1,1)=fu(nf1,1)+NwfT*(tauU*uDfg);
    end
    
    if not(isDirichlet)
      fq(nf1,1)=fq(nf1,1)+sqrt(kappa)*NwfT*(nx.*Ufg);
      fq(nf2,1)=fq(nf2,1)+sqrt(kappa)*NwfT*(ny.*Ufg);
      if nsd==3
        fq(nf3,1)=fq(nf3,1)+sqrt(kappa)*NwfT*(nz.*Ufg);
      end
      
      fu(nf1,1)=fu(nf1,1)+NwfT*(tauU*Ufg);
      
      fU(nef1,1)=fU(nef1,1)...
                 +NwfT*(sqrt(kappa)*(qxfg.*nx+qyfg.*ny)...
                        +tauU*(ufg-Ufg));
      if nsd==3
        fU(nef1,1)=fU(nef1,1)+NwfT*(sqrt(kappa)*qzfg.*nz);
      end
    end
    
    if isInterface
      fU(nef1,1)=fU(nef1,1)-Nw12fT*(kappa*(Du2xfg.*(-n12x)+Du2yfg.*(-n12y))...
                                   -gamma/h*u2fg);
      if nsd==3
        fU(nef1,1)=fU(nef1,1)-Nw12fT*(kappa*Du2zfg.*(-n12z));
      end
      
      fU(nef1,1)=fU(nef1,1)-NwfT*(gamma/h*Ufg);
    end
    
    if isNeumann
      fU(nef1,1)=fU(nef1,1)+NwfT*(fNfg);
    end

    if isRobin
      fU(nef1,1)=fU(nef1,1)+NwfT*(Hfg.*(uinffg-ufg));
    end
    
    % Remove undetermination
    if isDirichlet
      KUU(nef1,nef1)=eye(NumFaceNodes);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iq=1:nsd*NumElementNodes;
iu=iq(end)+(1:NumElementNodes);
iU=1:NumElementFaces*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=MP*zeros((nsd+1)*NumElementNodes,(nsd+1)*NumElementNodes);
LhsLG=MP*zeros((nsd+1)*NumElementNodes,NumElementFaces*NumFaceNodes);
LhsGL=MP*zeros(NumElementFaces*NumFaceNodes,(nsd+1)*NumElementNodes);
LhsGG=MP*zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
RhsL=MP*zeros((nsd+1)*NumElementNodes,1);
RhsG=MP*zeros(NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iq,iq)=Kqq;
LhsLL(iq,iu)=Kqu;
LhsLL(iu,iq)=Kqu';
LhsLL(iu,iu)=Kuu;

% Lhs local-global
LhsLG(iq,iU)=KqU;
LhsLG(iu,iU)=KuU;

% Rhs local
RhsL(iq,1)=fq;
RhsL(iu,1)=fu;

% Lhs global-local
LhsGL(iU,iq)=KqU';
LhsGL(iU,iu)=KuU';

% Lhs global-global
LhsGG(iU,iU)=KUU;

% Rhs global
RhsG(iU,1)=fU;

% Add missing terms for axisymmetric formulation for projection correction only
if isProjectionCorrection && isAxisymmetric
  Kuq_axisym=MP*zeros(NumElementNodes,nsd*NumElementNodes);
  Kuq_axisym(ne1,ne2)=+sqrt(kappa)*NweT*((1./Xeg(:,2)).*Ne);
  LhsLL(iu,iq)=LhsLL(iu,iq)+Kuq_axisym;

  fu_axisym=zeros(NumElementNodes,1);
  fu_axisym(ne1,1)=-sqrt(kappa)*NweT*(1./Xeg(:,2).*qyeg);
  RhsL(iu,1)=RhsL(iu,1)+fu_axisym;
end

% Matrix and vector for local problem
MatVecLocal=LhsLL\[LhsLG,RhsL];

% Extract matrix for local problem
MatLocal=MatVecLocal(:,1:end-1);

% Extract vector for local problem
VecLocal=MatVecLocal(:,end);

% Lhs for global problem
LhsGlobal=LhsGG-LhsGL*MatLocal;
if not(isProjectionCorrection && isAxisymmetric)
  LhsGlobal=(LhsGlobal+LhsGlobal.')/2;
end

% Rhs for global problem
RhsGlobal=RhsG-LhsGL*VecLocal;

end

%% Do coupling element
function [LhsCoup]=doCouplingElement(...
  iFaceInterface,iD1,iD2,Parameters,Mesh,Faces,RefElement,Sizes)       

% Get general parameters
iElem1=Faces(iD1,iD2).Interface(iFaceInterface,1);
iElem2=Faces(iD1,iD2).Interface(iFaceInterface,3);
iFace1=Faces(iD1,iD2).Interface(iFaceInterface,2);
iFace2=Faces(iD1,iD2).Interface(iFaceInterface,4);
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes2=Sizes(iD2).NumElementNodes;
NumElementFaces1=Sizes(iD1).NumElementFaces;
NumFaceNodes1=Sizes(iD1).NumFaceNodes;
C1e=Mesh(iD1).Elements(:,iElem1)';
C2e=Mesh(iD2).Elements(:,iElem2)';
X1e=Mesh(iD1).Nodes(:,C1e)';
X2e=Mesh(iD2).Nodes(:,C2e)';
kappa=Parameters(iD1).ThermalConductivity;
gamma=Parameters(iD2).NitschePenalty;
MP=Parameters(iD1).MP;

% Initialize lhs
KU1u2=MP*zeros(NumElementFaces1*NumFaceNodes1,NumElementNodes2);

% Compute weights at Gauss points
[~,N2ex,N2ey,N2ez,~,~,pinvN2e]=mapShapeFunctions(1,RefElement(iD2,iD2),...
                                                   RefElement(iD2,iD2),X2e,nsd);
  
% Indices
n2e1=1:NumElementNodes2;

% Compute weights at Gauss points
FaceNodes1=RefElement(iD1,iD1).FaceNodesElem;
FaceNodes2=RefElement(iD2,iD2).FaceNodesElem;
X1f=X1e(FaceNodes1(iFace1,:),:);
X2f=X2e(FaceNodes2(iFace2,:),:);
[N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions(0,RefElement(iD1,iD2),RefElement(iD1,iD2),X1f,nsd);
N21f=RefElement(iD2,iD1).ShapeFunctionsFace;
[~,~,~,~,w2fg]=mapShapeFunctions(0,RefElement(iD2,iD2),RefElement(iD2,iD2),X2f,nsd);

% Compute characteristic element size
h=sum(w2fg);

% Indices
n1ef1=(iFace1-1)*NumFaceNodes1+(1:NumFaceNodes1);
n2f1=FaceNodes2(iFace2,:);

% Flip face
Node2Match1stNode1=Faces(iD1,iD2).Interface(iFaceInterface,5);
order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
n2f1=n2f1(order);

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

% Compute basic matrices
Nw12fT=(w12fg.*N12f)';
M12f=Nw12fT*N21f;
Nw21xfT=(w12fg.*N21xf)';
Nw21yfT=(w12fg.*N21yf)';
if nsd==3
  Nw21zfT=(w12fg.*N21zf)';
end
C21xnxf=Nw21xfT*(n12x.*N12f);
C21ynyf=Nw21yfT*(n12y.*N12f);
if nsd==3
  C21znzf=Nw21zfT*(n12z.*N12f);
end

% Compute lhs
KU1u2(n1ef1,n2e1)=KU1u2(n1ef1,n2e1)-kappa*C21xnxf'...
                                   -kappa*C21ynyf';
if nsd==3
  KU1u2(n1ef1,n2e1)=KU1u2(n1ef1,n2e1)-kappa*C21znzf';
end

KU1u2(n1ef1,n2f1)=KU1u2(n1ef1,n2f1)-gamma/h*M12f;

% Compute elemental contributions to lhs
LhsCoup=KU1u2;

end

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  Nodes,SolutionLocal,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
kappa=Parameters.ThermalConductivity;
MP=Parameters.MP;
Xe=Nodes';

% Get solution
qe=reshape(SolutionLocal(:,1:nsd),[],1);
ue=reshape(SolutionLocal(:,nsd+1),[],1);

% Initialize lhs
Kpp=MP*zeros(NumElementNodesPost,NumElementNodesPost);
Ktp=MP*zeros(1,NumElementNodesPost);

% Initialize rhs
fp=MP*zeros(NumElementNodesPost,1);
ft=MP*zeros(1,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions(1,RefElement.PostLow,RefElement.Post,Xe,nsd);
N1e=ones(length(weg),1);

% Indices
ne1=1:NumElementNodesPost;
nle1=1:NumElementNodes;
nle2=nle1+NumElementNodes;
nle3=nle2+NumElementNodes;

% Compute variables at nodes
qxe=qe(nle1);
qye=qe(nle2);
if nsd==3
    qze=qe(nle3);
end

% Compute variables at Gauss points
qxeg=Nle*qxe;
qyeg=Nle*qye;
if nsd==3
  qzeg=Nle*qze;
end
ueg=Nle*ue;

% Compute basic matrices
Nw1eT=(weg.*N1e)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
Kxxe=NwexT*Nex;
Kyye=NweyT*Ney;
if nsd==3
  Kzze=NwezT*Nez;
end

% Compute lhs
Kpp(ne1,ne1)=-sqrt(kappa)*Kxxe...
             -sqrt(kappa)*Kyye;
if nsd==3
  Kpp(ne1,ne1)=Kpp(ne1,ne1)-sqrt(kappa)*Kzze;
end

Ktp(1,ne1)=Nw1eT*(Ne);

% Compute rhs
fp(ne1,1)=NwexT*(qxeg)...
         +NweyT*(qyeg);
if nsd==3
  fp(ne1,1)=fp(ne1,1)+NwezT*(qzeg);
end

ft(1,1)=Nw1eT*(ueg);

% Indices
ip=1:NumElementNodesPost;
it=ip(end)+1;

% Initialization of lhs and rhs
LhsPost=MP*zeros(NumElementNodesPost+1,NumElementNodesPost);
RhsPost=MP*zeros(NumElementNodesPost+1,1);

% Lhs for post-processing
LhsPost(ip,ip)=Kpp;
LhsPost(it,ip)=Ktp;

% Rhs for post-processing
RhsPost(ip,1)=fp;
RhsPost(it,1)=ft;

end