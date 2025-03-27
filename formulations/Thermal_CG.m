classdef Thermal_CG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1;

    % Discretization type
    DiscretizationType='CG';

    % Time derivative order
    TimeDerOrder=1;
    
    % Time/frequency domain
    Domain='Time';
    
  end
  
  methods
    
    %% Initialize unknowns
    function [Block]=initializeUnknowns(~,iD,Block,Parameters,Time,Sizes)
      Block(iD,iD).SolutionGlobal=MP*zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,1);
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionOld=MP*zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,...
          Parameters(iD).TimeDerOrder);
      end
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,~)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionGlobal=Parameters(iD).Temperature(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime);
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionGlobal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
              Block(iD,iD).SolutionOld(:,:,1+iBDF)=Parameters(iD).Temperature(....
                Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld);
          end
        end
      end
    end

    %% Evaluate solution at fixed DOFs
    function [Block]=evaluateSolutionFixedDofs(~,iD,Block,Parameters,Mesh,Time,Sizes)
      Block(iD,iD).SolutionGlobal(Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp),:)=...
        Parameters(iD).Temperature(...
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
      LhsCoef=MP*zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=MP*zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      
      % Extract coupling data
      SolutionGlobalElemCoupled=MP*double.empty(Sizes.NumElements,0);
      if matchField(Faces(iD1,iD1),'Interface') && not(isempty(Faces(iD1,iD1).Interface))
        iD2=setdiff(1:2,iD1);
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
          Parameters,Time,RefElement.Value,Sizes); %#ok
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
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
    function [Results]=storeResults(~,iD,iST,Results,Block,~,~,~,Time,~,~)
      if iST==1
        Results(iD).Time=MP*[];
        Results(iD).Temperature=MP*[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).Temperature(:,:,iST)=Block(iD,iD).SolutionGlobal;
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal]=buildBlockElement(...
  iD1,Nodes,Faces,SolutionGlobal,SolutionOld,SolutionGlobalCoupled,...
  Parameters,Time,RefElement,Sizes)

% Get general parameters
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes=Sizes(iD1).NumElementNodes;
NumElementFaces=Sizes(iD1).NumElementFaces;
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
gamma=Parameters(iD1).NitschePenalty;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
ue=SolutionGlobal;
if isTimeDependent
  uolde=reshape(SolutionOld,[],BDFo);
end

% Initialize lhs
Kuu=MP*zeros(NumElementNodes,NumElementNodes);

% Initialize rhs
fu=MP*zeros(NumElementNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,~,pinvNe]=mapShapeFunctions(1,RefElement(iD1,iD1),RefElement(iD1,iD1),Xe,nsd);

% Indices
ne1=1:NumElementNodes;

% Compute variables at nodes
Duxe=pinvNe*(Nex*ue);
Duye=pinvNe*(Ney*ue);
if nsd==3
  Duze=pinvNe*(Nez*ue);
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
ueg=Ne*ue;
if isTimeDependent
  uoldeg=Ne*uolde;
end
Duxeg=Ne*Duxe;
Duyeg=Ne*Duye;
if nsd==3
  Duzeg=Ne*Duze;
end
seg=s(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);

% Compute basic matrices
NweT=(weg.*Ne)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
Me=NweT*Ne;
Me=(Me+Me')/2;
Kxxe=NwexT*Nex;
Kxxe=(Kxxe+Kxxe')/2;
Kyye=NweyT*Ney;
Kyye=(Kyye+Kyye')/2;
if nsd==3
  Kzze=NwezT*Nez;
  Kzze=(Kzze+Kzze')/2;
end

% Compute lhs
if isTimeDependent
  Kuu(ne1,ne1)=rho*cp*alpha(1)/dt*Me;
end

Kuu(ne1,ne1)=Kuu(ne1,ne1)+kappa*Kxxe...
                         +kappa*Kyye;
if nsd==3
  Kuu(ne1,ne1)=Kuu(ne1,ne1)+kappa*Kzze;
end

% Compute rhs
if isTimeDependent
  fu(ne1,1)=-NweT*(rho*cp/dt*ueg*alpha(1)...
                  +rho*cp/dt*uoldeg*alpha(2:BDFo+1,1));
end

fu(ne1,1)=fu(ne1,1)-kappa*NwexT*(Duxeg)...
                   -kappa*NweyT*(Duyeg);
if nsd==3
  fu(ne1,1)=fu(ne1,1)-kappa*NwezT*(Duzeg);
end

fu(ne1,1)=fu(ne1,1)+NweT*(seg);

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=Faces.Exterior(iFace);
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement(iD1,iD1).FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement(iD1,iD1),RefElement(iD1,iD1),Xf,nsd);
    
    % Compute characteristic element size
    h=sum(wfg);
    
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
    uf=ue(nf1);
    Duxf=Duxe(nf1);
    Duyf=Duye(nf1);
    if nsd==3
      Duzf=Duze(nf1);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    ufg=Nf*uf;
    if isDirichlet || isInterface
      Duxfg=Nf*Duxf;
      Duyfg=Nf*Duyf;
      if nsd==3
        Duzfg=Nf*Duzf;
      end
    end
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
    if isDirichlet || isInterface
      NwxfT=(wfg.*Nxf)';
      NwyfT=(wfg.*Nyf)';
      if nsd==3
        NwzfT=(wfg.*Nzf)';
      end
      Cxnxf=NwxfT*(nx.*Nf);
      Cynyf=NwyfT*(ny.*Nf);
      if nsd==3
        Cznzf=NwzfT*(nz.*Nf);
      end
    end
    
    % Get quantities for coupling ------------------------------------------------------------------
    if isInterface
      % Get general parameters
      iD2=setdiff(1:2,iD1);
      iFace2=Faces.Interface(2,iFace);
      NumFaceNodes2=Sizes(iD2).NumFaceNodes;
      
      % Get solution
      U2e=SolutionGlobalCoupled{iFace};
      
      % Compute weights at Gauss points
      [N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions(0,RefElement(iD1,iD2),...
                                                      RefElement(iD1,iD2),Xf,nsd);
      N21f=RefElement(iD2,iD1).ShapeFunctionsFace;
      
      % Indices
      n2ef1=(iFace2-1)*NumFaceNodes2+(1:NumFaceNodes2);
      
      % Flip face
      Node2Match1stNode1=Faces.Interface(3,iFace);
      order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
      n2ef1=n2ef1(order);
      
      % Compute derivatives of shape functions
      N12xf=N12f*Nxe(nf1,:);
      N12yf=N12f*Nye(nf1,:);
      if nsd==3
        N12zf=N12f*Nze(nf1,:);
      end
      
      % Compute variables at nodes
      U2f=U2e(n2ef1);
      
      % Compute variables at Gauss points
      U2fg=N21f*U2f;
      
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
      Kuu(nf1,ne1)=Kuu(nf1,ne1)-kappa*Cxnxf'...
                               -kappa*Cynyf';
      if nsd==3
        Kuu(nf1,ne1)=Kuu(nf1,ne1)-kappa*Cznzf';
      end
      
      Kuu(nf1,nf1)=Kuu(nf1,nf1)+gamma/h*Mf;
      
      Kuu(ne1,nf1)=Kuu(ne1,nf1)-kappa*Cxnxf...
                               -kappa*Cynyf;
      if nsd==3
        Kuu(ne1,nf1)=Kuu(ne1,nf1)-kappa*Cznzf;
      end
    end

    if isRobin
      Kuu(nf1,nf1)=Kuu(nf1,nf1)+NwfT*(Hfg.*Nf);
    end
    
    % Compute rhs
    if isDirichlet
      fu(nf1,1)=fu(nf1,1)+NwfT*(+kappa*(Duxfg.*nx+Duyfg.*ny)...
                                -gamma/h*(ufg-uDfg));
      if nsd==3
        fu(nf1,1)=fu(nf1,1)+NwfT*(kappa*Duzfg.*nz);
      end
      
      fu(ne1,1)=fu(ne1,1)+kappa*NwxfT*(nx.*(ufg-uDfg))...
                         +kappa*NwyfT*(ny.*(ufg-uDfg));
      if nsd==3
        fu(ne1,1)=fu(ne1,1)+kappa*NwzfT*(nz.*(ufg-uDfg));
      end
    end
    
    if isInterface
      fu(nf1,1)=fu(nf1,1)+NwfT*(+kappa*(Duxfg.*nx+Duyfg.*ny)...
                                -gamma/h*ufg);
      if nsd==3
        fu(nf1,1)=fu(nf1,1)+NwfT*(kappa*Duzfg.*nz);
      end
      
      fu(nf1,1)=fu(nf1,1)+Nw12fT*(gamma/h*U2fg);
      
      fu(ne1,1)=fu(ne1,1)+kappa*NwxfT*(nx.*ufg)...
                         +kappa*NwyfT*(ny.*ufg);
      if nsd==3
        fu(ne1,1)=fu(ne1,1)+kappa*NwzfT*(nz.*ufg);
      end
      
      fu(ne1,1)=fu(ne1,1)-kappa*Nw12xfT*(n12x.*U2fg)...
                         -kappa*Nw12yfT*(n12y.*U2fg);
      if nsd==3
        fu(ne1,1)=fu(ne1,1)-kappa*Nw12zfT*(n12z.*U2fg);
      end
    end
    
    if isNeumann
      fu(nf1,1)=fu(nf1,1)+NwfT*(fNfg);
    end

    if isRobin
      fu(nf1,1)=fu(nf1,1)+NwfT*(Hfg.*(uinffg-ufg));
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Lhs for global problem
LhsGlobal=Kuu;
LhsGlobal=(LhsGlobal+LhsGlobal.')/2;

% Rhs for global problem
RhsGlobal=fu;

end

%% Do coupling element
function [LhsCoup]=doCouplingElement(...
  iFaceInterface,iD1,iD2,Parameters,Mesh,Faces,RefElement,Sizes)

% Get general parameters
iElem1=Faces(iD1,iD2).Interface(iFaceInterface,1);
iFace1=Faces(iD1,iD2).Interface(iFaceInterface,2);
iFace2=Faces(iD1,iD2).Interface(iFaceInterface,4);
nsd=Sizes(iD1).NumSpaceDim;
NumElementNodes1=Sizes(iD1).NumElementNodes;
NumElementFaces2=Sizes(iD2).NumElementFaces;
NumFaceNodes2=Sizes(iD2).NumFaceNodes;
C1e=Mesh(iD1).Elements(:,iElem1)';
X1e=Mesh(iD1).Nodes(:,C1e)';
kappa=Parameters(iD1).ThermalConductivity;
gamma=Parameters(iD1).NitschePenalty;
MP=Parameters(iD1).MP;

% Initialize lhs
Ku1U2=MP*zeros(NumElementNodes1,NumElementFaces2*NumFaceNodes2);

% Compute weights at Gauss points
[~,N1ex,N1ey,N1ez,~,~,pinvN1e]=mapShapeFunctions(1,RefElement(iD1,iD1),RefElement(iD1,iD1),X1e,nsd);
  
% Indices
n1e1=1:NumElementNodes1;

% Compute weights at Gauss points
FaceNodes1=RefElement(iD1,iD1).FaceNodesElem;
X1f=X1e(FaceNodes1(iFace1,:),:);
[~,~,~,~,w1fg]=mapShapeFunctions(0,RefElement(iD1,iD1),RefElement(iD1,iD1),X1f,nsd);
[N12f,n12x,n12y,n12z,w12fg]=mapShapeFunctions(0,RefElement(iD1,iD2),RefElement(iD1,iD2),X1f,nsd);
N21f=RefElement(iD2,iD1).ShapeFunctionsFace;

% Compute characteristic element size
h=sum(w1fg);

% Indices
n1f1=FaceNodes1(iFace1,:);
n2ef1=(iFace2-1)*NumFaceNodes2+(1:NumFaceNodes2);

% Flip face
Node2Match1stNode1=Faces(iD1,iD2).Interface(iFaceInterface,5);
order=flipFace(nsd,Parameters(iD2).Degree,Node2Match1stNode1);
n2ef1=n2ef1(order);

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

% Compute basic matrices
Nw12fT=(w12fg.*N12f)';
M12f=Nw12fT*N21f;
Nw12xfT=(w12fg.*N12xf)';
Nw12yfT=(w12fg.*N12yf)';
if nsd==3
  Nw12zfT=(w12fg.*N12zf)';
end
C12xnxf=Nw12xfT*(n12x.*N21f);
C12ynyf=Nw12yfT*(n12y.*N21f);
if nsd==3
  C12znzf=Nw12zfT*(n12z.*N21f);
end

% Compute lhs
Ku1U2(n1f1,n2ef1)=Ku1U2(n1f1,n2ef1)-gamma/h*M12f;

Ku1U2(n1e1,n2ef1)=Ku1U2(n1e1,n2ef1)+kappa*C12xnxf...
                                   +kappa*C12ynyf;
if nsd==3
  Ku1U2(n1e1,n2ef1)=Ku1U2(n1e1,n2ef1)+kappa*C12znzf;
end

% Compute elemental contributions to lhs
LhsCoup=Ku1U2;

end