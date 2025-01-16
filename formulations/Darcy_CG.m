classdef Darcy_CG < Formulation  
  
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
    function [Block]=initializeUnknowns(~,iD,Block,~,~,Sizes)
      Block(iD,iD).SolutionGlobal=zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,1);
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,~,Block,~,~,~,~,~,~)
      
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,~,Parameters,~,~,Time,RefElement,...
        Sizes)
      NodesElem=Elements(iD1).Nodes;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem]=...
          buildBlockElement(iElem,NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},...
          Parameters,Time,RefElement.Value,Sizes); %#ok
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
      end
      Block(iD1,iD1).LhsGlobal=fsparse(Block(iD1,iD1).LhsRowIndices,...
                                       Block(iD1,iD1).LhsColIndices,LhsCoef(:));
      Block(iD1,iD1).RhsGlobal=fsparse(Block(iD1,iD1).RhsRowIndices,1,RhsCoef(:));
    end
    
    %% Store results
    function [Results]=storeResults(~,iD,iST,Results,Block,~,~,~,Time,~,~)
      if iST==1
        Results(iD).Time=[];
        Results(iD).Pressure=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).Pressure(:,:,iST)=Block(iD,iD).SolutionGlobal;
    end
    
    %% Data for Paraview
    function [PointData,CellData]=dataForParaview(~,Results,~,~,~,~)
      
      % Write pressure
      p=Results.Pressure(:,:,end);
      PointData=[sprintf('\nSCALARS Pressure float\n'),...
                 sprintf('LOOKUP_TABLE default\n'),...
                 sprintf('%.12f\n',p')];
      
      CellData='';
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal]=buildBlockElement(...
  iElem,Nodes,Faces,SolutionGlobal,Parameters,Time,RefElement,Sizes)

% Get general parameters
isNavierStokes=strcmp(Parameters.UnitCellProblem,'Navier-Stokes');
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
gamma=Parameters.NitschePenalty;
dDp=Parameters.RelativeFiniteDifference;
pD=Parameters.Pressure;
fN=Parameters.Flux;
s=Parameters.Source;
KappaLxx=Parameters.PermeabilityLinear(1);
KappaLyy=Parameters.PermeabilityLinear(2);
Kappa=Parameters.Permeability;
Xe=Nodes';
t=Time.Time;

% Get solution
pe=SolutionGlobal;

% Initialize lhs
Kpp=zeros(NumElementNodes,NumElementNodes);

% Initialize rhs
fp=zeros(NumElementNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,~,weg]=mapShapeFunctions('Element',RefElement,RefElement,Xe,nsd);

% Reduce to only 1 Gauss point
Ne=mean(Ne);
Nex=mean(Nex);
Ney=mean(Ney);
weg=sum(weg);

% Indices
ne1=1:NumElementNodes;

% Compute variables at nodes
Dpxe=Ne\(Nex*pe);
Dpye=Ne\(Ney*pe);

% Compute variables at Gauss points
Xeg=Ne*Xe;
Dpxeg=Ne*Dpxe;
Dpyeg=Ne*Dpye;
seg=s(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);

% Compute basic matrices
NweT=(weg.*Ne)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
Kxxe=NwexT*Nex;
Kxxe=(Kxxe+Kxxe')/2;
Kyye=NweyT*Ney;
Kyye=(Kyye+Kyye')/2;

% Initialize permeability and its derivatives
Kappaxx=KappaLxx;
Kappayy=KappaLyy;
dKappaxxdDpx=0;
dKappayydDpx=0;
dKappaxxdDpy=0;
dKappayydDpy=0;

% Compute permeability and its derivatives
if isNavierStokes && not(Dpxeg==0 && Dpyeg==0)

  % Retrieve unit cell
  UnitCell=Parameters.UnitCell.Value;
  
  % Compute permeability
  [Kappaxx,Kappayy]=Kappa(Dpxeg,Dpyeg,Parameters,UnitCell);

  % Compute derivative of permeability wrt pressure gradient in x
  if not(norm(Dpxeg)<norm(Dpyeg)*1e-6 || dDp==0)
    [Kappaxx_dDpx,Kappayy_dDpx]=Kappa(Dpxeg*(1+dDp),Dpyeg,Parameters,UnitCell);
    dKappaxxdDpx=(Kappaxx_dDpx-Kappaxx)/(Dpxeg*dDp);
    dKappayydDpx=(Kappayy_dDpx-Kappayy)/(Dpxeg*dDp);
  else
    dKappaxxdDpx=0;
    dKappayydDpx=0;
  end

  % Compute derivative of permeability wrt pressure gradient in y
  if not(norm(Dpyeg)<norm(Dpxeg)*1e-6 || dDp==0)
    [Kappaxx_dDpy,Kappayy_dDpy]=Kappa(Dpxeg,Dpyeg*(1+dDp),Parameters,UnitCell);
    dKappaxxdDpy=(Kappaxx_dDpy-Kappaxx)/(Dpyeg*dDp);
    dKappayydDpy=(Kappayy_dDpy-Kappayy)/(Dpyeg*dDp);
  else
    dKappaxxdDpy=0;
    dKappayydDpy=0;
  end

end

% Print permeability and its derivatives
fprintf('\nElem = %d',iElem);
fprintf('\t[Kxx,Kyy]=[%.2e,%.2e]',Kappaxx,Kappayy);
fprintf('\t[dKxxdDpx,dKxxdDpy,dKyydDpx,dKyydDpy]=[%.2e,%.2e,%.2e,%.2e]',...
  dKappaxxdDpx,dKappaxxdDpy,dKappayydDpx,dKappayydDpy);

% Compute lhs
Kpp(ne1,ne1)=+Kappaxx*Kxxe...
             +Kappayy*Kyye...
             +NwexT*((dKappaxxdDpx.*Dpxeg).*Nex+(dKappaxxdDpy.*Dpxeg).*Ney)...
             +NweyT*((dKappayydDpx.*Dpyeg).*Nex+(dKappayydDpy.*Dpyeg).*Ney);

% Compute rhs
fp(ne1,1)=-Kappaxx*NwexT*(Dpxeg)...
          -Kappayy*NweyT*(Dpyeg)...
          +NweT*(seg);

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=Faces.Exterior(iFace);
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement.FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,~,wfg]=mapShapeFunctions('Face',RefElement,RefElement,Xf,nsd);
    
    % Reduce to only 1 Gauss point
    Nf=mean(Nf);
    nx=mean(nx);
    ny=mean(ny);
    wfg=sum(wfg);
    
    % Compute characteristic element size
    h=sum(wfg);
    
    % Check boundary
    isDirichlet=Faces.Dirichlet(iFace);
    isNeumann=Faces.Neumann(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    
    % Compute derivatives of shape functions
    if isDirichlet
      Nxe=Ne\Nex;
      Nye=Ne\Ney;
      Nxf=Nf*Nxe(nf1,:);
      Nyf=Nf*Nye(nf1,:);
    end
    
    % Compute variables at nodes
    pf=pe(nf1);
    Dpxf=Dpxe(nf1);
    Dpyf=Dpye(nf1);
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    pfg=Nf*pf;
    if isDirichlet
      Dpxfg=Nf*Dpxf;
      Dpyfg=Nf*Dpyf;
    end
    if isDirichlet
      pDfg=pD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    elseif isNeumann
      fNfg=fN(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    Mf=NwfT*Nf;
    Mf=(Mf+Mf')/2;
    if isDirichlet
      NwxfT=(wfg.*Nxf)';
      NwyfT=(wfg.*Nyf)';
      Cxnxf=NwxfT*(nx.*Nf);
      Cynyf=NwyfT*(ny.*Nf);
    end
    
    % Compute lhs
    if isDirichlet
      Kpp(nf1,ne1)=Kpp(nf1,ne1)-Kappaxx*Cxnxf'...
                               -Kappayy*Cynyf';
      
      Kpp(nf1,nf1)=Kpp(nf1,nf1)+gamma/h*Mf;
      
      Kpp(ne1,nf1)=Kpp(ne1,nf1)-Kappaxx*Cxnxf...
                               -Kappayy*Cynyf;
    end
    
    % Compute rhs
    if isDirichlet
      fp(nf1,1)=fp(nf1,1)+NwfT*(+(Kappaxx.*Dpxfg.*nx+Kappayy.*Dpyfg.*ny)...
                                -gamma/h*(pfg-pDfg));
      
      fp(ne1,1)=fp(ne1,1)+Kappaxx*NwxfT*(nx.*(pfg-pDfg))...
                         +Kappayy*NwyfT*(ny.*(pfg-pDfg));
    end
    
    if isNeumann
      fp(nf1,1)=fp(nf1,1)+NwfT*(fNfg);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Lhs for global problem
LhsGlobal=Kpp;
% LhsGlobal=(LhsGlobal+LhsGlobal.')/2;

% Rhs for global problem
RhsGlobal=fp;

end