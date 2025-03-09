classdef Darcy2Phase_CG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1+1;

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
      Block(iD,iD).SolutionGlobal=zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,1);
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionOld=zeros(Sizes(iD).NumGlobalNodes,Sizes(iD).NumGlobalComp,...
          Parameters(iD).TimeDerOrder);
      end
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,~)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionGlobal=[...
          Parameters(iD).Pressure(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).Saturation(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionGlobal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
              Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              Parameters(iD).Pressure(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).Saturation(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end

    %% Evaluate solution at fixed DOFs
    function [Block]=evaluateSolutionFixedDofs(~,iD,Block,Parameters,Mesh,Time,Sizes)
      Block(iD,iD).SolutionGlobal(Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp),:)=[...
        Parameters(iD).Pressure(...
        Mesh(iD).Nodes(1,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(2,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(3,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',Time.Time),...
        Parameters(iD).Saturation(...
        Mesh(iD).Nodes(1,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(2,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',...
        Mesh(iD).Nodes(3,Block(iD,iD).DOFsFixed(1:end/Sizes(iD).NumGlobalComp))',Time.Time)];
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,~,Parameters,~,~,Time,RefElement,...
        Sizes)
      NodesElem=Elements(iD1).Nodes;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem]=...
          buildBlockElement(iElem,NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionOldElem{iElem},...
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
        Results(iD).Saturation=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).Pressure(:,:,iST)=Block(iD,iD).SolutionGlobal(:,1);
      Results(iD).Saturation(:,:,iST)=Block(iD,iD).SolutionGlobal(:,2);
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal]=buildBlockElement(...
  iElem,Nodes,Faces,SolutionGlobal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
dS=Parameters.RelativeFiniteDifference;
ep=Parameters.Porosity;
r1=Parameters.Density1;
r2=Parameters.Density2;
g=Parameters.Gravity;
q1=Parameters.MassFlowRate1;
q2=Parameters.MassFlowRate2;
f1N=Parameters.NormalMomentum1;
f2N=Parameters.NormalMomentum2;
solveUnitCell=Parameters.solveUnitCell;
Xe=Nodes';
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
Pe=SolutionGlobal(:,1);
Se=SolutionGlobal(:,2);
if isTimeDependent
  Solde=reshape(SolutionOld(:,2,:),[],BDFo);
end

% Initialize lhs
KPP=zeros(NumElementNodes,NumElementNodes);
KPS=zeros(NumElementNodes,NumElementNodes);
KSP=zeros(NumElementNodes,NumElementNodes);
KSS=zeros(NumElementNodes,NumElementNodes);

% Initialize rhs
fP=zeros(NumElementNodes,1);
fS=zeros(NumElementNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe(1:nsd+1,:),nsd);

% Reduce to only 1 Gauss point
Ne=mean(Ne);
Nex=mean(Nex);
Ney=mean(Ney);
if nsd==3
  Nez=mean(Nez);
end
weg=sum(weg);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;

% Compute variables at Gauss points
Xeg=Ne*Xe;
DPxeg=Nex*Pe;
DPyeg=Ney*Pe;
if nsd==3
  DPzeg=Nez*Pe;
end
Seg=Ne*Se;
DSxeg=Nex*Se;
DSyeg=Ney*Se;
if nsd==3
  DSzeg=Nez*Se;
end
if isTimeDependent
  Soldeg=Ne*Solde;
end
q1eg=q1(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
q2eg=q2(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);

% Compute basic matrices
NweT=(weg.*Ne)';
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
Me=NweT*Ne;
Me=(Me+Me')/2;

% Retrieve unit cell
UnitCell=Parameters.UnitCell.Value;

% Solve unit cell
if nsd==2
  [Geg,K11xxeg,K11yyeg,K12xxeg,K12yyeg,...
       K21xxeg,K21yyeg,K22xxeg,K22yyeg]=...
    solveUnitCell(Seg,iElem,Parameters,UnitCell);
elseif nsd==3
  [Geg,K11xxeg,K11yyeg,K11zzeg,K12xxeg,K12yyeg,K12zzeg,...
       K21xxeg,K21yyeg,K21zzeg,K22xxeg,K22yyeg,K22zzeg]=...
    solveUnitCell(Seg,iElem,Parameters,UnitCell);
end

% Compute derivatives wrt saturation
if nsd==2
  [GdSeg,K11xxdSeg,K11yydSeg,K12xxdSeg,K12yydSeg,...
         K21xxdSeg,K21yydSeg,K22xxdSeg,K22yydSeg]=...
    solveUnitCell(Seg*(1+dS)+(Seg==0)*dS,iElem,Parameters,UnitCell);
  dGdSeg=(GdSeg-Geg)./(Seg*dS+(Seg==0)*dS);
  dK11xxdSeg=(K11xxdSeg-K11xxeg)./(Seg*dS+(Seg==0)*dS);
  dK11yydSeg=(K11yydSeg-K11yyeg)./(Seg*dS+(Seg==0)*dS);
  dK12xxdSeg=(K12xxdSeg-K12xxeg)./(Seg*dS+(Seg==0)*dS);
  dK12yydSeg=(K12yydSeg-K12yyeg)./(Seg*dS+(Seg==0)*dS);
  dK21xxdSeg=(K21xxdSeg-K21xxeg)./(Seg*dS+(Seg==0)*dS);
  dK21yydSeg=(K21yydSeg-K21yyeg)./(Seg*dS+(Seg==0)*dS);
  dK22xxdSeg=(K22xxdSeg-K22xxeg)./(Seg*dS+(Seg==0)*dS);
  dK22yydSeg=(K22yydSeg-K22yyeg)./(Seg*dS+(Seg==0)*dS);
elseif nsd==3
  [GdSeg,K11xxdSeg,K11yydSeg,K11zzdSeg,K12xxdSeg,K12yydSeg,K12zzdSeg,...
         K21xxdSeg,K21yydSeg,K21zzdSeg,K22xxdSeg,K22yydSeg,K22zzdSeg]=...
    solveUnitCell(Seg*(1+dS)+(Seg==0)*dS,iElem,Parameters,UnitCell);
  dGdSeg=(GdSeg-Geg)./(Seg*dS+(Seg==0)*dS);
  dK11xxdSeg=(K11xxdSeg-K11xxeg)./(Seg*dS+(Seg==0)*dS);
  dK11yydSeg=(K11yydSeg-K11yyeg)./(Seg*dS+(Seg==0)*dS);
  dK11zzdSeg=(K11zzdSeg-K11zzeg)./(Seg*dS+(Seg==0)*dS);
  dK12xxdSeg=(K12xxdSeg-K12xxeg)./(Seg*dS+(Seg==0)*dS);
  dK12yydSeg=(K12yydSeg-K12yyeg)./(Seg*dS+(Seg==0)*dS);
  dK12zzdSeg=(K12zzdSeg-K12zzeg)./(Seg*dS+(Seg==0)*dS);
  dK21xxdSeg=(K21xxdSeg-K21xxeg)./(Seg*dS+(Seg==0)*dS);
  dK21yydSeg=(K21yydSeg-K21yyeg)./(Seg*dS+(Seg==0)*dS);
  dK21zzdSeg=(K21zzdSeg-K21zzeg)./(Seg*dS+(Seg==0)*dS);
  dK22xxdSeg=(K22xxdSeg-K22xxeg)./(Seg*dS+(Seg==0)*dS);
  dK22yydSeg=(K22yydSeg-K22yyeg)./(Seg*dS+(Seg==0)*dS);
  dK22zzdSeg=(K22zzdSeg-K22zzeg)./(Seg*dS+(Seg==0)*dS);
end

% Print effective parameters and derivatives
% fprintf('\nElem = %d',iElem);
% if nsd==2
%   fprintf(['\t[G,K11xx,K11yy,K12xx,K12yy,',...
%                 'K21xx,K21yy,K22xx,K22yy]=',...
%            '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
%            Geg,K11xxeg,K11yyeg,K12xxeg,K12yyeg,...
%                K21xxeg,K21yyeg,K22xxeg,K22yyeg);
%   fprintf(['\t[dGdS,dK11xxdS,dK11yydS,dK12xxdS,dK12yydS,',...
%                    'dK21xxdS,dK21yydS,dK22xxdS,dK22yydS]=',...
%            '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
%            dGdSeg,dK11xxdSeg,dK11yydSeg,dK12xxdSeg,dK12yydSeg,...
%                   dK21xxdSeg,dK21yydSeg,dK22xxdSeg,dK22yydSeg);
% elseif nsd==3
%   fprintf(['\t[G,K11xx,K11yy,K11zz,K12xx,K12yy,K12zz,',...
%                 'K21xx,K21yy,K21zz,K22xx,K22yy,K22zz]=',...
%            '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
%            Geg,K11xxeg,K11yyeg,K11zzeg,K12xxeg,K12yyeg,K12zzeg,...
%                K21xxeg,K21yyeg,K21zzeg,K22xxeg,K22yyeg,K22zzeg);
%   fprintf(['\t[dGdS,dK11xxdS,dK11yydS,dK11zzdS,dK12xxdS,dK12yydS,dK12zzdS,',...
%                    'dK21xxdS,dK21yydS,dK21zzdS,dK22xxdS,dK22yydS,dK22zzdS]=',...
%            '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
%            dGdSeg,dK11xxdSeg,dK11yydSeg,dK11zzdSeg,dK12xxdSeg,dK12yydSeg,dK12zzdSeg,...
%                   dK21xxdSeg,dK21yydSeg,dK21zzdSeg,dK22xxdSeg,dK22yydSeg,dK22zzdSeg);
% end

% Compute terms for speedup
DGxeg=dGdSeg.*DSxeg; % MAYBE I NEED TO LINEARIZE dG/dS BY COMPUTING d^2G/dS^2!
DGyeg=dGdSeg.*DSyeg;
if nsd==3
  DGzeg=Nez*dGdSeg.*DSzeg;
end
gx=g(1);
gy=g(2);
if nsd==3
  gz=g(3);
end
v1xeg=K11xxeg.*(-DPxeg+DGxeg+r1*gx)+K12xxeg.*(-DPxeg-DGxeg+r2*gx);
v1yeg=K11yyeg.*(-DPyeg+DGyeg+r1*gy)+K12yyeg.*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  v1zeg=K11zzeg.*(-DPzeg+DGzeg+r1*gz)+K12zzeg.*(-DPzeg-DGzeg+r2*gz);
end
v2xeg=K21xxeg.*(-DPxeg+DGxeg+r1*gx)+K22xxeg.*(-DPxeg-DGxeg+r2*gx);
v2yeg=K21yyeg.*(-DPyeg+DGyeg+r1*gy)+K22yyeg.*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  v2zeg=K21zzeg.*(-DPzeg+DGzeg+r1*gz)+K22zzeg.*(-DPzeg-DGzeg+r2*gz);
end
dv1xdPeg=-K11xxeg-K12xxeg;
dv1ydPeg=-K11yyeg-K12yyeg;
if nsd==3
  dv1zdPeg=-K11zzeg-K12zzeg;
end
dv2xdPeg=-K21xxeg-K22xxeg;
dv2ydPeg=-K21yyeg-K22yyeg;
if nsd==3
  dv2zdPeg=-K21zzeg-K22zzeg;
end
dv1xdSKeg=dK11xxdSeg.*(-DPxeg+DGxeg+r1*gx)+dK12xxdSeg.*(-DPxeg-DGxeg+r2*gx);
dv1ydSKeg=dK11yydSeg.*(-DPyeg+DGyeg+r1*gy)+dK12yydSeg.*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  dv1zdSKeg=dK11zzdSeg.*(-DPzeg+DGzeg+r1*gz)+dK12zzdSeg.*(-DPzeg-DGzeg+r2*gz);
end
dv2xdSKeg=dK21xxdSeg.*(-DPxeg+DGxeg+r1*gx)+dK22xxdSeg.*(-DPxeg-DGxeg+r2*gx);
dv2ydSKeg=dK21yydSeg.*(-DPyeg+DGyeg+r1*gy)+dK22yydSeg.*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  dv2zdSKeg=dK21zzdSeg.*(-DPzeg+DGzeg+r1*gz)+dK22zzdSeg.*(-DPzeg-DGzeg+r2*gz);
end
dv1xdSGeg=(K11xxeg-K12xxeg).*dGdSeg;
dv1ydSGeg=(K11yyeg-K12yyeg).*dGdSeg;
if nsd==3
  dv1zdSGeg=(K11zzeg-K12zzeg).*dGdSeg;
end
dv2xdSGeg=(K21xxeg-K22xxeg).*dGdSeg;
dv2ydSGeg=(K21yyeg-K22yyeg).*dGdSeg;
if nsd==3
  dv2zdSGeg=(K21zzeg-K22zzeg).*dGdSeg;
end

% Compute lhs
KPP(ne1,ne1)=-NwexT*((r1*dv1xdPeg).*Nex)...
             -NweyT*((r1*dv1ydPeg).*Ney);
if nsd==3
  KPP(ne1,ne1)=KPP(ne1,ne1)-NwezT*((r1*dv1zdPeg).*Nez);
end

if isTimeDependent
  KPS(ne1,ne1)=-ep*r1*alpha(1)/dt*Me;
end

KPS(ne1,ne1)=KPS(ne1,ne1)-NwexT*((r1*dv1xdSKeg).*Ne+(r1*dv1xdSGeg).*Nex)...
                         -NweyT*((r1*dv1ydSKeg).*Ne+(r1*dv1ydSGeg).*Ney);
if nsd==3
  KPS(ne1,ne1)=KPS(ne1,ne1)-NwezT*((r1*dv1zdSKeg).*Ne+(r1*dv1zdSGeg).*Nez);
end

KSP(ne1,ne1)=-NwexT*((r2*dv2xdPeg).*Nex)...
             -NweyT*((r2*dv2ydPeg).*Ney);
if nsd==3
  KSP(ne1,ne1)=KSP(ne1,ne1)-NwezT*((r2*dv2zdPeg).*Nez);
end

if isTimeDependent
  KSS(ne1,ne1)=+ep*r2*alpha(1)/dt*Me;
end

KSS(ne1,ne1)=KSS(ne1,ne1)-NwexT*((r2*dv2xdSKeg).*Ne+(r2*dv2xdSGeg).*Nex)...
                         -NweyT*((r2*dv2ydSKeg).*Ne+(r2*dv2ydSGeg).*Ney);
if nsd==3
  KSS(ne1,ne1)=KSS(ne1,ne1)-NwezT*((r2*dv2zdSKeg).*Ne+(r2*dv2zdSGeg).*Nez);
end

% Compute rhs
if isTimeDependent
  fP(ne1,1)=+NweT*(ep*r1/dt*Seg*alpha(1)...
                  +ep*r1/dt*Soldeg*alpha(2:BDFo+1,1));
end

fP(ne1,1)=fP(ne1,1)+NwexT*(r1*v1xeg)...
                   +NweyT*(r1*v1yeg);
if nsd==3
  fP(ne1,1)=fP(ne1,1)+NwezT*(r1*v1zeg);
end

fP(ne1,1)=fP(ne1,1)+NweT*(q1eg);

if isTimeDependent
  fS(ne1,1)=-NweT*(ep*r2/dt*Seg*alpha(1)...
                  +ep*r2/dt*Soldeg*alpha(2:BDFo+1,1));
end

fS(ne1,1)=fS(ne1,1)+NwexT*(r2*v2xeg)...
                   +NweyT*(r2*v2yeg);
if nsd==3
  fS(ne1,1)=fS(ne1,1)+NwezT*(r2*v2zeg);
end

fS(ne1,1)=fS(ne1,1)+NweT*(q2eg);

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=Faces.Neumann(iFace);
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement.FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,~,~,~,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf(1:nsd,:),nsd);
    
    % Reduce to only 1 Gauss point
    Nf=mean(Nf);
    wfg=sum(wfg);
    
    % Check boundary
    isNeumann=Faces.Neumann(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    if isNeumann
      f1Nfg=f1N(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      f2Nfg=f2N(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    
    % Compute rhs
    if isNeumann
      fP(nf1,1)=fP(nf1,1)+NwfT*(f1Nfg);

      fS(nf1,1)=fS(nf1,1)+NwfT*(f2Nfg);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Initialization of lhs and rhs
LhsGlobal=zeros(2*NumElementNodes,2*NumElementNodes);
RhsGlobal=zeros(2*NumElementNodes,1);

% Lhs for global problem
LhsGlobal(ne1,ne1)=KPP;
LhsGlobal(ne1,ne2)=KPS;
LhsGlobal(ne2,ne1)=KSP;
LhsGlobal(ne2,ne2)=KSS;

% Rhs for global problem
RhsGlobal(ne1,1)=fP;
RhsGlobal(ne2,1)=fS;

end