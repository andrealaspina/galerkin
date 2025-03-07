classdef Darcy2PhaseRichards_CG < Formulation  
  
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
ReduceTo1GaussPoint=strcmp(Parameters.ReduceTo1GaussPoint,'yes');
PrintGaussData=strcmp(Parameters.PrintGaussData,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
gamma=Parameters.NitschePenalty;
dS=Parameters.RelativeFiniteDifference;
ep=Parameters.Porosity;
r1=Parameters.Density1;
r2=Parameters.Density2;
g=Parameters.Gravity;
PD=Parameters.Pressure;
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

% Retrieve unit cell
UnitCell=Parameters.UnitCell.Value;

% Get solution
Pe=SolutionGlobal(:,1);
Se=SolutionGlobal(:,2);
if isTimeDependent
  Solde=reshape(SolutionOld(:,2,:),[],BDFo);
end

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,~,pinvNe]=mapShapeFunctions('Element',RefElement,RefElement,Xe(1:nsd+1,:),nsd);

% Reduce to only 1 Gauss point
if ReduceTo1GaussPoint
  Ne=mean(Ne,1);
  Nex=mean(Nex,1);
  Ney=mean(Ney,1);
  if nsd==3
    Nez=mean(Nez,1);
  end
  weg=sum(weg,1);
end

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
if isTimeDependent
  Soldeg=Ne*Solde;
end
q1eg=q1(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
q2eg=q2(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);

% Compute variables at nodes
if any(Faces.Outflow)
  DPxe=pinvNe*DPxeg;
  DPye=pinvNe*DPyeg;
  if nsd==3
    DPze=pinvNe*DPzeg;
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

% Compute effective parameters
if nsd==2
  [Geg,K11xxeg,K11yyeg,K12xxeg,K12yyeg,...
       K21xxeg,K21yyeg,K22xxeg,K22yyeg]=...
    solveUnitCell(Seg,iElem,Parameters,UnitCell);
elseif nsd==3
  [Geg,K11xxeg,K11yyeg,K11zzeg,K12xxeg,K12yyeg,K12zzeg,...
       K21xxeg,K21yyeg,K21zzeg,K22xxeg,K22yyeg,K22zzeg]=...
    solveUnitCell(Seg,iElem,Parameters,UnitCell);
end

% Compute derivatives of effective parameters wrt saturation
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
if ReduceTo1GaussPoint && PrintGaussData
  fprintf('\n    Elem = %d',iElem);
  if nsd==2
    fprintf(['\t[S,G,K11xx,K11yy,K12xx,K12yy,',...
                    'K21xx,K21yy,K22xx,K22yy]=',...
             '[%.3f,%.2f,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             Seg,Geg,K11xxeg,K11yyeg,K12xxeg,K12yyeg,...
                     K21xxeg,K21yyeg,K22xxeg,K22yyeg);
    fprintf(['\t[dGdS,dK11xxdS,dK11yydS,dK12xxdS,dK12yydS,',...
                     'dK21xxdS,dK21yydS,dK22xxdS,dK22yydS]=',...
             '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             dGdSeg,dK11xxdSeg,dK11yydSeg,dK12xxdSeg,dK12yydSeg,...
                    dK21xxdSeg,dK21yydSeg,dK22xxdSeg,dK22yydSeg);
  elseif nsd==3
    fprintf(['\t[S,G,K11xx,K11yy,K11zz,K12xx,K12yy,K12zz,',...
                    'K21xx,K21yy,K21zz,K22xx,K22yy,K22zz]=',...
             '[%.3f,%.2f,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             Seg,Geg,K11xxeg,K11yyeg,K11zzeg,K12xxeg,K12yyeg,K12zzeg,...
                     K21xxeg,K21yyeg,K21zzeg,K22xxeg,K22yyeg,K22zzeg);
    fprintf(['\t[dGdS,dK11xxdS,dK11yydS,dK11zzdS,dK12xxdS,dK12yydS,dK12zzdS,',...
                     'dK21xxdS,dK21yydS,dK21zzdS,dK22xxdS,dK22yydS,dK22zzdS]=',...
             '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             dGdSeg,dK11xxdSeg,dK11yydSeg,dK11zzdSeg,dK12xxdSeg,dK12yydSeg,dK12zzdSeg,...
                    dK21xxdSeg,dK21yydSeg,dK21zzdSeg,dK22xxdSeg,dK22yydSeg,dK22zzdSeg);
  end
end

% Compute terms for speedup
gx=g(1);
gy=g(2);
if nsd==3
  gz=g(3);
end
Ge=pinvNe*Geg;
DGxeg=Nex*Ge;
DGyeg=Ney*Ge;
if nsd==3
  DGzeg=Nez*Ge;
end
DGxe=pinvNe*DGxeg;
DGye=pinvNe*DGyeg;
if nsd==3
  DGze=pinvNe*DGzeg;
end
vxeg=(K11xxeg+K21xxeg).*(-DPxeg+DGxeg+r1*gx)+(K12xxeg+K22xxeg).*(-DPxeg-DGxeg+r2*gx);
vyeg=(K11yyeg+K21yyeg).*(-DPyeg+DGyeg+r1*gy)+(K12yyeg+K22yyeg).*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  vzeg=(K11zzeg+K21zzeg).*(-DPzeg+DGzeg+r1*gz)+(K12zzeg+K22zzeg).*(-DPzeg-DGzeg+r2*gz);
end
v2xeg=K21xxeg.*(-DPxeg+DGxeg+r1*gx)+K22xxeg.*(-DPxeg-DGxeg+r2*gx);
v2yeg=K21yyeg.*(-DPyeg+DGyeg+r1*gy)+K22yyeg.*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  v2zeg=K21zzeg.*(-DPzeg+DGzeg+r1*gz)+K22zzeg.*(-DPzeg-DGzeg+r2*gz);
end
dvxdPeg=-(K11xxeg+K21xxeg)-(K12xxeg+K22xxeg);
dvydPeg=-(K11yyeg+K21yyeg)-(K12yyeg+K22yyeg);
if nsd==3
  dvzdPeg=-(K11zzeg+K21zzeg)-(K12zzeg+K22zzeg);
end
dv2xdPeg=-K21xxeg-K22xxeg;
dv2ydPeg=-K21yyeg-K22yyeg;
if nsd==3
  dv2zdPeg=-K21zzeg-K22zzeg;
end
dvxdSKeg=(dK11xxdSeg+dK21xxdSeg).*(-DPxeg+DGxeg+r1*gx)...
        +(dK12xxdSeg+dK22xxdSeg).*(-DPxeg-DGxeg+r2*gx);
dvydSKeg=(dK11yydSeg+dK21yydSeg).*(-DPyeg+DGyeg+r1*gy)...
        +(dK12yydSeg+dK22yydSeg).*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  dvzdSKeg=(dK11zzdSeg+dK21zzdSeg).*(-DPzeg+DGzeg+r1*gz)...
          +(dK12zzdSeg+dK22zzdSeg).*(-DPzeg-DGzeg+r2*gz);
end
dv2xdSKeg=dK21xxdSeg.*(-DPxeg+DGxeg+r1*gx)+dK22xxdSeg.*(-DPxeg-DGxeg+r2*gx);
dv2ydSKeg=dK21yydSeg.*(-DPyeg+DGyeg+r1*gy)+dK22yydSeg.*(-DPyeg-DGyeg+r2*gy);
if nsd==3
  dv2zdSKeg=dK21zzdSeg.*(-DPzeg+DGzeg+r1*gz)+dK22zzdSeg.*(-DPzeg-DGzeg+r2*gz);
end
dvxdSGeg=(K11xxeg+K21xxeg-K12xxeg-K22xxeg).*dGdSeg;
dvydSGeg=(K11yyeg+K21yyeg-K12yyeg-K22yyeg).*dGdSeg;
if nsd==3
  dvzdSGeg=(K11zzeg+K21zzeg-K12zzeg-K22zzeg).*dGdSeg;
end
dv2xdSGeg=(K21xxeg-K22xxeg).*dGdSeg;
dv2ydSGeg=(K21yyeg-K22yyeg).*dGdSeg;
if nsd==3
  dv2zdSGeg=(K21zzeg-K22zzeg).*dGdSeg;
end

% Compute lhs
KPP=-NwexT*(dvxdPeg.*Nex)...
    -NweyT*(dvydPeg.*Ney);
if nsd==3
  KPP=KPP-NwezT*(dvzdPeg.*Nez);
end

KPS=-NwexT*(dvxdSKeg.*Ne+dvxdSGeg.*Nex)...
    -NweyT*(dvydSKeg.*Ne+dvydSGeg.*Ney);
if nsd==3
  KPS=KPS-NwezT*(dvzdSKeg.*Ne+dvzdSGeg.*Nez);
end

KSP=-NwexT*(dv2xdPeg.*Nex)...
    -NweyT*(dv2ydPeg.*Ney);
if nsd==3
  KSP=KSP-NwezT*(dv2zdPeg.*Nez);
end

KSS=-NwexT*(dv2xdSKeg.*Ne+dv2xdSGeg.*Nex)...
    -NweyT*(dv2ydSKeg.*Ne+dv2ydSGeg.*Ney);
if nsd==3
  KSS=KSS-NwezT*(dv2zdSKeg.*Ne+dv2zdSGeg.*Nez);
end

if isTimeDependent
  KSS=KSS+ep*alpha(1)/dt*Me;
end

% Compute rhs
fP=+NwexT*(vxeg)...
   +NweyT*(vyeg);
if nsd==3
  fP=fP+NwezT*(vzeg);
end

fP=fP+NweT*(q1eg/r1+q2eg/r2);

fS=+NwexT*(v2xeg)...
   +NweyT*(v2yeg);
if nsd==3
  fS=fS+NwezT*(v2zeg);
end

if isTimeDependent
  fS=fS-NweT*(ep*[Seg,Soldeg]*alpha/dt);
end

fS=fS+NweT*(q2eg/r2);

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=Faces.Neumann(iFace) || Faces.Outflow(iFace);
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement.FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions('Face',RefElement,RefElement,Xf(1:nsd,:),nsd);
    
    % Reduce to only 1 Gauss point
    if ReduceTo1GaussPoint
      Nf=mean(Nf,1);
      wfg=sum(wfg,1);
    end

    % Compute characteristic element size
    h=sum(wfg);
    
    % Check boundary
    isNeumann=Faces.Neumann(iFace);
    isOutflow=Faces.Outflow(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);

    % Compute derivatives of shape functions
    if isOutflow
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
    Pf=Pe(nf1);
    Sf=Se(nf1);
    if isOutflow
      DPxf=DPxe(nf1);
      DPyf=DPye(nf1);
      if nsd==3
        DPzf=DPze(nf1);
      end
      DGxf=DGxe(nf1);
      DGyf=DGye(nf1);
      if nsd==3
        DGzf=DGze(nf1);
      end
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    Pfg=Nf*Pf;
    Sfg=Nf*Sf;
    if isOutflow
      PDfg=PD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isNeumann
      f1Nfg=f1N(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      f2Nfg=f2N(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isOutflow
      DPxfg=Nf*DPxf;
      DPyfg=Nf*DPyf;
      if nsd==3
        DPzfg=Nf*DPzf;
      end
      DGxfg=Nf*DGxf;
      DGyfg=Nf*DGyf;
      if nsd==3
        DGzfg=Nf*DGzf;
      end
    end
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    Mf=NwfT*Nf;
    Mf=(Mf+Mf')/2;
    
    % Compute effective parameters
    if isOutflow
      if nsd==2
        [Gfg,~,~,~,~,K21xxfg,K21yyfg,K22xxfg,K22yyfg]=...
          solveUnitCell(Sfg,iElem,Parameters,UnitCell);
      elseif nsd==3
        [Gfg,~,~,~,~,~,~,K21xxfg,K21yyfg,K21zzfg,K22xxfg,K22yyfg,K22zzfg]=...
          solveUnitCell(Sfg,iElem,Parameters,UnitCell);
      end
    end
    
    % Compute derivatives of effective parameters wrt saturation
    if isOutflow
      if nsd==2
        [GdSfg,~,~,~,~,K21xxdSfg,K21yydSfg,K22xxdSfg,K22yydSfg]=...
          solveUnitCell(Sfg*(1+dS)+(Sfg==0)*dS,iElem,Parameters,UnitCell);
        dGdSfg=(GdSfg-Gfg)./(Sfg*dS+(Sfg==0)*dS);
        dK21xxdSfg=(K21xxdSfg-K21xxfg)./(Sfg*dS+(Sfg==0)*dS);
        dK21yydSfg=(K21yydSfg-K21yyfg)./(Sfg*dS+(Sfg==0)*dS);
        dK22xxdSfg=(K22xxdSfg-K22xxfg)./(Sfg*dS+(Sfg==0)*dS);
        dK22yydSfg=(K22yydSfg-K22yyfg)./(Sfg*dS+(Sfg==0)*dS);
      elseif nsd==3
        [GdSfg,~,~,~,~,~,~,K21xxdSfg,K21yydSfg,K21zzdSfg,K22xxdSfg,K22yydSfg,K22zzdSfg]=...
          solveUnitCell(Sfg*(1+dS)+(Sfg==0)*dS,iElem,Parameters,UnitCell);
        dGdSfg=(GdSfg-Gfg)./(Sfg*dS+(Sfg==0)*dS);
        dK21xxdSfg=(K21xxdSfg-K21xxfg)./(Sfg*dS+(Sfg==0)*dS);
        dK21yydSfg=(K21yydSfg-K21yyfg)./(Sfg*dS+(Sfg==0)*dS);
        dK21zzdSfg=(K21zzdSfg-K21zzfg)./(Sfg*dS+(Sfg==0)*dS);
        dK22xxdSfg=(K22xxdSfg-K22xxfg)./(Sfg*dS+(Sfg==0)*dS);
        dK22yydSfg=(K22yydSfg-K22yyfg)./(Sfg*dS+(Sfg==0)*dS);
        dK22zzdSfg=(K22zzdSfg-K22zzfg)./(Sfg*dS+(Sfg==0)*dS);
      end
    end

    % Compute terms for speedup
    if isOutflow
      v2xfg=K21xxfg.*(-DPxfg+DGxfg+r1*gx)+K22xxfg.*(-DPxfg-DGxfg+r2*gx);
      v2yfg=K21yyfg.*(-DPyfg+DGyfg+r1*gy)+K22yyfg.*(-DPyfg-DGyfg+r2*gy);
      if nsd==3
        v2zfg=K21zzfg.*(-DPzfg+DGzfg+r1*gz)+K22zzfg.*(-DPzfg-DGzfg+r2*gz);
      end
      dv2xdPfg=-K21xxfg-K22xxfg;
      dv2ydPfg=-K21yyfg-K22yyfg;
      if nsd==3
        dv2zdPfg=-K21zzfg-K22zzfg;
      end
      dv2xdSKfg=dK21xxdSfg.*(-DPxfg+DGxfg+r1*gx)+dK22xxdSfg.*(-DPxfg-DGxfg+r2*gx);
      dv2ydSKfg=dK21yydSfg.*(-DPyfg+DGyfg+r1*gy)+dK22yydSfg.*(-DPyfg-DGyfg+r2*gy);
      if nsd==3
        dv2zdSKfg=dK21zzdSfg.*(-DPzfg+DGzfg+r1*gz)+dK22zzdSfg.*(-DPzfg-DGzfg+r2*gz);
      end
      dv2xdSGfg=(K21xxfg-K22xxfg).*dGdSfg;
      dv2ydSGfg=(K21yyfg-K22yyfg).*dGdSfg;
      if nsd==3
        dv2zdSGfg=(K21zzfg-K22zzfg).*dGdSfg;
      end
      if nsd==2
        v2nfg=v2xfg.*nx+v2yfg.*ny;
      elseif nsd==3
        v2nfg=v2xfg.*nx+v2yfg.*ny+v2zfg.*nz;
      end
    end

    % Compute lhs
    if isOutflow
      KPP(nf1,nf1)=KPP(nf1,nf1)+gamma/h*Mf;

      KSP(nf1,ne1)=KSP(nf1,ne1)+NwfT*((dv2xdPfg.*nx).*Nxf...
                                     +(dv2ydPfg.*ny).*Nyf);
      if nsd==3
        KSP(nf1,ne1)=KSP(nf1,ne1)+NwfT*((dv2zdPfg.*nz).*Nzf);
      end

      KSS(nf1,nf1)=KSS(nf1,nf1)+NwfT*((dv2xdSKfg.*nx).*Nf...
                                     +(dv2ydSKfg.*ny).*Nf);
      if nsd==3
        KSS(nf1,nf1)=KSS(nf1,nf1)+NwfT*((dv2zdSKfg.*nz).*Nf);
      end

      KSS(nf1,ne1)=KSS(nf1,ne1)+NwfT*((dv2xdSGfg.*nx).*Nxf...
                                     +(dv2ydSGfg.*ny).*Nyf);
      if nsd==3
        KSS(nf1,ne1)=KSS(nf1,ne1)+NwfT*((dv2zdSGfg.*nz).*Nzf);
      end
    end
    
    % Compute rhs
    if isNeumann
      fP(nf1,1)=fP(nf1,1)+NwfT*(f1Nfg/r1+f2Nfg/r2);

      fS(nf1,1)=fS(nf1,1)+NwfT*(f2Nfg/r2);
    end

    if isOutflow
      fP(nf1,1)=fP(nf1,1)+NwfT*(-gamma/h*(Pfg-PDfg));
      
      fS(nf1,1)=fS(nf1,1)-NwfT*(v2nfg);
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