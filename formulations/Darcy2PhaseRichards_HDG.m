classdef Darcy2PhaseRichards_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1+1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim+NumSpaceDim+1+1;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) 1+1;

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
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,Faces,Time,RefElement,...
        Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
          zeros(Sizes(iD).NumLocalNodes,1),...
          Parameters(iD).Saturation(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;

        % Impose initial conditions on the saturation trace
        for iFace=1:Sizes(iD).NumFaces
          [iElem,iElemFace]=find(Faces(iD,iD).Connectivity==iFace);
          [iElem,Order]=min(iElem);
          iElemFace=iElemFace(Order);
          FaceNodes=RefElement(iD,iD).FaceNodesElem;
          Ce=Mesh(iD).Elements(:,iElem)';
          Xe=Mesh(iD).Nodes(:,Ce)';
          Xf=Xe(FaceNodes(iElemFace,:),:);
          Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp*Sizes(iD).NumFaceNodes+...
            (1:(1+1)*Sizes(iD).NumFaceNodes))=reshape(...
            [zeros(Sizes(iD).NumFaceNodes,1),...
             Parameters(iD).Saturation(Xf(:,1),Xf(:,2),Xf(:,3),Time.InitialTime)],[],1);
        end
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
              zeros(Sizes(iD).NumLocalNodes,1),...
              Parameters(iD).Saturation(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,~,Parameters,~,~,Time,RefElement,Sizes)
      NodesElem=Elements(iD1).Nodes;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionLocalElem=Elements(iD1).SolutionLocal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      MatLocal=cell(Sizes(iD1).NumElements,1);
      VecLocal=cell(Sizes(iD1).NumElements,1);
      Fill=0;
      Volume=0;
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem,MatLocalElem,VecLocalElem,FillElem,VolumeElem]=...
          buildBlockElement(iElem,NodesElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          Parameters,Time,RefElement.Value,Sizes); %#ok
        LhsCoef(:,iElem)=reshape(LhsGlobalElem',[],1);
        RhsCoef(:,iElem)=reshape(RhsGlobalElem',[],1);
        MatLocal{iElem}=MatLocalElem;
        VecLocal{iElem}=VecLocalElem;
        Fill=Fill+FillElem;
        Volume=Volume+VolumeElem;
      end
      Block(iD1).FillFactor=Fill/Volume;
      Block(iD1,iD1).LhsGlobal=fsparse(Block(iD1,iD1).LhsRowIndices,...
                                       Block(iD1,iD1).LhsColIndices,LhsCoef(:));
      Block(iD1,iD1).RhsGlobal=fsparse(Block(iD1,iD1).RhsRowIndices,1,RhsCoef(:));
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
          doPostProcessElement(iElem,NodesElem{iElem},...
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
        Results(iD).PressureGradient=[];
        Results(iD).ChemicalPotentialGradient=[];
        Results(iD).Pressure=[];
        Results(iD).Saturation=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).PressureGradient(:,:,iST)=Block(iD,iD).SolutionLocal(:,1:Sizes(iD).NumSpaceDim);
      Results(iD).ChemicalPotentialGradient(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumSpaceDim+(1:Sizes(iD).NumSpaceDim));
      Results(iD).Pressure(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumSpaceDim+...
        Sizes(iD).NumSpaceDim+1);
      Results(iD).Saturation(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumSpaceDim+...
        Sizes(iD).NumSpaceDim+1+1);
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).PressurePost=Block(iD,iD).SolutionPost(:,1);
        Results(iD).SaturationPost=Block(iD,iD).SolutionPost(:,1+1);
      end
      Results(iD).FillFactor(iST)=Block(iD,iD).FillFactor;
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal,Fill,Volume]=buildBlockElement(...
  iElem,Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
ReduceTo1GaussPoint=strcmp(Parameters.ReduceTo1GaussPoint,'yes');
PrintGaussData=strcmp(Parameters.PrintGaussData,'yes');
SolvePressureEquation=strcmp(Parameters.SolvePressureEquation,'yes');
EquilibrateLocalProblems=strcmp(Parameters.EquilibrateLocalProblems,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
ds=Parameters.RelativeFiniteDifference;
ep=Parameters.Porosity;
r1=Parameters.Density1;
r2=Parameters.Density2;
g=Parameters.Gravity;
pD=Parameters.Pressure;
sD=Parameters.Saturation;
q1=Parameters.MassFlowRate1;
q2=Parameters.MassFlowRate2;
f1N=Parameters.NormalMomentum1;
f2N=Parameters.NormalMomentum2;
solveUnitCell=Parameters.solveUnitCell;
Xe=Nodes';
tauP=Parameters.StabPressure;
tauS=Parameters.StabSaturation;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Retrieve unit cell
UnitCell=Parameters.UnitCell.Value;

% Get solution
Qe=reshape(SolutionLocal(:,1:nsd),[],1);
He=reshape(SolutionLocal(:,nsd+(1:nsd)),[],1);
pe=reshape(SolutionLocal(:,nsd+nsd+1),[],1);
se=reshape(SolutionLocal(:,nsd+nsd+1+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  solde=reshape(SolutionOld(:,nsd+nsd+1+1,:),[],BDFo);
end

% Initialize lhs
KQQ=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KQp=zeros(nsd*NumElementNodes,NumElementNodes);
KQP=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
KHH=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KHs=zeros(nsd*NumElementNodes,NumElementNodes);
KHS=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
KpQ=zeros(NumElementNodes,nsd*NumElementNodes);
KpH=zeros(NumElementNodes,nsd*NumElementNodes);
Kpp=zeros(NumElementNodes,NumElementNodes);
Kps=zeros(NumElementNodes,NumElementNodes);
KpP=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KpS=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KsQ=zeros(NumElementNodes,nsd*NumElementNodes);
KsH=zeros(NumElementNodes,nsd*NumElementNodes);
Kss=zeros(NumElementNodes,NumElementNodes);
KsS=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KPQ=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KPH=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KPp=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KPP=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KPS=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KSQ=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KSH=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KSs=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KSS=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fQ=zeros(nsd*NumElementNodes,1);
fH=zeros(nsd*NumElementNodes,1);
fp=zeros(NumElementNodes,1);
fs=zeros(NumElementNodes,1);
fP=zeros(NumElementFaces*NumFaceNodes,1);
fS=zeros(NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctionsLinear(1,RefElement,Xe(1:nsd+1,:),nsd);

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
ne3=ne2+NumElementNodes;

% Compute variables at nodes
Qxe=Qe(ne1);
Qye=Qe(ne2);
if nsd==3
  Qze=Qe(ne3);
end
Hxe=He(ne1);
Hye=He(ne2);
if nsd==3
  Hze=He(ne3);
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
Qxeg=Ne*Qxe;
Qyeg=Ne*Qye;
if nsd==3
  Qzeg=Ne*Qze;
end
Hxeg=Ne*Hxe;
Hyeg=Ne*Hye;
if nsd==3
  Hzeg=Ne*Hze;
end
peg=Ne*pe;
seg=Ne*se;
if isTimeDependent
  soldeg=Ne*solde;
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
Cxe=NwexT*Ne;
Cye=NweyT*Ne;
if nsd==3
  Cze=NwezT*Ne;
end

% Compute effective parameters
if nsd==2
  [Geg,K11xxeg,K11yyeg,K12xxeg,K12yyeg,...
       K21xxeg,K21yyeg,K22xxeg,K22yyeg]=...
    solveUnitCell(seg,iElem,Parameters,UnitCell);
elseif nsd==3
  [Geg,K11xxeg,K11yyeg,K11zzeg,K12xxeg,K12yyeg,K12zzeg,...
       K21xxeg,K21yyeg,K21zzeg,K22xxeg,K22yyeg,K22zzeg]=...
    solveUnitCell(seg,iElem,Parameters,UnitCell);
end

% Compute derivatives of effective parameters wrt saturation
if nsd==2
  [Gdseg,K11xxdseg,K11yydseg,K12xxdseg,K12yydseg,...
         K21xxdseg,K21yydseg,K22xxdseg,K22yydseg]=...
    solveUnitCell(seg*(1+ds)+(seg==0)*ds,iElem,Parameters,UnitCell);
  dGdseg=(Gdseg-Geg)./(seg*ds+(seg==0)*ds);
  dK11xxdseg=(K11xxdseg-K11xxeg)./(seg*ds+(seg==0)*ds);
  dK11yydseg=(K11yydseg-K11yyeg)./(seg*ds+(seg==0)*ds);
  dK12xxdseg=(K12xxdseg-K12xxeg)./(seg*ds+(seg==0)*ds);
  dK12yydseg=(K12yydseg-K12yyeg)./(seg*ds+(seg==0)*ds);
  dK21xxdseg=(K21xxdseg-K21xxeg)./(seg*ds+(seg==0)*ds);
  dK21yydseg=(K21yydseg-K21yyeg)./(seg*ds+(seg==0)*ds);
  dK22xxdseg=(K22xxdseg-K22xxeg)./(seg*ds+(seg==0)*ds);
  dK22yydseg=(K22yydseg-K22yyeg)./(seg*ds+(seg==0)*ds);
elseif nsd==3
  [Gdseg,K11xxdseg,K11yydseg,K11zzdseg,K12xxdseg,K12yydseg,K12zzdseg,...
         K21xxdseg,K21yydseg,K21zzdseg,K22xxdseg,K22yydseg,K22zzdseg]=...
    solveUnitCell(seg*(1+ds)+(seg==0)*ds,iElem,Parameters,UnitCell);
  dGdseg=(Gdseg-Geg)./(seg*ds+(seg==0)*ds);
  dK11xxdseg=(K11xxdseg-K11xxeg)./(seg*ds+(seg==0)*ds);
  dK11yydseg=(K11yydseg-K11yyeg)./(seg*ds+(seg==0)*ds);
  dK11zzdseg=(K11zzdseg-K11zzeg)./(seg*ds+(seg==0)*ds);
  dK12xxdseg=(K12xxdseg-K12xxeg)./(seg*ds+(seg==0)*ds);
  dK12yydseg=(K12yydseg-K12yyeg)./(seg*ds+(seg==0)*ds);
  dK12zzdseg=(K12zzdseg-K12zzeg)./(seg*ds+(seg==0)*ds);
  dK21xxdseg=(K21xxdseg-K21xxeg)./(seg*ds+(seg==0)*ds);
  dK21yydseg=(K21yydseg-K21yyeg)./(seg*ds+(seg==0)*ds);
  dK21zzdseg=(K21zzdseg-K21zzeg)./(seg*ds+(seg==0)*ds);
  dK22xxdseg=(K22xxdseg-K22xxeg)./(seg*ds+(seg==0)*ds);
  dK22yydseg=(K22yydseg-K22yyeg)./(seg*ds+(seg==0)*ds);
  dK22zzdseg=(K22zzdseg-K22zzeg)./(seg*ds+(seg==0)*ds);
end

% Print effective parameters and derivatives
if ReduceTo1GaussPoint && PrintGaussData
  fprintf('\n    Elem = %d',iElem);
  if nsd==2
    fprintf(['\t[s,G,K11xx,K11yy,K12xx,K12yy,',...
                    'K21xx,K21yy,K22xx,K22yy]=',...
             '[%.3f,%.2f,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             seg,Geg,K11xxeg,K11yyeg,K12xxeg,K12yyeg,...
                     K21xxeg,K21yyeg,K22xxeg,K22yyeg);
    fprintf(['\t[dGds,dK11xxds,dK11yyds,dK12xxds,dK12yyds,',...
                     'dK21xxds,dK21yyds,dK22xxds,dK22yyds]=',...
             '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             dGdseg,dK11xxdseg,dK11yydseg,dK12xxdseg,dK12yydseg,...
                    dK21xxdseg,dK21yydseg,dK22xxdseg,dK22yydseg);
  elseif nsd==3
    fprintf(['\t[s,G,K11xx,K11yy,K11zz,K12xx,K12yy,K12zz,',...
                    'K21xx,K21yy,K21zz,K22xx,K22yy,K22zz]=',...
             '[%.3f,%.2f,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             seg,Geg,K11xxeg,K11yyeg,K11zzeg,K12xxeg,K12yyeg,K12zzeg,...
                     K21xxeg,K21yyeg,K21zzeg,K22xxeg,K22yyeg,K22zzeg);
    fprintf(['\t[dGds,dK11xxds,dK11yyds,dK11zzds,dK12xxds,dK12yyds,dK12zzds,',...
                     'dK21xxds,dK21yyds,dK21zzds,dK22xxds,dK22yyds,dK22zzds]=',...
             '[%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e]'],...
             dGdseg,dK11xxdseg,dK11yydseg,dK11zzdseg,dK12xxdseg,dK12yydseg,dK12zzdseg,...
                    dK21xxdseg,dK21yydseg,dK21zzdseg,dK22xxdseg,dK22yydseg,dK22zzdseg);
  end
end

% Compute terms for speedup
gx=g(1);
gy=g(2);
if nsd==3
  gz=g(3);
end
vxeg=(K11xxeg+K21xxeg).*(-Qxeg+Hxeg+r1*gx)+(K12xxeg+K22xxeg).*(-Qxeg-Hxeg+r2*gx);
vyeg=(K11yyeg+K21yyeg).*(-Qyeg+Hyeg+r1*gy)+(K12yyeg+K22yyeg).*(-Qyeg-Hyeg+r2*gy);
if nsd==3
  vzeg=(K11zzeg+K21zzeg).*(-Qzeg+Hzeg+r1*gz)+(K12zzeg+K22zzeg).*(-Qzeg-Hzeg+r2*gz);
end
v2xeg=K21xxeg.*(-Qxeg+Hxeg+r1*gx)+K22xxeg.*(-Qxeg-Hxeg+r2*gx);
v2yeg=K21yyeg.*(-Qyeg+Hyeg+r1*gy)+K22yyeg.*(-Qyeg-Hyeg+r2*gy);
if nsd==3
  v2zeg=K21zzeg.*(-Qzeg+Hzeg+r1*gz)+K22zzeg.*(-Qzeg-Hzeg+r2*gz);
end
dvxdseg=(dK11xxdseg+dK21xxdseg).*(-Qxeg+Hxeg+r1*gx)+(dK12xxdseg+dK22xxdseg).*(-Qxeg-Hxeg+r2*gx);
dvydseg=(dK11yydseg+dK21yydseg).*(-Qyeg+Hyeg+r1*gy)+(dK12yydseg+dK22yydseg).*(-Qyeg-Hyeg+r2*gy);
if nsd==3
  dvzdseg=(dK11zzdseg+dK21zzdseg).*(-Qzeg+Hzeg+r1*gz)+(dK12zzdseg+dK22zzdseg).*(-Qzeg-Hzeg+r2*gz);
end
dv2xdseg=dK21xxdseg.*(-Qxeg+Hxeg+r1*gx)+dK22xxdseg.*(-Qxeg-Hxeg+r2*gx);
dv2ydseg=dK21yydseg.*(-Qyeg+Hyeg+r1*gy)+dK22yydseg.*(-Qyeg-Hyeg+r2*gy);
if nsd==3
  dv2zdseg=dK21zzdseg.*(-Qzeg+Hzeg+r1*gz)+dK22zzdseg.*(-Qzeg-Hzeg+r2*gz);
end
dvxdQxeg=-(K11xxeg+K21xxeg)-(K12xxeg+K22xxeg);
dvydQyeg=-(K11yyeg+K21yyeg)-(K12yyeg+K22yyeg);
if nsd==3
  dvzdQzeg=-(K11zzeg+K21zzeg)-(K12zzeg+K22zzeg);
end
dv2xdQxeg=-K21xxeg-K22xxeg;
dv2ydQyeg=-K21yyeg-K22yyeg;
if nsd==3
  dv2zdQzeg=-K21zzeg-K22zzeg;
end
dvxdHxeg=(K11xxeg+K21xxeg)-(K12xxeg+K22xxeg);
dvydHyeg=(K11yyeg+K21yyeg)-(K12yyeg+K22yyeg);
if nsd==3
  dvzdHzeg=(K11zzeg+K21zzeg)-(K12zzeg+K22zzeg);
end
dv2xdHxeg=K21xxeg-K22xxeg;
dv2ydHyeg=K21yyeg-K22yyeg;
if nsd==3
  dv2zdHzeg=K21zzeg-K22zzeg;
end

% Compute lhs
KQQ(ne1,ne1)=Me;
KQQ(ne2,ne2)=Me;
if nsd==3
  KQQ(ne3,ne3)=Me;
end

KQp(ne1,ne1)=Cxe;
KQp(ne2,ne1)=Cye;
if nsd==3
  KQp(ne3,ne1)=Cze;
end

KHH(ne1,ne1)=Me;
KHH(ne2,ne2)=Me;
if nsd==3
  KHH(ne3,ne3)=Me;
end

KHs(ne1,ne1)=NwexT*((dGdseg).*Ne);
KHs(ne2,ne1)=NweyT*((dGdseg).*Ne);
if nsd==3
  KHs(ne3,ne1)=NwezT*((dGdseg).*Ne);
end

Kps(ne1,ne1)=-NwexT*((dvxdseg).*Ne)...
             -NweyT*((dvydseg).*Ne);
if nsd==3
  Kps(ne1,ne1)=Kps(ne1,ne1)-NwezT*((dvzdseg).*Ne);
end

KpQ(ne1,ne1)=-NwexT*((dvxdQxeg).*Ne);
KpQ(ne1,ne2)=-NweyT*((dvydQyeg).*Ne);
if nsd==3
  KpQ(ne1,ne3)=-NwezT*((dvzdQzeg).*Ne);
end

KpH(ne1,ne1)=-NwexT*((dvxdHxeg).*Ne);
KpH(ne1,ne2)=-NweyT*((dvydHyeg).*Ne);
if nsd==3
  KpH(ne1,ne3)=-NwezT*((dvzdHzeg).*Ne);
end

if isTimeDependent
  Kss(ne1,ne1)=ep*alpha(1)/dt*Me;
end

Kss(ne1,ne1)=Kss(ne1,ne1)-NwexT*((dv2xdseg).*Ne)...
                         -NweyT*((dv2ydseg).*Ne);
if nsd==3
  Kss(ne1,ne1)=Kss(ne1,ne1)-NwezT*((dv2zdseg).*Ne);
end

KsQ(ne1,ne1)=-NwexT*((dv2xdQxeg).*Ne);
KsQ(ne1,ne2)=-NweyT*((dv2ydQyeg).*Ne);
if nsd==3
  KsQ(ne1,ne3)=-NwezT*((dv2zdQzeg).*Ne);
end

KsH(ne1,ne1)=-NwexT*((dv2xdHxeg).*Ne);
KsH(ne1,ne2)=-NweyT*((dv2ydHyeg).*Ne);
if nsd==3
  KsH(ne1,ne3)=-NwezT*((dv2zdHzeg).*Ne);
end

% Compute rhs
fQ(ne1,1)=-NweT*(Qxeg)...
          -NwexT*(peg);
fQ(ne2,1)=-NweT*(Qyeg)...
          -NweyT*(peg);
if nsd==3
  fQ(ne3,1)=-NweT*(Qzeg)...
            -NwezT*(peg);
end

fH(ne1,1)=-NweT*(Hxeg)...
          -NwexT*(Geg);
fH(ne2,1)=-NweT*(Hyeg)...
          -NweyT*(Geg);
if nsd==3
  fH(ne3,1)=-NweT*(Hzeg)...
            -NwezT*(Geg);
end

fp=fp+NwexT*(vxeg)...
     +NweyT*(vyeg);
if nsd==3
  fp=fp+NwezT*(vzeg);
end

fp=fp+NweT*(q1eg/r1+q2eg/r2);

if isTimeDependent
  fs=-NweT*(ep*[seg,soldeg]*alpha/dt);
end

fs=fs+NwexT*(v2xeg)...
     +NweyT*(v2yeg);
if nsd==3
  fs=fs+NwezT*(v2zeg);
end

fs=fs+NweT*(q2eg/r2);

% Faces loop
for iFace=1:NumElementFaces
  
  % Check need to compute face
  ComputeFace=true;
  
  % Compute face
  if ComputeFace
    % Compute weights at Gauss points
    FaceNodes=RefElement.FaceNodesElem;
    Xf=Xe(FaceNodes(iFace,:),:);
    [Nf,nx,ny,nz,wfg]=mapShapeFunctionsLinear(0,RefElement,Xf(1:nsd,:),nsd);

    % Reduce to only 1 Gauss point
    if ReduceTo1GaussPoint
      Nf=mean(Nf,1);
      wfg=sum(wfg,1);
    end
    
    % Check boundary
    isInterior=Faces.Interior(1,iFace);
    isDirichlet_p=Faces.Dirichlet_p(iFace);
    isDirichlet_s=Faces.Dirichlet_s(iFace);
    isNeumann_p=Faces.Neumann_p(iFace);
    isNeumann_s=Faces.Neumann_s(iFace);
    isDoNothing_p=Faces.DoNothing_p(iFace);
    isDoNothing_s=Faces.DoNothing_s(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nefU1=(iFace-1)*(1+1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefP1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefS1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
      nefU1=nefU1(order);
      nefU2=nefU2(order);
      nefP1=nefP1(order);
      nefS1=nefS1(order);
    end
    
    % Flip face (periodic (slave))
    if not(isInterior) && matchField(Faces,'Periodic') && Faces.Periodic(1,iFace)
      isInterior=true;
      Node2Match1stNode1=Faces.Periodic(2,iFace);
      FlipFace=max(Node2Match1stNode1);
      if FlipFace
        order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
        nefU1=nefU1(order);
        nefU2=nefU2(order);
        nefP1=nefP1(order);
        nefS1=nefS1(order);
      end
    end
    
    % Compute variables at nodes
    Qxf=Qxe(nf1);
    Qyf=Qye(nf1);
    if nsd==3
      Qzf=Qze(nf1);
    end
    Hxf=Hxe(nf1);
    Hyf=Hye(nf1);
    if nsd==3
      Hzf=Hze(nf1);
    end
    pf=pe(nf1);
    sf=se(nf1);
    Pf=Ue(nefU1);
    Sf=Ue(nefU2);
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    Qxfg=Nf*Qxf;
    Qyfg=Nf*Qyf;
    if nsd==3
      Qzfg=Nf*Qzf;
    end
    Hxfg=Nf*Hxf;
    Hyfg=Nf*Hyf;
    if nsd==3
      Hzfg=Nf*Hzf;
    end
    pfg=Nf*pf;
    sfg=Nf*sf;
    Pfg=Nf*Pf;
    Sfg=Nf*Sf;
    if isDirichlet_p
      pDfg=pD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isDirichlet_s
      sDfg=sD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isNeumann_p || isNeumann_s
      f1Nfg=f1N(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      f2Nfg=f2N(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
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
    
    % Compute effective parameters
    if nsd==2
      [Gfg,K11xxfg,K11yyfg,K12xxfg,K12yyfg,...
           K21xxfg,K21yyfg,K22xxfg,K22yyfg]=...
        solveUnitCell(Sfg,iElem,Parameters,UnitCell);
    elseif nsd==3
      [Gfg,K11xxfg,K11yyfg,K11zzfg,K12xxfg,K12yyfg,K12zzfg,...
           K21xxfg,K21yyfg,K21zzfg,K22xxfg,K22yyfg,K22zzfg]=...
        solveUnitCell(Sfg,iElem,Parameters,UnitCell);
    end
    
    % Compute derivatives of effective parameters wrt saturation
    if nsd==2
      [GdSfg,K11xxdSfg,K11yydSfg,K12xxdSfg,K12yydSfg,...
             K21xxdSfg,K21yydSfg,K22xxdSfg,K22yydSfg]=...
        solveUnitCell(Sfg*(1+ds)+(Sfg==0)*ds,iElem,Parameters,UnitCell);
      dGdSfg=(GdSfg-Gfg)./(Sfg*ds+(Sfg==0)*ds);
      dK11xxdSfg=(K11xxdSfg-K11xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK11yydSfg=(K11yydSfg-K11yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK12xxdSfg=(K12xxdSfg-K12xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK12yydSfg=(K12yydSfg-K12yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK21xxdSfg=(K21xxdSfg-K21xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK21yydSfg=(K21yydSfg-K21yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK22xxdSfg=(K22xxdSfg-K22xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK22yydSfg=(K22yydSfg-K22yyfg)./(Sfg*ds+(Sfg==0)*ds);
    elseif nsd==3
      [GdSfg,K11xxdSfg,K11yydSfg,K11zzdSfg,K12xxdSfg,K12yydSfg,K12zzdSfg,...
             K21xxdSfg,K21yydSfg,K21zzdSfg,K22xxdSfg,K22yydSfg,K22zzdSfg]=...
        solveUnitCell(Sfg*(1+ds)+(Sfg==0)*ds,iElem,Parameters,UnitCell);
      dGdSfg=(GdSfg-Gfg)./(Sfg*ds+(Sfg==0)*ds);
      dK11xxdSfg=(K11xxdSfg-K11xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK11yydSfg=(K11yydSfg-K11yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK11zzdSfg=(K11zzdSfg-K11zzfg)./(Sfg*ds+(Sfg==0)*ds);
      dK12xxdSfg=(K12xxdSfg-K12xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK12yydSfg=(K12yydSfg-K12yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK12zzdSfg=(K12zzdSfg-K12zzfg)./(Sfg*ds+(Sfg==0)*ds);
      dK21xxdSfg=(K21xxdSfg-K21xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK21yydSfg=(K21yydSfg-K21yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK21zzdSfg=(K21zzdSfg-K21zzfg)./(Sfg*ds+(Sfg==0)*ds);
      dK22xxdSfg=(K22xxdSfg-K22xxfg)./(Sfg*ds+(Sfg==0)*ds);
      dK22yydSfg=(K22yydSfg-K22yyfg)./(Sfg*ds+(Sfg==0)*ds);
      dK22zzdSfg=(K22zzdSfg-K22zzfg)./(Sfg*ds+(Sfg==0)*ds);
    end
    
    % Compute terms for speedup
    Vxfg=(K11xxfg+K21xxfg).*(-Qxfg+Hxfg+r1*gx)+(K12xxfg+K22xxfg).*(-Qxfg-Hxfg+r2*gx);
    Vyfg=(K11yyfg+K21yyfg).*(-Qyfg+Hyfg+r1*gy)+(K12yyfg+K22yyfg).*(-Qyfg-Hyfg+r2*gy);
    if nsd==3
      Vzfg=(K11zzfg+K21zzfg).*(-Qzfg+Hzfg+r1*gz)+(K12zzfg+K22zzfg).*(-Qzfg-Hzfg+r2*gz);
    end
    V2xfg=K21xxfg.*(-Qxfg+Hxfg+r1*gx)+K22xxfg.*(-Qxfg-Hxfg+r2*gx);
    V2yfg=K21yyfg.*(-Qyfg+Hyfg+r1*gy)+K22yyfg.*(-Qyfg-Hyfg+r2*gy);
    if nsd==3
      V2zfg=K21zzfg.*(-Qzfg+Hzfg+r1*gz)+K22zzfg.*(-Qzfg-Hzfg+r2*gz);
    end
    dVxdSfg=(dK11xxdSfg+dK21xxdSfg).*(-Qxfg+Hxfg+r1*gx)+(dK12xxdSfg+dK22xxdSfg).*(-Qxfg-Hxfg+r2*gx);
    dVydSfg=(dK11yydSfg+dK21yydSfg).*(-Qyfg+Hyfg+r1*gy)+(dK12yydSfg+dK22yydSfg).*(-Qyfg-Hyfg+r2*gy);
    if nsd==3
      dVzdSfg=(dK11zzdSfg+dK21zzdSfg).*(-Qzfg+Hzfg+r1*gz)+(dK12zzdSfg+dK22zzdSfg).*(-Qzfg-Hzfg+r2*gz);
    end
    dV2xdSfg=dK21xxdSfg.*(-Qxfg+Hxfg+r1*gx)+dK22xxdSfg.*(-Qxfg-Hxfg+r2*gx);
    dV2ydSfg=dK21yydSfg.*(-Qyfg+Hyfg+r1*gy)+dK22yydSfg.*(-Qyfg-Hyfg+r2*gy);
    if nsd==3
      dV2zdSfg=dK21zzdSfg.*(-Qzfg+Hzfg+r1*gz)+dK22zzdSfg.*(-Qzfg-Hzfg+r2*gz);
    end
    dVxdQxfg=-(K11xxfg+K21xxfg)-(K12xxfg+K22xxfg);
    dVydQyfg=-(K11yyfg+K21yyfg)-(K12yyfg+K22yyfg);
    if nsd==3
      dVzdQzfg=-(K11zzfg+K21zzfg)-(K12zzfg+K22zzfg);
    end
    dV2xdQxfg=-K21xxfg-K22xxfg;
    dV2ydQyfg=-K21yyfg-K22yyfg;
    if nsd==3
      dV2zdQzfg=-K21zzfg-K22zzfg;
    end
    dVxdHxfg=(K11xxfg+K21xxfg)-(K12xxfg+K22xxfg);
    dVydHyfg=(K11yyfg+K21yyfg)-(K12yyfg+K22yyfg);
    if nsd==3
      dVzdHzfg=(K11zzfg+K21zzfg)-(K12zzfg+K22zzfg);
    end
    dV2xdHxfg=K21xxfg-K22xxfg;
    dV2ydHyfg=K21yyfg-K22yyfg;
    if nsd==3
      dV2zdHzfg=K21zzfg-K22zzfg;
    end
    if nsd==2
      Vnfg=Vxfg.*nx+Vyfg.*ny;
      V2nfg=V2xfg.*nx+V2yfg.*ny;
      dVndSfg=dVxdSfg.*nx+dVydSfg.*ny;
      dV2ndSfg=dV2xdSfg.*nx+dV2ydSfg.*ny;
    elseif nsd==3
      Vnfg=Vxfg.*nx+Vyfg.*ny+Vzfg.*nz;
      V2nfg=V2xfg.*nx+V2yfg.*ny+V2zfg.*nz;
      dVndSfg=dVxdSfg.*nx+dVydSfg.*ny+dVzdSfg.*nz;
      dV2ndSfg=dV2xdSfg.*nx+dV2ydSfg.*ny+dV2zdSfg.*nz;
    end
    
    % Compute lhs
    KQP(nf1,nefP1)=KQP(nf1,nefP1)-Mfnx;
    KQP(nf2,nefP1)=KQP(nf2,nefP1)-Mfny;
    if nsd==3
      KQP(nf3,nefP1)=KQP(nf3,nefP1)-Mfnz;
    end

    KHS(nf1,nefS1)=KHS(nf1,nefS1)-NwfT*((nx.*dGdSfg).*Nf);
    KHS(nf2,nefS1)=KHS(nf2,nefS1)-NwfT*((ny.*dGdSfg).*Nf);
    if nsd==3
      KHS(nf3,nefS1)=KHS(nf3,nefS1)-NwfT*((nz.*dGdSfg).*Nf);
    end
    
    KpS(nf1,nefS1)=KpS(nf1,nefS1)+NwfT*((dVndSfg).*Nf);

    KpQ(nf1,nf1)=KpQ(nf1,nf1)+NwfT*((dVxdQxfg.*nx).*Nf);
    KpQ(nf1,nf2)=KpQ(nf1,nf2)+NwfT*((dVydQyfg.*ny).*Nf);
    if nsd==3
      KpQ(nf1,nf3)=KpQ(nf1,nf3)+NwfT*((dVzdQzfg.*nz).*Nf);
    end

    KpH(nf1,nf1)=KpH(nf1,nf1)+NwfT*((dVxdHxfg.*nx).*Nf);
    KpH(nf1,nf2)=KpH(nf1,nf2)+NwfT*((dVydHyfg.*ny).*Nf);
    if nsd==3
      KpH(nf1,nf3)=KpH(nf1,nf3)+NwfT*((dVzdHzfg.*nz).*Nf);
    end

    Kpp(nf1,nf1)=Kpp(nf1,nf1)+tauP*Mf;

    KpP(nf1,nefP1)=KpP(nf1,nefP1)-tauP*Mf;

    KsS(nf1,nefS1)=KsS(nf1,nefS1)+NwfT*((dV2ndSfg).*Nf);

    KsQ(nf1,nf1)=KsQ(nf1,nf1)+NwfT*((dV2xdQxfg.*nx).*Nf);
    KsQ(nf1,nf2)=KsQ(nf1,nf2)+NwfT*((dV2ydQyfg.*ny).*Nf);
    if nsd==3
      KsQ(nf1,nf3)=KsQ(nf1,nf3)+NwfT*((dV2zdQzfg.*nz).*Nf);
    end

    KsH(nf1,nf1)=KsH(nf1,nf1)+NwfT*((dV2xdHxfg.*nx).*Nf);
    KsH(nf1,nf2)=KsH(nf1,nf2)+NwfT*((dV2ydHyfg.*ny).*Nf);
    if nsd==3
      KsH(nf1,nf3)=KsH(nf1,nf3)+NwfT*((dV2zdHzfg.*nz).*Nf);
    end

    Kss(nf1,nf1)=Kss(nf1,nf1)+tauS*Mf;

    KsS(nf1,nefS1)=KsS(nf1,nefS1)-tauS*Mf;

    if isInterior || isNeumann_p
      KPS(nefP1,nefS1)=-NwfT*((dVndSfg).*Nf);

      KPQ(nefP1,nf1)=KPQ(nefP1,nf1)-NwfT*((dVxdQxfg.*nx).*Nf);
      KPQ(nefP1,nf2)=KPQ(nefP1,nf2)-NwfT*((dVydQyfg.*ny).*Nf);
      if nsd==3
        KPQ(nefP1,nf3)=KPQ(nefP1,nf3)-NwfT*((dVzdQzfg.*nz).*Nf);
      end

      KPH(nefP1,nf1)=KPH(nefP1,nf1)-NwfT*((dVxdHxfg.*nx).*Nf);
      KPH(nefP1,nf2)=KPH(nefP1,nf2)-NwfT*((dVydHyfg.*ny).*Nf);
      if nsd==3
        KPH(nefP1,nf3)=KPH(nefP1,nf3)-NwfT*((dVzdHzfg.*nz).*Nf);
      end

      KPp(nefP1,nf1)=KPp(nefP1,nf1)-tauP*Mf;

      KPP(nefP1,nefP1)=tauP*Mf;
    elseif isDirichlet_p
      KPP(nefP1,nefP1)=Mf;
    elseif isDoNothing_p
      KPp(nefP1,nf1)=KPp(nefP1,nf1)-tauP*Mf;

      KPP(nefP1,nefP1)=tauP*Mf;
    end

    if isInterior || isNeumann_s
      KSS(nefS1,nefS1)=KSS(nefS1,nefS1)-NwfT*((dV2ndSfg).*Nf);

      KSQ(nefS1,nf1)=KSQ(nefS1,nf1)-NwfT*((dV2xdQxfg.*nx).*Nf);
      KSQ(nefS1,nf2)=KSQ(nefS1,nf2)-NwfT*((dV2ydQyfg.*ny).*Nf);
      if nsd==3
        KSQ(nefS1,nf3)=KSQ(nefS1,nf3)-NwfT*((dV2zdQzfg.*nz).*Nf);
      end

      KSH(nefS1,nf1)=KSH(nefS1,nf1)-NwfT*((dV2xdHxfg.*nx).*Nf);
      KSH(nefS1,nf2)=KSH(nefS1,nf2)-NwfT*((dV2ydHyfg.*ny).*Nf);
      if nsd==3
        KSH(nefS1,nf3)=KSH(nefS1,nf3)-NwfT*((dV2zdHzfg.*nz).*Nf);
      end

      KSs(nefS1,nf1)=KSs(nefS1,nf1)-tauS*Mf;

      KSS(nefS1,nefS1)=KSS(nefS1,nefS1)+tauS*Mf;
    elseif isDirichlet_s
      KSS(nefS1,nefS1)=Mf;
    elseif isDoNothing_s
      KSs(nefS1,nf1)=KSs(nefS1,nf1)-tauS*Mf;

      KSS(nefS1,nefS1)=KSS(nefS1,nefS1)+tauS*Mf;
    end
    
    % Compute rhs
    fQ(nf1,1)=fQ(nf1,1)+NwfT*(nx.*Pfg);
    fQ(nf2,1)=fQ(nf2,1)+NwfT*(ny.*Pfg);
    if nsd==3
      fQ(nf3,1)=fQ(nf3,1)+NwfT*(nz.*Pfg);
    end

    fH(nf1,1)=fH(nf1,1)+NwfT*(nx.*Gfg);
    fH(nf2,1)=fH(nf2,1)+NwfT*(ny.*Gfg);
    if nsd==3
      fH(nf3,1)=fH(nf3,1)+NwfT*(nz.*Gfg);
    end

    fp(nf1,1)=fp(nf1,1)-NwfT*(Vnfg+tauP*(pfg-Pfg));

    fs(nf1,1)=fs(nf1,1)-NwfT*(V2nfg+tauS*(sfg-Sfg));

    if isInterior
      fP(nefP1,1)=+NwfT*(Vnfg+tauP*(pfg-Pfg));
    elseif isNeumann_p
      fP(nefP1,1)=+NwfT*(Vnfg+tauP*(pfg-Pfg)+f1Nfg/r1+f2Nfg/r2);
    elseif isDirichlet_p
      fP(nefP1,1)=-NwfT*(Pfg-pDfg);
    elseif isDoNothing_p
      fP(nefP1,1)=+NwfT*(tauP*(pfg-Pfg));
    end

    if isInterior
      fS(nefS1,1)=+NwfT*(V2nfg+tauS*(sfg-Sfg));
    elseif isNeumann_s
      fS(nefS1,1)=+NwfT*(V2nfg+tauS*(sfg-Sfg)+f2Nfg/r2);
    elseif isDirichlet_s
      fS(nefS1,1)=-NwfT*(Sfg-sDfg);
    elseif isDoNothing_s
      fS(nefS1,1)=+NwfT*(tauS*(sfg-Sfg));
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iQ=1:nsd*NumElementNodes;
iH=iQ(end)+(1:nsd*NumElementNodes);
ip=iH(end)+(1:NumElementNodes);
is=ip(end)+(1:NumElementNodes);
iP=reshape((0:NumElementFaces-1)*(1+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iS=reshape((0:NumElementFaces-1)*(1+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((nsd+nsd+1+1)*NumElementNodes,(nsd+nsd+1+1)*NumElementNodes);
LhsLG=zeros((nsd+nsd+1+1)*NumElementNodes,(1+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((1+1)*NumElementFaces*NumFaceNodes,(nsd+nsd+1+1)*NumElementNodes);
LhsGG=zeros((1+1)*NumElementFaces*NumFaceNodes,(1+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((nsd+nsd+1+1)*NumElementNodes,1);
RhsG=zeros((1+1)*NumElementFaces*NumFaceNodes,1);

% Remove pressure equation
if not(SolvePressureEquation)
  KQQ=eye(size(KQQ));
  KQp=zeros(size(KQp));
  KpQ=zeros(size(KpQ));
  KpH=zeros(size(KpH));
  Kpp=eye(size(Kpp));
  Kps=zeros(size(Kps));
  KQP=zeros(size(KQP));
  KpP=zeros(size(KpP));
  KpS=zeros(size(KpS));
  fQ=zeros(size(fQ));
  fp=zeros(size(fp));
  KPQ=zeros(size(KPQ));
  KPH=zeros(size(KPH));
  KPp=zeros(size(KPp));
  KSQ=zeros(size(KSQ));
  KPP=eye(size(KPP));
  KPS=zeros(size(KPS));
  fP=zeros(size(fP));
end

% Lhs local-local
LhsLL(iQ,iQ)=KQQ;
LhsLL(iQ,ip)=KQp;
LhsLL(iH,iH)=KHH;
LhsLL(iH,is)=KHs;
LhsLL(ip,iQ)=KpQ;
LhsLL(ip,iH)=KpH;
LhsLL(ip,ip)=Kpp;
LhsLL(ip,is)=Kps;
LhsLL(is,iQ)=KsQ;
LhsLL(is,iH)=KsH;
LhsLL(is,is)=Kss;

% Lhs local-global
LhsLG(iQ,iP)=KQP;
LhsLG(iH,iS)=KHS;
LhsLG(ip,iP)=KpP;
LhsLG(ip,iS)=KpS;
LhsLG(is,iS)=KsS;

% Rhs local
RhsL(iQ,1)=fQ;
RhsL(iH,1)=fH;
RhsL(ip,1)=fp;
RhsL(is,1)=fs;

% Lhs global-local
LhsGL(iP,iQ)=KPQ;
LhsGL(iP,iH)=KPH;
LhsGL(iP,ip)=KPp;
LhsGL(iS,iQ)=KSQ;
LhsGL(iS,iH)=KSH;
LhsGL(iS,is)=KSs;

% Lhs global-global
LhsGG(iP,iP)=KPP;
LhsGG(iP,iS)=KPS;
LhsGG(iS,iS)=KSS;

% Rhs global
RhsG(iP,1)=fP;
RhsG(iS,1)=fS;

% Equilibrate
if EquilibrateLocalProblems
  [EquilP,EquilR,EquilC]=equilibrate(LhsLL);
  LhsLL=EquilR*EquilP*LhsLL*EquilC;
  LhsLGRhsL=EquilR*EquilP*[LhsLG,RhsL];
else
  LhsLGRhsL=[LhsLG,RhsL];
end

% Matrix and vector for local problem
MatVecLocal=LhsLL\LhsLGRhsL;

% Un-equilibrate
if EquilibrateLocalProblems
  MatVecLocal=EquilC*MatVecLocal;
end

% Extract matrix for local problem
MatLocal=MatVecLocal(:,1:end-1);

% Extract vector for local problem
VecLocal=MatVecLocal(:,end);

% Lhs for global problem
LhsGlobal=LhsGG-LhsGL*MatLocal;

% Rhs for global problem
RhsGlobal=RhsG-LhsGL*VecLocal;

% Compute fill factor
Fill=sum(seg.*weg);
Volume=sum(weg);

end

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  iElem,Nodes,SolutionLocal,Parameters,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
ds=Parameters.RelativeFiniteDifference;
solveUnitCell=Parameters.solveUnitCell;
Xe=Nodes';

% Retrieve unit cell
UnitCell=Parameters.UnitCell.Value;

% Get solution
Qe=reshape(SolutionLocal(:,1:nsd),[],1);
He=reshape(SolutionLocal(:,nsd+(1:nsd)),[],1);
pe=reshape(SolutionLocal(:,nsd+nsd+1),[],1);
se=reshape(SolutionLocal(:,nsd+nsd+1+1),[],1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions('Element',RefElement.PostLow,RefElement.Post,Xe,nsd);
N1e=ones(length(weg),1);

% Indices
nle1=1:NumElementNodes;
nle2=nle1+NumElementNodes;
nle3=nle2+NumElementNodes;

% Compute variables at nodes
Qxe=Qe(nle1);
Qye=Qe(nle2);
if nsd==3
  Qze=Qe(nle3);
end
Hxe=He(nle1);
Hye=He(nle2);
if nsd==3
  Hze=He(nle3);
end

% Compute variables at Gauss points
Qxeg=Nle*Qxe;
Qyeg=Nle*Qye;
if nsd==3
  Qzeg=Nle*Qze;
end
Hxeg=Nle*Hxe;
Hyeg=Nle*Hye;
if nsd==3
  Hzeg=Nle*Hze;
end
peg=Nle*pe;
seg=Nle*se;

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

% Compute effective parameters
[Geg]=solveUnitCell(seg,iElem,Parameters,UnitCell);

% Compute derivatives of effective parameters wrt saturation
[Gdseg]=solveUnitCell(seg*(1+ds)+(seg==0)*ds,iElem,Parameters,UnitCell);
dGdseg=(Gdseg-Geg)./(seg*ds+(seg==0)*ds);
if all(dGdseg==0)
  dGdseg=ones(size(dGdseg));
end

% Compute lhs
KPpp=Kxxe...
    +Kyye;
if nsd==3
  KPpp=KPpp+Kzze;
end

KSpp=NwexT*((dGdseg).*Nex)...
    +NweyT*((dGdseg).*Ney);
if nsd==3
  KSpp=KSpp+NwezT*((dGdseg).*Nez);
end

KPtp=Nw1eT*(Ne);

KStp=Nw1eT*(Ne);

% Compute rhs
fPp=NwexT*(Qxeg)...
   +NweyT*(Qyeg);
if nsd==3
  fPp=fPp+NwezT*(Qzeg);
end

fSp=NwexT*(Hxeg)...
   +NweyT*(Hyeg);
if nsd==3
  fSp=fSp+NwezT*(Hzeg);
end

fPt=Nw1eT*(peg);

fSt=Nw1eT*(seg);

% Indices
iPp=1:NumElementNodesPost;
iPt=iPp(end)+1;
iSp=iPt(end)+(1:NumElementNodesPost); iSp2=iPp(end)+(1:NumElementNodesPost);
iSt=iSp(end)+1;

% Initialization of lhs and rhs
LhsPost=zeros(NumElementNodesPost+1+NumElementNodesPost+1,NumElementNodesPost+NumElementNodesPost);
RhsPost=zeros(NumElementNodesPost+1+NumElementNodesPost+1,1);

% Lhs for post-processing
LhsPost(iPp,iPp)=KPpp;
LhsPost(iPt,iPp)=KPtp;
LhsPost(iSp,iSp2)=KSpp;
LhsPost(iSt,iSp2)=KStp;

% Rhs for post-processing
RhsPost(iPp,1)=fPp;
RhsPost(iPt,1)=fPt;
RhsPost(iSp,1)=fSp;
RhsPost(iSt,1)=fSt;

end