classdef WeaklyCompressibleFlowVP_HDG < Formulation  
  
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
      
      % Consider reference pressure instead of zero
      for iFace=1:Sizes(iD).NumFaces
        Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp.*Sizes(iD).NumFaceNodes+...
          Sizes(iD).NumSpaceDim*Sizes(iD).NumFaceNodes+(1:Sizes(iD).NumFaceNodes))=...
          Parameters(iD).ReferencePressure;
      end
      Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim)=...
        Parameters(iD).ReferencePressure;
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
          Parameters(iD).Velocity(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).Pressure(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
              Parameters(iD).Velocity(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).Pressure(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
          end
        end
      end
    end
    
    %% Build block
    function [Block,Elements]=buildBlock(~,iD1,Block,Elements,~,Parameters,~,~,Time,RefElement,...
        Sizes)
      NodesElem=Elements(iD1).Nodes;
      NodesInitialElem=Elements(iD1).NodesInitial;
      FacesElem=Elements(iD1).Faces;
      SolutionGlobalElem=Elements(iD1).SolutionGlobal;
      SolutionLocalElem=Elements(iD1).SolutionLocal;
      SolutionOldElem=Elements(iD1).SolutionOld;
      if length(Parameters)>1
        iD2=setdiff(1:2,iD1);
        SolutionGlobalElemCoupled=Elements(iD2).SolutionGlobal;
        SolutionOldElemCoupled=Elements(iD2).SolutionOld;
      else
        SolutionGlobalElemCoupled=cell(Sizes(iD1).NumElements,1);
        SolutionOldElemCoupled=cell(Sizes(iD1).NumElements,1);
      end
      LhsCoef=zeros(Sizes(iD1).NumElementLhsCoef,Sizes(iD1).NumElements);
      RhsCoef=zeros(Sizes(iD1).NumElementRhsCoef,Sizes(iD1).NumElements);
      MatLocal=cell(Sizes(iD1).NumElements,1);
      VecLocal=cell(Sizes(iD1).NumElements,1);
      parfor iElem=1:Sizes(iD1).NumElements
        [LhsGlobalElem,RhsGlobalElem,MatLocalElem,VecLocalElem]=...
          buildBlockElement(iD1,NodesElem{iElem},NodesInitialElem{iElem},FacesElem(iElem),...
          SolutionGlobalElem{iElem},SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          SolutionGlobalElemCoupled{iElem},SolutionOldElemCoupled{iElem},...
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
    
    %% Do coupling
    function [Block]=doCoupling(~,iD1,iD2,Block,~,~,~,~,~,~,~,~)
        Block(iD1,iD2).LhsGlobal=fsparse(Block(iD1,iD2).LhsRowIndices,Block(iD1,iD2).LhsColIndices,0);
    end
    
    %% Store results
    function [Results]=storeResults(~,iD,iST,Results,Block,~,Parameters,~,Time,~,Sizes)
      if iST==1
        Results(iD).Time=[];
        Results(iD).ScaledStrainRate=[];
        Results(iD).Velocity=[];
        Results(iD).Pressure=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledStrainRate(:,:,iST)=Block(iD,iD).SolutionLocal(:,1:Sizes(iD).NumVoigtComp);
      Results(iD).Velocity(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        (1:Sizes(iD).NumSpaceDim));
      Results(iD).Pressure(:,:,iST)=Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+...
        Sizes(iD).NumSpaceDim+1);
      if strcmp(Parameters(iD).PostProcessingHDG,'yes')
        Results(iD).VelocityPost=Block(iD,iD).SolutionPost;
      end
    end
    
    %% Data for Paraview
    function [PointData,CellData]=dataForParaview(~,Results,Parameters,~,Sizes,isPostProcess)
      
      if not(isPostProcess)
        
        % Write scaled strain rate
        L=Results.ScaledStrainRate(:,:,end);
        PointData=sprintf('\nTENSORS Scaled_strain_rate float\n');
        if Sizes.NumSpaceDim==2
          PointData=[PointData,sprintf(['%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f\n'],[L';zeros(6,size(L,1))])];
        elseif Sizes.NumSpaceDim==3
          PointData=[PointData,sprintf(['%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f ',...
                                        '%.12f %.12f %.12f\n'],[L';zeros(3,size(L,1))])];
        end
        
        % Write velocity
        v=Results.Velocity(:,:,end);
        PointData=[PointData,sprintf('\nVECTORS Velocity float\n')];
        if Sizes.NumSpaceDim==2
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',[v';zeros(1,size(v,1))])];
        elseif Sizes.NumSpaceDim==3
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',v')];
        end
        
        % Write pressure
        p=Results.Pressure(:,:,end);
        PointData=[PointData,sprintf('\nSCALARS Pressure float\n'),...
                   sprintf('LOOKUP_TABLE default\n'),...
                   sprintf('%.12f\n',p')];
        
      elseif isPostProcess
        
        % Write postprocessed velocity
        vp=Results.VelocityPost(:,:,end);
        PointData=sprintf('\nVECTORS Velocity_post float\n');
        if Sizes.NumSpaceDim==2
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',[vp';zeros(1,size(vp,1))])];
        elseif Sizes.NumSpaceDim==3
          PointData=[PointData,sprintf('%.12f %.12f %.12f\n',vp')];
        end
        
      end
      
      % Write reference density
      rho_ref=Parameters.ReferenceDensity;
      rho_ref_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        rho_ref_DG(iElem,:)=rho_ref;
      end
      CellData=[sprintf('\nSCALARS Reference_density float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',rho_ref_DG')];
      
      % Write reference pressure
      p_ref=Parameters.ReferencePressure;
      p_ref_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        p_ref_DG(iElem,:)=p_ref;
      end
      CellData=[CellData,sprintf('\nSCALARS Reference_pressure float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',p_ref_DG')];
      
      % Write compressibility coefficient
      epsilon=Parameters.CompressibilityCoeff;
      epsilon_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        epsilon_DG(iElem,:)=epsilon;
      end
      CellData=[CellData,sprintf('\nSCALARS Compressibility_coefficient float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',epsilon_DG')];
      
      % Write dynamic viscosity
      mu=Parameters.DynamicViscosity;
      mu_DG=zeros(Sizes.NumElements,1);
      for iElem=1:Sizes.NumElements
        mu_DG(iElem,:)=mu;
      end
      CellData=[CellData,sprintf('\nSCALARS Dynamic_viscosity float\n'),...
                sprintf('LOOKUP_TABLE default\n'),...
                sprintf('%.12f\n',mu_DG')];
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  iD1,Nodes,NodesInitial,Faces,SolutionGlobal,SolutionLocal,SolutionOld,...
  SolutionGlobalCoupled,SolutionOldCoupled,Parameters,Time,RefElement,Sizes)

% Get general parameters
isConvectiveFlow=strcmp(Parameters(iD1).ConvectiveFlow,'yes');
isArbitraryLagrangianEulerian=strcmp(Parameters(iD1).ArbitraryLagrangianEulerian,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes(iD1).NumSpaceDim;
msd=nsd*(nsd+1)/2;
NumElementNodes=Sizes(iD1).NumElementNodes;
NumElementFaces=Sizes(iD1).NumElementFaces;
NumFaceNodes=Sizes(iD1).NumFaceNodes;
mu=Parameters(iD1).DynamicViscosity;
lambda=-2/3*Parameters(iD1).DynamicViscosity;
r=Parameters(iD1).Density;
drdp=Parameters(iD1).DDensityDPressure;
vD=Parameters(iD1).Velocity;
pD=Parameters(iD1).Pressure;
tN=Parameters(iD1).Traction;
Rc=Parameters(iD1).ResidualContinuity;
f=Parameters(iD1).Force;
Xe=Nodes';
X0e=NodesInitial';
tauV=Parameters(iD1).StabVelocity;
tauP=Parameters(iD1).StabPressure;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
ve=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
pe=reshape(SolutionLocal(:,msd+nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  volde=reshape(SolutionOld(:,msd+(1:nsd),:),[],BDFo);
  polde=reshape(SolutionOld(:,msd+nsd+1,:),[],BDFo);
end
if isArbitraryLagrangianEulerian && isTimeDependent
  if length(Parameters)==1
    ue=reshape(Parameters(iD1).Displacement(X0e(:,1),X0e(:,2),X0e(:,3),t),[],1);
    uolde=zeros(nsd*NumElementNodes,BDFo);
    for iBDF=1:BDFo
      uolde(:,iBDF)=reshape(Parameters(iD1).Displacement(X0e(:,1),X0e(:,2),X0e(:,3),t-iBDF*dt),...
        [],1);
    end
  else
    ue=reshape(SolutionGlobalCoupled,[],1);
    uolde=reshape(SolutionOldCoupled(:,:,1:BDFo),[],BDFo);
  end
  ae=1/dt*ue*alpha(1)+1/dt*uolde*alpha(2:BDFo+1,1);
else
  ae=ve*0;
end

% Initialize lhs
KLL=zeros(msd*NumElementNodes,msd*NumElementNodes);
KLv=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KLV=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kvv=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kvp=zeros(nsd*NumElementNodes,NumElementNodes);
KvV=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KvP=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Kpv=zeros(NumElementNodes,nsd*NumElementNodes);
Kpp=zeros(NumElementNodes,NumElementNodes);
KpV=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KpP=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KVL=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KVv=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KVV=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KVP=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KPp=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KPP=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fL=zeros(msd*NumElementNodes,1);
fv=zeros(nsd*NumElementNodes,1);
fp=zeros(NumElementNodes,1);
fV=zeros(nsd*NumElementFaces*NumFaceNodes,1);
fP=zeros(NumElementFaces*NumFaceNodes,1);

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
  Lxxe=Le(ne1);
  Lyye=Le(ne2);
  Lxye=Le(ne3);
elseif nsd==3
  Lxxe=Le(ne1);
  Lyye=Le(ne2);
  Lzze=Le(ne3);
  Lxye=Le(ne4);
  Lxze=Le(ne5);
  Lyze=Le(ne6);
end
vxe=ve(ne1);
vye=ve(ne2);
if nsd==3
  vze=ve(ne3);
end
if isTimeDependent
  voldxe=volde(ne1,:);
  voldye=volde(ne2,:);
  if nsd==3
    voldze=volde(ne3,:);
  end
end
axe=ae(ne1);
aye=ae(ne2);
if nsd==3
  aze=ae(ne3);
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
if nsd==2
  Lxxeg=Ne*Lxxe;
  Lyyeg=Ne*Lyye;
  Lxyeg=Ne*Lxye;
elseif nsd==3
  Lxxeg=Ne*Lxxe;
  Lyyeg=Ne*Lyye;
  Lzzeg=Ne*Lzze;
  Lxyeg=Ne*Lxye;
  Lxzeg=Ne*Lxze;
  Lyzeg=Ne*Lyze;
end
vxeg=Ne*vxe;
vyeg=Ne*vye;
if nsd==3
  vzeg=Ne*vze;
end
peg=Ne*pe;
if isTimeDependent
  voldxeg=Ne*voldxe;
  voldyeg=Ne*voldye;
  if nsd==3
    voldzeg=Ne*voldze;
  end
  poldeg=Ne*polde;
end
Rceg=Rc(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
feg=f(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
fxeg=feg(:,1);
fyeg=feg(:,2);
if nsd==3
  fzeg=feg(:,3);
end
reg=r(peg);
drdpeg=drdp(peg);
if isTimeDependent
  dvxdteg=1/dt*vxeg*alpha(1)+1/dt*voldxeg*alpha(2:BDFo+1,1);
  dvydteg=1/dt*vyeg*alpha(1)+1/dt*voldyeg*alpha(2:BDFo+1,1);
  if nsd==3
    dvzdteg=1/dt*vzeg*alpha(1)+1/dt*voldzeg*alpha(2:BDFo+1,1);
  end
  dpdteg=1/dt*peg*alpha(1)+1/dt*poldeg*alpha(2:BDFo+1,1);
end
axeg=Ne*axe;
ayeg=Ne*aye;
if nsd==3
  azeg=Ne*aze;
end
if isArbitraryLagrangianEulerian
  if nsd==2
    divaeg=Nex*axe+Ney*aye;
  else
    divaeg=Nex*axe+Ney*aye+Nez*aze;
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

% Compute linearization of viscous stress
if nsd==2
  Voigt1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  Voigt2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  Voigt3=sqrt(mu);
elseif nsd==3
  Voigt1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  Voigt2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  Voigt3=sqrt(mu);
end

% Compute lhs
KLL(ne1,ne1)=-Me;
KLL(ne2,ne2)=-Me;
KLL(ne3,ne3)=-Me;
if nsd==3
  KLL(ne4,ne4)=-Me;
  KLL(ne5,ne5)=-Me;
  KLL(ne6,ne6)=-Me;
end

if nsd==2
  KLv(ne1,ne1)=Voigt1*Cxe;
  KLv(ne2,ne1)=Voigt2*Cxe;
  KLv(ne3,ne1)=Voigt3*Cye;
  KLv(ne1,ne2)=Voigt2*Cye;
  KLv(ne2,ne2)=Voigt1*Cye;
  KLv(ne3,ne2)=Voigt3*Cxe;
elseif nsd==3
  KLv(ne1,ne1)=Voigt1*Cxe;
  KLv(ne2,ne1)=Voigt2*Cxe;
  KLv(ne3,ne1)=Voigt2*Cxe;
  KLv(ne4,ne1)=Voigt3*Cye;
  KLv(ne5,ne1)=Voigt3*Cze;
  KLv(ne1,ne2)=Voigt2*Cye;
  KLv(ne2,ne2)=Voigt1*Cye;
  KLv(ne3,ne2)=Voigt2*Cye;
  KLv(ne4,ne2)=Voigt3*Cxe;
  KLv(ne6,ne2)=Voigt3*Cze;
  KLv(ne1,ne3)=Voigt2*Cze;
  KLv(ne2,ne3)=Voigt2*Cze;
  KLv(ne3,ne3)=Voigt1*Cze;
  KLv(ne5,ne3)=Voigt3*Cxe;
  KLv(ne6,ne3)=Voigt3*Cye;
end

if isTimeDependent
  Kvv(ne1,ne1)=alpha(1)/dt*NweT*((reg).*Ne);
  Kvv(ne2,ne2)=alpha(1)/dt*NweT*((reg).*Ne);
  if nsd==3
    Kvv(ne3,ne3)=alpha(1)/dt*NweT*((reg).*Ne);
  end
  
  Kvv(ne1,ne1)=Kvv(ne1,ne1)+NweT*((drdpeg.*dpdteg).*Ne);
  Kvv(ne2,ne2)=Kvv(ne2,ne2)+NweT*((drdpeg.*dpdteg).*Ne);
  if nsd==3
    Kvv(ne3,ne3)=Kvv(ne3,ne3)+NweT*((drdpeg.*dpdteg).*Ne);
  end
  
  Kvp(ne1,ne1)=NweT*((drdpeg.*dvxdteg).*Ne);
  Kvp(ne2,ne1)=NweT*((drdpeg.*dvydteg).*Ne);
  if nsd==3
    Kvp(ne3,ne1)=NweT*((drdpeg.*dvzdteg).*Ne);
  end
  
  Kvp(ne1,ne1)=Kvp(ne1,ne1)+alpha(1)/dt*NweT*((drdpeg.*vxeg).*Ne);
  Kvp(ne2,ne1)=Kvp(ne2,ne1)+alpha(1)/dt*NweT*((drdpeg.*vyeg).*Ne);
  if nsd==3
    Kvp(ne3,ne1)=Kvp(ne3,ne1)+alpha(1)/dt*NweT*((drdpeg.*vzeg).*Ne);
  end
end

if isArbitraryLagrangianEulerian
  Kvv(ne1,ne1)=Kvv(ne1,ne1)+NweT*((reg.*divaeg).*Ne);
  Kvv(ne2,ne2)=Kvv(ne2,ne2)+NweT*((reg.*divaeg).*Ne);
  if nsd==3
    Kvv(ne3,ne3)=Kvv(ne3,ne3)+NweT*((reg.*divaeg).*Ne);
  end
  
  Kvp(ne1,ne1)=Kvp(ne1,ne1)+NweT*((drdpeg.*vxeg.*divaeg).*Ne);
  Kvp(ne2,ne1)=Kvp(ne2,ne1)+NweT*((drdpeg.*vyeg.*divaeg).*Ne);
  if nsd==3
    Kvp(ne3,ne1)=Kvp(ne3,ne1)+NweT*((drdpeg.*vzeg.*divaeg).*Ne);
  end
end

if isConvectiveFlow
  if nsd==2
    Kvv(ne1,ne1)=Kvv(ne1,ne1)-NwexT*((reg.*(2*vxeg-axeg)).*Ne)...
                             -NweyT*((reg.*(  vyeg-ayeg)).*Ne);
    Kvv(ne2,ne1)=Kvv(ne2,ne1)-NwexT*((reg.*vyeg).*Ne);
    Kvv(ne1,ne2)=Kvv(ne1,ne2)-NweyT*((reg.*vxeg).*Ne);
    Kvv(ne2,ne2)=Kvv(ne2,ne2)-NwexT*((reg.*(  vxeg-axeg)).*Ne)...
                             -NweyT*((reg.*(2*vyeg-ayeg)).*Ne);
  elseif nsd==3
    Kvv(ne1,ne1)=Kvv(ne1,ne1)-NwexT*((reg.*(2*vxeg-axeg)).*Ne)...
                             -NweyT*((reg.*(  vyeg-ayeg)).*Ne)...
                             -NwezT*((reg.*(  vzeg-azeg)).*Ne);
    Kvv(ne2,ne1)=Kvv(ne2,ne1)-NwexT*((reg.*vyeg).*Ne);
    Kvv(ne3,ne1)=Kvv(ne3,ne1)-NwezT*((reg.*vzeg).*Ne);
    Kvv(ne1,ne2)=Kvv(ne1,ne2)-NweyT*((reg.*vxeg).*Ne);
    Kvv(ne2,ne2)=Kvv(ne2,ne2)-NwexT*((reg.*(  vxeg-axeg)).*Ne)...
                             -NweyT*((reg.*(2*vyeg-ayeg)).*Ne)...
                             -NwezT*((reg.*(  vzeg-azeg)).*Ne);
    Kvv(ne3,ne2)=Kvv(ne3,ne2)-NwezT*((reg.*vzeg).*Ne);
    Kvv(ne1,ne3)=Kvv(ne1,ne3)-NwezT*((reg.*vxeg).*Ne);
    Kvv(ne2,ne3)=Kvv(ne2,ne3)-NwezT*((reg.*vyeg).*Ne);
    Kvv(ne3,ne3)=Kvv(ne3,ne3)-NwexT*((reg.*(  vxeg-axeg)).*Ne)...
                             -NweyT*((reg.*(  vyeg-ayeg)).*Ne)...
                             -NwezT*((reg.*(2*vzeg-azeg)).*Ne);
  end
  
  if nsd==2
    Kvp(ne1,ne1)=Kvp(ne1,ne1)-NwexT*((drdpeg.*vxeg.*(vxeg-axeg)).*Ne)...
                             -NweyT*((drdpeg.*vxeg.*(vyeg-ayeg)).*Ne);
    Kvp(ne2,ne1)=Kvp(ne2,ne1)-NwexT*((drdpeg.*vyeg.*(vxeg-axeg)).*Ne)...
                             -NweyT*((drdpeg.*vyeg.*(vyeg-ayeg)).*Ne);
  elseif nsd==3
    Kvp(ne1,ne1)=Kvp(ne1,ne1)-NwexT*((drdpeg.*vxeg.*(vxeg-axeg)).*Ne)...
                             -NweyT*((drdpeg.*vxeg.*(vyeg-ayeg)).*Ne)...
                             -NwezT*((drdpeg.*vxeg.*(vzeg-azeg)).*Ne);
    Kvp(ne2,ne1)=Kvp(ne2,ne1)-NwexT*((drdpeg.*vyeg.*(vxeg-axeg)).*Ne)...
                             -NweyT*((drdpeg.*vyeg.*(vyeg-ayeg)).*Ne)...
                             -NwezT*((drdpeg.*vyeg.*(vzeg-azeg)).*Ne);
    Kvp(ne3,ne1)=Kvp(ne3,ne1)-NwexT*((drdpeg.*vzeg.*(vxeg-axeg)).*Ne)...
                             -NweyT*((drdpeg.*vzeg.*(vyeg-ayeg)).*Ne)...
                             -NwezT*((drdpeg.*vzeg.*(vzeg-azeg)).*Ne);
  end
end

Kvp(ne1,ne1)=Kvp(ne1,ne1)+CxeT;
Kvp(ne2,ne1)=Kvp(ne2,ne1)+CyeT;
if nsd==3
  Kvp(ne3,ne1)=Kvp(ne3,ne1)+CzeT;
end

if isTimeDependent
  Kpp(ne1,ne1)=alpha(1)/dt*NweT*((drdpeg).*Ne);
end

if isArbitraryLagrangianEulerian
  Kpp(ne1,ne1)=Kpp(ne1,ne1)+NweT*((drdpeg.*divaeg).*Ne);
end

Kpv(ne1,ne1)=-NwexT*((reg).*Ne);
Kpv(ne1,ne2)=-NweyT*((reg).*Ne);
if nsd==3
  Kpv(ne1,ne3)=-NwezT*((reg).*Ne);
end

Kpp(ne1,ne1)=Kpp(ne1,ne1)-NwexT*((drdpeg.*(vxeg-axeg)).*Ne)...
                         -NweyT*((drdpeg.*(vyeg-ayeg)).*Ne);
if nsd==3
  Kpp(ne1,ne1)=Kpp(ne1,ne1)-NwezT*((drdpeg.*(vzeg-azeg)).*Ne);
end

% Compute rhs
if nsd==2
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(Voigt1*vxeg)...
            -NweyT*(Voigt2*vyeg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(Voigt2*vxeg)...
            -NweyT*(Voigt1*vyeg);
  fL(ne3,1)=+NweT*(Lxyeg)...
            -NweyT*(Voigt3*vxeg)...
            -NwexT*(Voigt3*vyeg);
elseif nsd==3
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(Voigt1*vxeg)...
            -NweyT*(Voigt2*vyeg)...
            -NwezT*(Voigt2*vzeg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(Voigt2*vxeg)...
            -NweyT*(Voigt1*vyeg)...
            -NwezT*(Voigt2*vzeg);
  fL(ne3,1)=+NweT*(Lzzeg)...
            -NwexT*(Voigt2*vxeg)...
            -NweyT*(Voigt2*vyeg)...
            -NwezT*(Voigt1*vzeg);
  fL(ne4,1)=+NweT*(Lxyeg)...
            -NwexT*(Voigt3*vyeg)...
            -NweyT*(Voigt3*vxeg);
  fL(ne5,1)=+NweT*(Lxzeg)...
            -NwexT*(Voigt3*vzeg)...
            -NwezT*(Voigt3*vxeg);
  fL(ne6,1)=+NweT*(Lyzeg)...
            -NweyT*(Voigt3*vzeg)...
            -NwezT*(Voigt3*vyeg);
end

if isTimeDependent
  fv(ne1,1)=-NweT*(reg.*(1/dt*vxeg*alpha(1)...
                        +1/dt*voldxeg*alpha(2:BDFo+1,1))...
               +drdpeg.*(1/dt*peg*alpha(1)...
                        +1/dt*poldeg*alpha(2:BDFo+1,1)).*vxeg);
  fv(ne2,1)=-NweT*(reg.*(1/dt*vyeg*alpha(1)...
                        +1/dt*voldyeg*alpha(2:BDFo+1,1))...
               +drdpeg.*(1/dt*peg*alpha(1)...
                        +1/dt*poldeg*alpha(2:BDFo+1,1)).*vyeg);
  if nsd==3
    fv(ne3,1)=-NweT*(reg.*(1/dt*vzeg*alpha(1)...
                          +1/dt*voldzeg*alpha(2:BDFo+1,1))...
                 +drdpeg.*(1/dt*peg*alpha(1)...
                          +1/dt*poldeg*alpha(2:BDFo+1,1)).*vzeg);
  end
end

if isArbitraryLagrangianEulerian
  fv(ne1,1)=fv(ne1,1)-NweT*(reg.*vxeg.*divaeg);
  fv(ne2,1)=fv(ne2,1)-NweT*(reg.*vyeg.*divaeg);
  if nsd==3
    fv(ne3,1)=fv(ne3,1)-NweT*(reg.*vzeg.*divaeg);
  end
end

if isConvectiveFlow
  if nsd==2
    fv(ne1,1)=fv(ne1,1)+NwexT*(reg.*vxeg.*(vxeg-axeg))...
                       +NweyT*(reg.*vxeg.*(vyeg-ayeg));
    fv(ne2,1)=fv(ne2,1)+NwexT*(reg.*vyeg.*(vxeg-axeg))...
                       +NweyT*(reg.*vyeg.*(vyeg-ayeg));
  elseif nsd==3
    fv(ne1,1)=fv(ne1,1)+NwexT*(reg.*vxeg.*(vxeg-axeg))...
                       +NweyT*(reg.*vxeg.*(vyeg-ayeg))...
                       +NwezT*(reg.*vxeg.*(vzeg-azeg));
    fv(ne2,1)=fv(ne2,1)+NwexT*(reg.*vyeg.*(vxeg-axeg))...
                       +NweyT*(reg.*vyeg.*(vyeg-ayeg))...
                       +NwezT*(reg.*vyeg.*(vzeg-azeg));
    fv(ne3,1)=fv(ne3,1)+NwexT*(reg.*vzeg.*(vxeg-axeg))...
                       +NweyT*(reg.*vzeg.*(vyeg-ayeg))...
                       +NwezT*(reg.*vzeg.*(vzeg-azeg));
  end
end

if nsd==2
  fv(ne1,1)=fv(ne1,1)-NweT*(+Voigt1*(Nex*Lxxe)+Voigt2*(Nex*Lyye)...
                            +Voigt3*(Ney*Lxye)...
                            +(Nex*pe));
  fv(ne2,1)=fv(ne2,1)-NweT*(+Voigt2*(Ney*Lxxe)+Voigt1*(Ney*Lyye)...
                            +Voigt3*(Nex*Lxye)...
                            +(Ney*pe));
elseif nsd==3
  fv(ne1,1)=fv(ne1,1)-NweT*(+Voigt1*(Nex*Lxxe)+Voigt2*(Nex*Lyye)+Voigt2*(Nex*Lzze)...
                            +Voigt3*(Ney*Lxye)+Voigt3*(Nez*Lxze)...
                            +(Nex*pe));
  fv(ne2,1)=fv(ne2,1)-NweT*(+Voigt2*(Ney*Lxxe)+Voigt1*(Ney*Lyye)+Voigt2*(Ney*Lzze)...
                            +Voigt3*(Nex*Lxye)+Voigt3*(Nez*Lyze)...
                            +(Ney*pe));
  fv(ne3,1)=fv(ne3,1)-NweT*(+Voigt2*(Nez*Lxxe)+Voigt2*(Nez*Lyye)+Voigt1*(Nez*Lzze)...
                            +Voigt3*(Nex*Lxze)+Voigt3*(Ney*Lyze)...
                            +(Nez*pe));
end

fv(ne1,1)=fv(ne1,1)+NweT*(fxeg);
fv(ne2,1)=fv(ne2,1)+NweT*(fyeg);
if nsd==3
  fv(ne3,1)=fv(ne3,1)+NweT*(fzeg);
end

if isTimeDependent
  fp(ne1,1)=-NweT*(drdpeg.*(1/dt*peg*alpha(1)...
                           +1/dt*poldeg*alpha(2:BDFo+1,1)));
end

if isArbitraryLagrangianEulerian
  fp(ne1,1)=fp(ne1,1)-NweT*(reg.*divaeg);
end

fp(ne1,1)=fp(ne1,1)+NwexT*(reg.*(vxeg-axeg))...
                   +NweyT*(reg.*(vyeg-ayeg));
if nsd==3
  fp(ne1,1)=fp(ne1,1)+NwezT*(reg.*(vzeg-azeg));
end

fp(ne1,1)=fp(ne1,1)+NweT*(Rceg);

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
    isExterior=Faces.Exterior(iFace);
    isDirichlet_v_x=Faces.Dirichlet_v_x(iFace);
    isDirichlet_v_y=Faces.Dirichlet_v_y(iFace);
    if nsd==3; isDirichlet_v_z=Faces.Dirichlet_v_z(iFace); end
    isDirichlet_p=Faces.Dirichlet_p(iFace);
    isNeumann_t_x=Faces.Neumann_t_x(iFace);
    isNeumann_t_y=Faces.Neumann_t_y(iFace);
    if nsd==3; isNeumann_t_z=Faces.Neumann_t_z(iFace); end
    isDirichlet_v=isDirichlet_v_x || isDirichlet_v_y || (nsd==3 && isDirichlet_v_z);
    isNeumann_t=isNeumann_t_x || isNeumann_t_y || (nsd==3 && isNeumann_t_z);
    
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
    nefV1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefV2=nefV1+NumFaceNodes;
    nefV3=nefV2+NumFaceNodes;
    nefP1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters(iD1).Degree,Node2Match1stNode1);
      nefU1=nefU1(order);
      nefU2=nefU2(order);
      nefU3=nefU3(order);
      nefU4=nefU4(order);
      nefV1=nefV1(order);
      nefV2=nefV2(order);
      nefV3=nefV3(order);
      nefP1=nefP1(order);
    end
    
    % Compute variables at nodes
    if nsd==2
      Lxxf=Lxxe(nf1);
      Lyyf=Lyye(nf1);
      Lxyf=Lxye(nf1);
    elseif nsd==3
      Lxxf=Lxxe(nf1);
      Lyyf=Lyye(nf1);
      Lzzf=Lzze(nf1);
      Lxyf=Lxye(nf1);
      Lxzf=Lxze(nf1);
      Lyzf=Lyze(nf1);
    end
    vxf=vxe(nf1);
    vyf=vye(nf1);
    if nsd==3
      vzf=vze(nf1);
    end
    pf=pe(nf1);
    axf=axe(nf1);
    ayf=aye(nf1);
    if nsd==3
      azf=aze(nf1);
    end
    Vxf=Ue(nefU1);
    Vyf=Ue(nefU2);
    if nsd==3
      Vzf=Ue(nefU3);
    end
    if nsd==2
      Pf=Ue(nefU3);
    elseif nsd==3
      Pf=Ue(nefU4);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    if nsd==2
      Lxxfg=Nf*Lxxf;
      Lyyfg=Nf*Lyyf;
      Lxyfg=Nf*Lxyf;
    elseif nsd==3
      Lxxfg=Nf*Lxxf;
      Lyyfg=Nf*Lyyf;
      Lzzfg=Nf*Lzzf;
      Lxyfg=Nf*Lxyf;
      Lxzfg=Nf*Lxzf;
      Lyzfg=Nf*Lyzf;
    end
    vxfg=Nf*vxf;
    vyfg=Nf*vyf;
    if nsd==3
      vzfg=Nf*vzf;
    end
    pfg=Nf*pf;
    axfg=Nf*axf;
    ayfg=Nf*ayf;
    if nsd==3
      azfg=Nf*azf;
    end
    Vxfg=Nf*Vxf;
    Vyfg=Nf*Vyf;
    if nsd==3
      Vzfg=Nf*Vzf;
    end
    Pfg=Nf*Pf;
    if isDirichlet_v
      vDfg=vD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      vDxfg=vDfg(:,1);
      vDyfg=vDfg(:,2);
      if nsd==3
        vDzfg=vDfg(:,3);
      end
    end
    if isDirichlet_p
      pDfg=pD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isNeumann_t
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
    
    % Compute common terms
    if isDirichlet_v_x
      Vxfg=vDxfg;
    end
    if isDirichlet_v_y
      Vyfg=vDyfg;
    end
    if nsd==3 && isDirichlet_v_z
      Vzfg=vDzfg;
    end
    if isDirichlet_p
      Pfg=pDfg;
    end
    rfg=r(Pfg);
    drdPfg=drdp(Pfg);
    
    % Compute lhs
    Kvv(nf1,nf1)=Kvv(nf1,nf1)+tauV*Mf;
    Kvv(nf2,nf2)=Kvv(nf2,nf2)+tauV*Mf;
    if nsd==3
      Kvv(nf3,nf3)=Kvv(nf3,nf3)+tauV*Mf;
    end
    
    Kpp(nf1,nf1)=Kpp(nf1,nf1)+tauP*Mf;
    
    if not(isDirichlet_v_x)
      if nsd==2
        KLV(nf1,nefV1)=KLV(nf1,nefV1)-Voigt1*Mfnx;
        KLV(nf2,nefV1)=KLV(nf2,nefV1)-Voigt2*Mfnx;
        KLV(nf3,nefV1)=KLV(nf3,nefV1)-Voigt3*Mfny;
      elseif nsd==3
        KLV(nf1,nefV1)=KLV(nf1,nefV1)-Voigt1*Mfnx;
        KLV(nf2,nefV1)=KLV(nf2,nefV1)-Voigt2*Mfnx;
        KLV(nf3,nefV1)=KLV(nf3,nefV1)-Voigt2*Mfnx;
        KLV(nf4,nefV1)=KLV(nf4,nefV1)-Voigt3*Mfny;
        KLV(nf5,nefV1)=KLV(nf5,nefV1)-Voigt3*Mfnz;
      end
    end
    
    if not(isDirichlet_v_y)
      if nsd==2
        KLV(nf1,nefV2)=KLV(nf1,nefV2)-Voigt2*Mfny;
        KLV(nf2,nefV2)=KLV(nf2,nefV2)-Voigt1*Mfny;
        KLV(nf3,nefV2)=KLV(nf3,nefV2)-Voigt3*Mfnx;
      elseif nsd==3
        KLV(nf1,nefV2)=KLV(nf1,nefV2)-Voigt2*Mfny;
        KLV(nf2,nefV2)=KLV(nf2,nefV2)-Voigt1*Mfny;
        KLV(nf3,nefV2)=KLV(nf3,nefV2)-Voigt2*Mfny;
        KLV(nf4,nefV2)=KLV(nf4,nefV2)-Voigt3*Mfnx;
        KLV(nf6,nefV2)=KLV(nf6,nefV2)-Voigt3*Mfnz;
      end
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      KLV(nf1,nefV3)=KLV(nf1,nefV3)-Voigt2*Mfnz;
      KLV(nf2,nefV3)=KLV(nf2,nefV3)-Voigt2*Mfnz;
      KLV(nf3,nefV3)=KLV(nf3,nefV3)-Voigt1*Mfnz;
      KLV(nf5,nefV3)=KLV(nf5,nefV3)-Voigt3*Mfnx;
      KLV(nf6,nefV3)=KLV(nf6,nefV3)-Voigt3*Mfny;
    end
    
    if not(isDirichlet_v_x) && isConvectiveFlow
      if nsd==2
        KvV(nf1,nefV1)=KvV(nf1,nefV1)+NwfT*((rfg.*((2*Vxfg-axfg).*nx...
                                                  +(  Vyfg-ayfg).*ny)).*Nf);
        KvV(nf2,nefV1)=KvV(nf2,nefV1)+NwfT*((rfg.*(Vyfg.*nx)).*Nf);
      elseif nsd==3
        KvV(nf1,nefV1)=KvV(nf1,nefV1)+NwfT*((rfg.*((2*Vxfg-axfg).*nx...
                                                  +(  Vyfg-ayfg).*ny...
                                                  +(  Vzfg-azfg).*nz)).*Nf);
        KvV(nf2,nefV1)=KvV(nf2,nefV1)+NwfT*((rfg.*(Vyfg.*nx)).*Nf);
        KvV(nf3,nefV1)=KvV(nf3,nefV1)+NwfT*((rfg.*(Vzfg.*nx)).*Nf);
      end
    end
    
    if not(isDirichlet_v_y) && isConvectiveFlow
      if nsd==2
        KvV(nf1,nefV2)=KvV(nf1,nefV2)+NwfT*((rfg.*(Vxfg.*ny)).*Nf);
        KvV(nf2,nefV2)=KvV(nf2,nefV2)+NwfT*((rfg.*((  Vxfg-axfg).*nx...
                                                  +(2*Vyfg-ayfg).*ny)).*Nf);
      elseif nsd==3
        KvV(nf1,nefV2)=KvV(nf1,nefV2)+NwfT*((rfg.*(Vxfg.*ny)).*Nf);
        KvV(nf2,nefV2)=KvV(nf2,nefV2)+NwfT*((rfg.*((  Vxfg-axfg).*nx...
                                                  +(2*Vyfg-ayfg).*ny...
                                                  +(  Vzfg-azfg).*nz)).*Nf);
        KvV(nf3,nefV2)=KvV(nf3,nefV2)+NwfT*((rfg.*(Vzfg.*ny)).*Nf);
      end
    end
    
    if nsd==3 && not(isDirichlet_v_z) && isConvectiveFlow
      KvV(nf1,nefV3)=KvV(nf1,nefV3)+NwfT*((rfg.*(Vxfg.*nz)).*Nf);
      KvV(nf2,nefV3)=KvV(nf2,nefV3)+NwfT*((rfg.*(Vyfg.*nz)).*Nf);
      KvV(nf3,nefV3)=KvV(nf3,nefV3)+NwfT*((rfg.*((  Vxfg-axfg).*nx...
                                                +(  Vyfg-ayfg).*ny...
                                                +(2*Vzfg-azfg).*nz)).*Nf);
    end
        
    if not(isDirichlet_p) && isConvectiveFlow
      if nsd==2
        KvP(nf1,nefP1)=KvP(nf1,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny).*Vxfg).*Nf);
        KvP(nf2,nefP1)=KvP(nf2,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny).*Vyfg).*Nf);
      elseif nsd==3
        KvP(nf1,nefP1)=KvP(nf1,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny...
                                                     +(Vzfg-azfg).*nz).*Vxfg).*Nf);
        KvP(nf2,nefP1)=KvP(nf2,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny...
                                                     +(Vzfg-azfg).*nz).*Vyfg).*Nf);
        KvP(nf3,nefP1)=KvP(nf3,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny...
                                                     +(Vzfg-azfg).*nz).*Vzfg).*Nf);
      end
    end
    
    if not(isDirichlet_v_x)
      KvV(nf1,nefV1)=KvV(nf1,nefV1)-tauV*Mf;
    end
    
    if not(isDirichlet_v_y)
      KvV(nf2,nefV2)=KvV(nf2,nefV2)-tauV*Mf;
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      KvV(nf3,nefV3)=KvV(nf3,nefV3)-tauV*Mf;
    end
    
    if not(isDirichlet_v_x)
      KpV(nf1,nefV1)=KpV(nf1,nefV1)+NwfT*((rfg.*nx).*Nf);
    end
    
    if not(isDirichlet_v_y)
      KpV(nf1,nefV2)=KpV(nf1,nefV2)+NwfT*((rfg.*ny).*Nf);
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      KpV(nf1,nefV3)=KpV(nf1,nefV3)+NwfT*((rfg.*nz).*Nf);
    end
    
    if not(isDirichlet_p)
      if nsd==2
        KpP(nf1,nefP1)=KpP(nf1,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny)).*Nf);
      elseif nsd==3
        KpP(nf1,nefP1)=KpP(nf1,nefP1)+NwfT*((drdPfg.*((Vxfg-axfg).*nx...
                                                     +(Vyfg-ayfg).*ny...
                                                     +(Vzfg-azfg).*nz)).*Nf);
      end
    end
    
    if not(isDirichlet_p)
      KpP(nf1,nefP1)=KpP(nf1,nefP1)-tauP*Mf;
    end
    
    if not(isExterior) || isNeumann_t_x
      if nsd==2
        KVL(nefV1,nf1)=KVL(nefV1,nf1)-Voigt1*Mfnx;
        KVL(nefV1,nf2)=KVL(nefV1,nf2)-Voigt2*Mfnx;
        KVL(nefV1,nf3)=KVL(nefV1,nf3)-Voigt3*Mfny;
      elseif nsd==3
        KVL(nefV1,nf1)=KVL(nefV1,nf1)-Voigt1*Mfnx;
        KVL(nefV1,nf2)=KVL(nefV1,nf2)-Voigt2*Mfnx;
        KVL(nefV1,nf3)=KVL(nefV1,nf3)-Voigt2*Mfnx;
        KVL(nefV1,nf4)=KVL(nefV1,nf4)-Voigt3*Mfny;
        KVL(nefV1,nf5)=KVL(nefV1,nf5)-Voigt3*Mfnz;
      end
    end
    
    if not(isExterior) || isNeumann_t_y
      if nsd==2
        KVL(nefV2,nf1)=KVL(nefV2,nf1)-Voigt2*Mfny;
        KVL(nefV2,nf2)=KVL(nefV2,nf2)-Voigt1*Mfny;
        KVL(nefV2,nf3)=KVL(nefV2,nf3)-Voigt3*Mfnx;
      elseif nsd==3
        KVL(nefV2,nf1)=KVL(nefV2,nf1)-Voigt2*Mfny;
        KVL(nefV2,nf2)=KVL(nefV2,nf2)-Voigt1*Mfny;
        KVL(nefV2,nf3)=KVL(nefV2,nf3)-Voigt2*Mfny;
        KVL(nefV2,nf4)=KVL(nefV2,nf4)-Voigt3*Mfnx;
        KVL(nefV2,nf6)=KVL(nefV2,nf6)-Voigt3*Mfnz;
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z)
      KVL(nefV3,nf1)=KVL(nefV3,nf1)-Voigt2*Mfnz;
      KVL(nefV3,nf2)=KVL(nefV3,nf2)-Voigt2*Mfnz;
      KVL(nefV3,nf3)=KVL(nefV3,nf3)-Voigt1*Mfnz;
      KVL(nefV3,nf5)=KVL(nefV3,nf5)-Voigt3*Mfnx;
      KVL(nefV3,nf6)=KVL(nefV3,nf6)-Voigt3*Mfny;
    end
    
    if not(isDirichlet_v_x)
      KVv(nefV1,nf1)=KVv(nefV1,nf1)-tauV*Mf;
    end
    
    if not(isDirichlet_v_y)
      KVv(nefV2,nf2)=KVv(nefV2,nf2)-tauV*Mf;
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      KVv(nefV3,nf3)=KVv(nefV3,nf3)-tauV*Mf;
    end
    
    if not(isDirichlet_v_x)
      KVV(nefV1,nefV1)=KVV(nefV1,nefV1)+tauV*Mf;
    end
    
    if not(isDirichlet_v_y)
      KVV(nefV2,nefV2)=KVV(nefV2,nefV2)+tauV*Mf;
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      KVV(nefV3,nefV3)=KVV(nefV3,nefV3)+tauV*Mf;
    end
    
    if (not(isExterior) || isNeumann_t_x) && not(isDirichlet_p)
      KVP(nefV1,nefP1)=KVP(nefV1,nefP1)-Mfnx;
    end
    
    if (not(isExterior) || isNeumann_t_y) && not(isDirichlet_p)
      KVP(nefV2,nefP1)=KVP(nefV2,nefP1)-Mfny;
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z) && not(isDirichlet_p)
      KVP(nefV3,nefP1)=KVP(nefV3,nefP1)-Mfnz;
    end
    
    if not(isDirichlet_p)
      KPp(nefP1,nf1)=KPp(nefP1,nf1)+tauP*Mf;
    end
    
    if not(isDirichlet_p)
      KPP(nefP1,nefP1)=KPP(nefP1,nefP1)-tauP*Mf;
    end
    
    % Compute rhs
    fv(nf1,1)=fv(nf1,1)-NwfT*(tauV*vxfg);
    fv(nf2,1)=fv(nf2,1)-NwfT*(tauV*vyfg);
    if nsd==3
      fv(nf3,1)=fv(nf3,1)-NwfT*(tauV*vzfg);
    end
    
    fp(nf1,1)=fp(nf1,1)-NwfT*(tauP*pfg);
    
    if nsd==2
      fL(nf1,1)=fL(nf1,1)+NwfT*(+Voigt1*nx.*Vxfg...
                                +Voigt2*ny.*Vyfg);
      fL(nf2,1)=fL(nf2,1)+NwfT*(+Voigt2*nx.*Vxfg...
                                +Voigt1*ny.*Vyfg);
      fL(nf3,1)=fL(nf3,1)+NwfT*(+Voigt3*ny.*Vxfg...
                                +Voigt3*nx.*Vyfg);
    elseif nsd==3
      fL(nf1,1)=fL(nf1,1)+NwfT*(+Voigt1*nx.*Vxfg...
                                +Voigt2*ny.*Vyfg...
                                +Voigt2*nz.*Vzfg);
      fL(nf2,1)=fL(nf2,1)+NwfT*(+Voigt2*nx.*Vxfg...
                                +Voigt1*ny.*Vyfg...
                                +Voigt2*nz.*Vzfg);
      fL(nf3,1)=fL(nf3,1)+NwfT*(+Voigt2*nx.*Vxfg...
                                +Voigt2*ny.*Vyfg...
                                +Voigt1*nz.*Vzfg);
      fL(nf4,1)=fL(nf4,1)+NwfT*(+Voigt3*nx.*Vyfg...
                                +Voigt3*ny.*Vxfg);
      fL(nf5,1)=fL(nf5,1)+NwfT*(+Voigt3*nx.*Vzfg...
                                +Voigt3*nz.*Vxfg);
      fL(nf6,1)=fL(nf6,1)+NwfT*(+Voigt3*ny.*Vzfg...
                                +Voigt3*nz.*Vyfg);
    end
    
    if isConvectiveFlow
      if nsd==2
        fv(nf1,1)=fv(nf1,1)-NwfT*(rfg.*Vxfg.*((Vxfg-axfg).*nx+(Vyfg-ayfg).*ny));
        fv(nf2,1)=fv(nf2,1)-NwfT*(rfg.*Vyfg.*((Vxfg-axfg).*nx+(Vyfg-ayfg).*ny));
      elseif nsd==3
        fv(nf1,1)=fv(nf1,1)-NwfT*(rfg.*Vxfg.*((Vxfg-axfg).*nx+(Vyfg-ayfg).*ny+(Vzfg-azfg).*nz));
        fv(nf2,1)=fv(nf2,1)-NwfT*(rfg.*Vyfg.*((Vxfg-axfg).*nx+(Vyfg-ayfg).*ny+(Vzfg-azfg).*nz));
        fv(nf3,1)=fv(nf3,1)-NwfT*(rfg.*Vzfg.*((Vxfg-axfg).*nx+(Vyfg-ayfg).*ny+(Vzfg-azfg).*nz));
      end
    end
    
    fv(nf1,1)=fv(nf1,1)+NwfT*(tauV*Vxfg);
    fv(nf2,1)=fv(nf2,1)+NwfT*(tauV*Vyfg);
    if nsd==3
      fv(nf3,1)=fv(nf3,1)+NwfT*(tauV*Vzfg);
    end
    
    fp(nf1,1)=fp(nf1,1)-NwfT*(rfg.*((Vxfg-axfg).*nx+(Vyfg-ayfg).*ny)-tauP*Pfg);
    if nsd==3
      fp(nf1,1)=fp(nf1,1)-NwfT*(rfg.*(Vzfg-azfg).*nz);
    end
    
    if not(isExterior) || isNeumann_t_x
      if nsd==2
        fV(nefV1,1)=fV(nefV1,1)+NwfT*(+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg...
                                      +Voigt3*ny.*Lxyfg...
                                      +Pfg.*nx);
      elseif nsd==3
        fV(nefV1,1)=fV(nefV1,1)+NwfT*(+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg+Voigt2*nx.*Lzzfg...
                                      +Voigt3*ny.*Lxyfg+Voigt3*nz.*Lxzfg...
                                      +Pfg.*nx);
      end
    end
    
    if not(isExterior) || isNeumann_t_y
      if nsd==2
        fV(nefV2,1)=fV(nefV2,1)+NwfT*(+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg...
                                      +Voigt3*nx.*Lxyfg...
                                      +Pfg.*ny);
      elseif nsd==3
        fV(nefV2,1)=fV(nefV2,1)+NwfT*(+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg+Voigt2*ny.*Lzzfg...
                                      +Voigt3*nx.*Lxyfg+Voigt3*nz.*Lyzfg...
                                      +Pfg.*ny);
      end
    end
    
    if nsd==3 && (not(isExterior) || isNeumann_t_z)
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(+Voigt2*nz.*Lxxfg+Voigt2*nz.*Lyyfg+Voigt1*nz.*Lzzfg...
                                    +Voigt3*nx.*Lxzfg+Voigt3*ny.*Lyzfg...
                                    +Pfg.*nz);
    end
    
    if not(isDirichlet_v_x)
      fV(nefV1,1)=fV(nefV1,1)+NwfT*(+tauV*(vxfg-Vxfg));
    end
    
    if not(isDirichlet_v_y)
      fV(nefV2,1)=fV(nefV2,1)+NwfT*(+tauV*(vyfg-Vyfg));
    end
    
    if nsd==3 && not(isDirichlet_v_z)
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(+tauV*(vzfg-Vzfg));
    end
    
    if not(isDirichlet_p)
      fP(nefP1,1)=fP(nefP1,1)-NwfT*(tauP*(pfg-Pfg));
    end
    
    if isNeumann_t_x
      fV(nefV1,1)=fV(nefV1,1)+NwfT*(tNxfg);
    end
    if isNeumann_t_y
      fV(nefV2,1)=fV(nefV2,1)+NwfT*(tNyfg);
    end
    if nsd==3 && isNeumann_t_z
      fV(nefV3,1)=fV(nefV3,1)+NwfT*(tNzfg);
    end
    
    % Remove undetermination
    if isDirichlet_v_x
      KVV(nefV1,nefV1)=eye(NumFaceNodes);
    end
    if isDirichlet_v_y
      KVV(nefV2,nefV2)=eye(NumFaceNodes);
    end
    if nsd==3 && isDirichlet_v_z
      KVV(nefV3,nefV3)=eye(NumFaceNodes);
    end
    if isDirichlet_p
      KPP(nefP1,nefP1)=eye(NumFaceNodes);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iL=1:msd*NumElementNodes;
iv=iL(end)+(1:nsd*NumElementNodes);
ip=iv(end)+(1:NumElementNodes);
iV=reshape((0:NumElementFaces-1)*(nsd+1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iP=reshape((0:NumElementFaces-1)*(nsd+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+nsd*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+nsd+1)*NumElementNodes,(msd+nsd+1)*NumElementNodes);
LhsLG=zeros((msd+nsd+1)*NumElementNodes,(nsd+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((nsd+1)*NumElementFaces*NumFaceNodes,(msd+nsd+1)*NumElementNodes);
LhsGG=zeros((nsd+1)*NumElementFaces*NumFaceNodes,(nsd+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+nsd+1)*NumElementNodes,1);
RhsG=zeros((nsd+1)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(iL,iL)=KLL;
LhsLL(iL,iv)=KLv;
LhsLL(iv,iL)=KLv';
LhsLL(iv,iv)=Kvv;
LhsLL(iv,ip)=Kvp;
LhsLL(ip,iv)=Kpv;
LhsLL(ip,ip)=Kpp;

% Lhs local-global
LhsLG(iL,iV)=KLV;
LhsLG(iv,iV)=KvV;
LhsLG(iv,iP)=KvP;
LhsLG(ip,iV)=KpV;
LhsLG(ip,iP)=KpP;

% Rhs local
RhsL(iL,1)=fL;
RhsL(iv,1)=fv;
RhsL(ip,1)=fp;

% Lhs global-local
LhsGL(iV,iL)=KVL;
LhsGL(iV,iv)=KVv;
LhsGL(iP,ip)=KPp;

% Lhs global-global
LhsGG(iV,iV)=KVV;
LhsGG(iV,iP)=KVP;
LhsGG(iP,iP)=KPP;

% Rhs global
RhsG(iV,1)=fV;
RhsG(iP,1)=fP;

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
mu=Parameters.DynamicViscosity;
lambda=-2/3*Parameters.DynamicViscosity;
vD=Parameters.Velocity;
Xe=Nodes';
t=Time.Time;

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
ve=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
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
  Lxxe=Le(nle1);
  Lyye=Le(nle2);
  Lxye=Le(nle3);
elseif nsd==3
  Lxxe=Le(nle1);
  Lyye=Le(nle2);
  Lzze=Le(nle3);
  Lxye=Le(nle4);
  Lxze=Le(nle5);
  Lyze=Le(nle6);
end
vxe=ve(nle1);
vye=ve(nle2);
if nsd==3
  vze=ve(nle3);
end

% Compute variables at Gauss points
if nsd==2
  Lxxeg=Nle*Lxxe;
  Lyyeg=Nle*Lyye;
  Lxyeg=Nle*Lxye;
elseif nsd==3
  Lxxeg=Nle*Lxxe;
  Lyyeg=Nle*Lyye;
  Lzzeg=Nle*Lzze;
  Lxyeg=Nle*Lxye;
  Lxzeg=Nle*Lxze;
  Lyzeg=Nle*Lyze;
end
vxeg=Nle*vxe;
vyeg=Nle*vye;
if nsd==3
  vzeg=Nle*vze;
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

% Compute linearization of viscous stress
if nsd==2
  Voigt1=(sqrt(2*(mu+lambda))+sqrt(2*mu))/2;
  Voigt2=(sqrt(2*(mu+lambda))-sqrt(2*mu))/2;
  Voigt3=sqrt(mu);
elseif nsd==3
  Voigt1=(sqrt(2*mu+3*lambda)+2*sqrt(2*mu))/3;
  Voigt2=(sqrt(2*mu+3*lambda)-1*sqrt(2*mu))/3;
  Voigt3=sqrt(mu);
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
  fp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg);
  fp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg);
elseif nsd==3
  fp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg)+NwezT*(Lxzeg);
  fp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg)+NwezT*(Lyzeg);
  fp(ne3,1)=NwexT*(Lxzeg)+NweyT*(Lyzeg)+NwezT*(Lzzeg);
end

ft(1,1)=Nw1eT*(vxeg);
ft(2,1)=Nw1eT*(vyeg);
if nsd==3
  ft(3,1)=Nw1eT*(vzeg);
end

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions('Face',RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  Xfg=Nlf*Xf;
  
  % Check boundary
  isDirichlet_v_x=Faces.Dirichlet_v_x(iFace);
  isDirichlet_v_y=Faces.Dirichlet_v_y(iFace);
  if nsd==3; isDirichlet_v_z=Faces.Dirichlet_v_z(iFace); end
  isDirichlet_v=isDirichlet_v_x || isDirichlet_v_y || (nsd==3 && isDirichlet_v_z);
  
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
  Vxf=Ue(nlefU1);
  Vyf=Ue(nlefU2);
  if nsd==3
    Vzf=Ue(nlefU3);
  end
  
  % Compute variables at Gauss points
  Vxfg=Nlf*Vxf;
  Vyfg=Nlf*Vyf;
  if nsd==3
    Vzfg=Nlf*Vzf;
  end
  if isDirichlet_v
    vDfg=vD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    vDxfg=vDfg(:,1);
    vDyfg=vDfg(:,2);
    if nsd==3
      vDzfg=vDfg(:,3);
    end
  end
  
  % Compute common terms
  if isDirichlet_v_x
    Vxfg=vDxfg;
  end
  if isDirichlet_v_y
    Vyfg=vDyfg;
  end
  if nsd==3 && isDirichlet_v_z
    Vzfg=vDzfg;
  end
  
  % Compute basic matrices
  Nw1fT=(wfg.*N1f)';
  
  % Compute rhs
  if nsd==2
    fr(1)=fr(1)+Nw1fT*(-Vxfg.*ny+Vyfg.*nx);
  elseif nsd==3
    fr(1)=fr(1)+Nw1fT*(-Vyfg.*nz+Vzfg.*ny);
    fr(2)=fr(2)+Nw1fT*(+Vxfg.*nz-Vzfg.*nx);
    fr(3)=fr(3)+Nw1fT*(-Vxfg.*ny+Vyfg.*nx);
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