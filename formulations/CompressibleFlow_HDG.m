classdef CompressibleFlow_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1+NumSpaceDim+1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) NumSpaceDim*(NumSpaceDim+1)/2+NumSpaceDim+1+NumSpaceDim+1;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) NumSpaceDim+1;

    % Discretization type
    DiscretizationType='HDG';

    % Time derivative order
    TimeDerOrder=1;
    
    % Domain
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
      
      % Consider one for the density instead of zero
      for iFace=1:Sizes(iD).NumFaces
        Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp.*Sizes(iD).NumFaceNodes+...
          (1:Sizes(iD).NumFaceNodes))=1;
      end
      Block(iD,iD).SolutionLocal(:,Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1)=1;
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,~,Time,~,Sizes)
      if strcmp(Time.TimeDependent,'yes')
        Block(iD,iD).SolutionLocal=[...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
          zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
          Parameters(iD).Density(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).Momentum(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime),...
          Parameters(iD).Energy(...
          Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.InitialTime)];
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Block(iD,iD).SolutionOld(:,:,1+iBDF)=[...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumVoigtComp),...
              zeros(Sizes(iD).NumLocalNodes,Sizes(iD).NumSpaceDim),...
              Parameters(iD).Density(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).Momentum(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld),...
              Parameters(iD).Energy(...
              Mesh(iD).Nodes(1,:)',Mesh(iD).Nodes(2,:)',Mesh(iD).Nodes(3,:)',Time.TimeOld)];
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
          Parameters,Time,RefElement,Sizes);
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
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,~,RefElement,Sizes)
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
        Results(iD).ScaledStrainRate=[];
        Results(iD).ScaledTemperatureGradient=[];
        Results(iD).Density=[];
        Results(iD).Momentum=[];
        Results(iD).Energy=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).ScaledStrainRate(:,:,iST)=Block(iD,iD).SolutionLocal(:,1:Sizes(iD).NumVoigtComp);
      Results(iD).ScaledTemperatureGradient(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+(1:Sizes(iD).NumSpaceDim));
      Results(iD).Density(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1);
      Results(iD).Momentum(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1+(1:Sizes(iD).NumSpaceDim));
      Results(iD).Energy(:,:,iST)=Block(iD,iD).SolutionLocal(:,...
        Sizes(iD).NumVoigtComp+Sizes(iD).NumSpaceDim+1+Sizes(iD).NumSpaceDim+1);
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).VelocityPost=Block(iD,iD).SolutionPost(:,1:Sizes(iD).NumSpaceDim);
        Results(iD).TemperaturePost=Block(iD,iD).SolutionPost(:,Sizes(iD).NumSpaceDim+1);
      end
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isConvectiveFlow=strcmp(Parameters.ConvectiveFlow,'yes');
SolveEnergyEquation=strcmp(Parameters.SolveEnergyEquation,'yes');
isTimeDependent=strcmp(Time.TimeDependent,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
mu=Parameters.DynamicViscosity;
lambda=-2/3*Parameters.DynamicViscosity;
kappa=Parameters.ThermalConductivity;
p=Parameters.ComputePressure;
dpdu=Parameters.ComputePressureLinearization;
T=Parameters.ComputeTemperature;
dTdu=Parameters.ComputeTemperatureLinearization;
rD=Parameters.Density;
wD=Parameters.Momentum;
eD=Parameters.Energy;
Rc=Parameters.ResidualContinuity;
f=Parameters.Force;
s=Parameters.HeatSource;
Xe=Nodes';
tauR=Parameters.StabDensity;
tauW=Parameters.StabMomentum;
tauE=Parameters.StabEnergy;
t=Time.Time;
if isTimeDependent
  dt=Time.TimeStepSize;
  BDFo=Time.BDFOrderEff;
  alpha=Time.BDF1stDerEff;
end

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
Qe=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
re=reshape(SolutionLocal(:,msd+nsd+1),[],1);
we=reshape(SolutionLocal(:,msd+nsd+1+(1:nsd)),[],1);
ee=reshape(SolutionLocal(:,msd+nsd+1+nsd+1),[],1);
Ue=SolutionGlobal;
if isTimeDependent
  rolde=reshape(SolutionOld(:,msd+nsd+1,:),[],BDFo);
  wolde=reshape(SolutionOld(:,msd+nsd+1+(1:nsd),:),[],BDFo);
  eolde=reshape(SolutionOld(:,msd+nsd+1+nsd+1,:),[],BDFo);
end

% Initialize lhs
KLL=zeros(msd*NumElementNodes,msd*NumElementNodes);
KLr=zeros(msd*NumElementNodes,NumElementNodes);
KLw=zeros(msd*NumElementNodes,nsd*NumElementNodes);
KLR=zeros(msd*NumElementNodes,NumElementFaces*NumFaceNodes);
KLW=zeros(msd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KQQ=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KQr=zeros(nsd*NumElementNodes,NumElementNodes);
KQw=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KQe=zeros(nsd*NumElementNodes,NumElementNodes);
KQR=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
KQW=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KQE=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Krr=zeros(NumElementNodes,NumElementNodes);
Krw=zeros(NumElementNodes,nsd*NumElementNodes);
KrR=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KrW=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KwL=zeros(nsd*NumElementNodes,msd*NumElementNodes);
Kwr=zeros(nsd*NumElementNodes,NumElementNodes);
Kww=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kwe=zeros(nsd*NumElementNodes,NumElementNodes);
KwR=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
KwW=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KeL=zeros(NumElementNodes,msd*NumElementNodes);
KeQ=zeros(NumElementNodes,nsd*NumElementNodes);
Ker=zeros(NumElementNodes,NumElementNodes);
Kew=zeros(NumElementNodes,nsd*NumElementNodes);
Kee=zeros(NumElementNodes,NumElementNodes);
KeR=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KeW=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KeE=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KRr=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KRR=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KRW=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KWL=zeros(nsd*NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KWw=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KWR=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KWW=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KWE=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KEL=zeros(NumElementFaces*NumFaceNodes,msd*NumElementNodes);
KEQ=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KEe=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KER=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KEW=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KEE=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);

% Initialize rhs
fL=zeros(msd*NumElementNodes,1);
fQ=zeros(nsd*NumElementNodes,1);
fr=zeros(NumElementNodes,1);
fw=zeros(nsd*NumElementNodes,1);
fe=zeros(NumElementNodes,1);
fR=zeros(NumElementFaces*NumFaceNodes,1);
fW=zeros(nsd*NumElementFaces*NumFaceNodes,1);
fE=zeros(NumElementFaces*NumFaceNodes,1);

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

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

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
Qxe=Qe(ne1);
Qye=Qe(ne2);
if nsd==3
  Qze=Qe(ne3);
end
wxe=we(ne1);
wye=we(ne2);
if nsd==3
  wze=we(ne3);
end
if isTimeDependent
  woldxe=wolde(ne1,:);
  woldye=wolde(ne2,:);
  if nsd==3
    woldze=wolde(ne3,:);
  end
end
if nsd==2
  wxyze=[wxe,wye];
elseif nsd==3
  wxyze=[wxe,wye,wze];
end
pe=p(re,wxyze,ee);
dpdue=dpdu(re,wxyze,ee);
dpdre=dpdue(:,1);
dpdwxe=dpdue(:,2);
dpdwye=dpdue(:,3);
if nsd==3
  dpdwze=dpdue(:,4);
end
if nsd==2
  dpde=dpdue(:,4);
elseif nsd==3
  dpde=dpdue(:,5);
end
if nsd==2
  sxxe=-Voigt1*Lxxe-Voigt2*Lyye;
  syye=-Voigt1*Lyye-Voigt2*Lxxe;
  sxye=-Voigt3*Lxye;
elseif nsd==3
  sxxe=-Voigt1*Lxxe-Voigt2*Lyye-Voigt2*Lzze;
  syye=-Voigt1*Lyye-Voigt2*Lxxe-Voigt2*Lzze;
  szze=-Voigt1*Lzze-Voigt2*Lxxe-Voigt2*Lyye;
  sxye=-Voigt3*Lxye;
  sxze=-Voigt3*Lxze;
  syze=-Voigt3*Lyze;
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
Lxxeg=Ne*Lxxe;
Lyyeg=Ne*Lyye;
Lxyeg=Ne*Lxye;
if nsd==3
  Lzzeg=Ne*Lzze;
  Lxzeg=Ne*Lxze;
  Lyzeg=Ne*Lyze;
end
Qxeg=Ne*Qxe;
Qyeg=Ne*Qye;
if nsd==3
  Qzeg=Ne*Qze;
end
reg=Ne*re;
wxeg=Ne*wxe;
wyeg=Ne*wye;
if nsd==3
  wzeg=Ne*wze;
end
eeg=Ne*ee;
peg=Ne*pe;
dpdueg=Ne*dpdue;
dpdreg=dpdueg(:,1);
dpdwxeg=dpdueg(:,2);
dpdwyeg=dpdueg(:,3);
if nsd==3
  dpdwzeg=dpdueg(:,4);
end
if nsd==2
  dpdeeg=dpdueg(:,4);
elseif nsd==3
  dpdeeg=dpdueg(:,5);
end
dpdxeg=Nex*pe;
dpdyeg=Ney*pe;
if nsd==3
  dpdzeg=Nez*pe;
end
ddpdrdxeg=Nex*dpdre;
ddpdrdyeg=Ney*dpdre;
if nsd==3
  ddpdrdzeg=Nez*dpdre;
end
ddpdwxdxeg=Nex*dpdwxe;
ddpdwxdyeg=Ney*dpdwxe;
ddpdwydxeg=Nex*dpdwye;
ddpdwydyeg=Ney*dpdwye;
if nsd==3
  ddpdwxdzeg=Nez*dpdwxe;
  ddpdwydzeg=Nez*dpdwye;
  ddpdwzdxeg=Nex*dpdwze;
  ddpdwzdyeg=Ney*dpdwze;
  ddpdwzdzeg=Nez*dpdwze;
end
ddpdedxeg=Nex*dpde;
ddpdedyeg=Ney*dpde;
if nsd==3
  ddpdedzeg=Nez*dpde;
end
if nsd==2
  wxyzeg=[wxeg,wyeg];
elseif nsd==3
  wxyzeg=[wxeg,wyeg,wzeg];
end
Teg=T(reg,wxyzeg,eeg);
dTdueg=dTdu(reg,wxyzeg,eeg);
dTdreg=dTdueg(:,1);
dTdwxeg=dTdueg(:,2);
dTdwyeg=dTdueg(:,3);
if nsd==3
  dTdwzeg=dTdueg(:,4);
end
if nsd==2
  dTdeeg=dTdueg(:,4);
elseif nsd==3
  dTdeeg=dTdueg(:,5);
end
dvxdxeg=Nex*(wxe./re);
dvxdyeg=Ney*(wxe./re);
dvydxeg=Nex*(wye./re);
dvydyeg=Ney*(wye./re);
if nsd==3
  dvxdzeg=Nez*(wxe./re);
  dvydzeg=Nez*(wye./re);
  dvzdxeg=Nex*(wze./re);
  dvzdyeg=Ney*(wze./re);
  dvzdzeg=Nez*(wze./re);
end
vxeg=wxeg./reg;
vyeg=wyeg./reg;
if nsd==3
  vzeg=wzeg./reg;
end
dwxrm2dxeg=Nex*(wxe./re.^2);
dwyrm2dyeg=Ney*(wye./re.^2);
if nsd==3
  dwzrm2dzeg=Nez*(wze./re.^2);
end
drm1dxeg=Nex*(1./re);
drm1dyeg=Ney*(1./re);
if nsd==3
  drm1dzeg=Nez*(1./re);
end
sxxeg=Ne*sxxe;
syyeg=Ne*syye;
sxyeg=Ne*sxye;
if nsd==3
  szzeg=Ne*szze;
  sxzeg=Ne*sxze;
  syzeg=Ne*syze;
end
if nsd==2
  drm2wsdxeg=Nex*(1./re.^2.*(wxe.*sxxe+wye.*sxye));
  drm2wsdyeg=Ney*(1./re.^2.*(wxe.*sxye+wye.*syye));
elseif nsd==3
  drm2wsdxeg=Nex*(1./re.^2.*(wxe.*sxxe+wye.*sxye+wze.*sxze));
  drm2wsdyeg=Ney*(1./re.^2.*(wxe.*sxye+wye.*syye+wze.*syze));
  drm2wsdzeg=Nez*(1./re.^2.*(wxe.*sxze+wye.*syze+wze.*szze));
end
if nsd==2
  rm2wsxeg=1./reg.^2.*(wxeg.*sxxeg+wyeg.*sxyeg);
  rm2wsyeg=1./reg.^2.*(wxeg.*sxyeg+wyeg.*syyeg);
elseif nsd==3
  rm2wsxeg=1./reg.^2.*(wxeg.*sxxeg+wyeg.*sxyeg+wzeg.*sxzeg);
  rm2wsyeg=1./reg.^2.*(wxeg.*sxyeg+wyeg.*syyeg+wzeg.*syzeg);
  rm2wszeg=1./reg.^2.*(wxeg.*sxzeg+wyeg.*syzeg+wzeg.*szzeg);
end
drm1sxxdxeg=Nex*(1./re.*sxxe);
drm1sxydxeg=Nex*(1./re.*sxye);
drm1sxydyeg=Ney*(1./re.*sxye);
drm1syydyeg=Ney*(1./re.*syye);
if nsd==3
  drm1sxzdxeg=Nex*(1./re.*sxze);
  drm1sxzdzeg=Nez*(1./re.*sxze);
  drm1syzdyeg=Ney*(1./re.*syze);
  drm1syzdzeg=Nez*(1./re.*syze);
  drm1szzdzeg=Nez*(1./re.*szze);
end
rm1sxxeg=1./reg.*sxxeg;
rm1sxyeg=1./reg.*sxyeg;
rm1syyeg=1./reg.*syyeg;
if nsd==3
  rm1sxzeg=1./reg.*sxzeg;
  rm1syzeg=1./reg.*syzeg;
  rm1szzeg=1./reg.*szzeg;
end
if isTimeDependent
  roldeg=Ne*rolde;
  woldxeg=Ne*woldxe;
  woldyeg=Ne*woldye;
  if nsd==3
    woldzeg=Ne*woldze;
  end
  eoldeg=Ne*eolde;
end
Rceg=Rc(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
feg=f(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
fxeg=feg(:,1);
fyeg=feg(:,2);
if nsd==3
  fzeg=feg(:,3);
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
  KLr(ne1,ne1)=-Voigt1*NwexT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt2*NweyT*((wyeg.*1./reg.^2).*Ne);
  KLr(ne2,ne1)=-Voigt2*NwexT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt1*NweyT*((wyeg.*1./reg.^2).*Ne);
  KLr(ne3,ne1)=-Voigt3*NweyT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt3*NwexT*((wyeg.*1./reg.^2).*Ne);
elseif nsd==3
  KLr(ne1,ne1)=-Voigt1*NwexT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt2*NweyT*((wyeg.*1./reg.^2).*Ne)...
               -Voigt2*NwezT*((wzeg.*1./reg.^2).*Ne);
  KLr(ne2,ne1)=-Voigt2*NwexT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt1*NweyT*((wyeg.*1./reg.^2).*Ne)...
               -Voigt2*NwezT*((wzeg.*1./reg.^2).*Ne);
  KLr(ne3,ne1)=-Voigt2*NwexT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt2*NweyT*((wyeg.*1./reg.^2).*Ne)...
               -Voigt1*NwezT*((wzeg.*1./reg.^2).*Ne);
  KLr(ne4,ne1)=-Voigt3*NweyT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt3*NwexT*((wyeg.*1./reg.^2).*Ne);
  KLr(ne5,ne1)=-Voigt3*NwezT*((wxeg.*1./reg.^2).*Ne)...
               -Voigt3*NwexT*((wzeg.*1./reg.^2).*Ne);
  KLr(ne6,ne1)=-Voigt3*NwezT*((wyeg.*1./reg.^2).*Ne)...
               -Voigt3*NweyT*((wzeg.*1./reg.^2).*Ne);
end

if nsd==2
  KLw(ne1,ne1)=Voigt1*NwexT*((1./reg).*Ne);
  KLw(ne2,ne1)=Voigt2*NwexT*((1./reg).*Ne);
  KLw(ne3,ne1)=Voigt3*NweyT*((1./reg).*Ne);
  KLw(ne1,ne2)=Voigt2*NweyT*((1./reg).*Ne);
  KLw(ne2,ne2)=Voigt1*NweyT*((1./reg).*Ne);
  KLw(ne3,ne2)=Voigt3*NwexT*((1./reg).*Ne);
elseif nsd==3
  KLw(ne1,ne1)=Voigt1*NwexT*((1./reg).*Ne);
  KLw(ne2,ne1)=Voigt2*NwexT*((1./reg).*Ne);
  KLw(ne3,ne1)=Voigt2*NwexT*((1./reg).*Ne);
  KLw(ne4,ne1)=Voigt3*NweyT*((1./reg).*Ne);
  KLw(ne5,ne1)=Voigt3*NwezT*((1./reg).*Ne);
  KLw(ne1,ne2)=Voigt2*NweyT*((1./reg).*Ne);
  KLw(ne2,ne2)=Voigt1*NweyT*((1./reg).*Ne);
  KLw(ne3,ne2)=Voigt2*NweyT*((1./reg).*Ne);
  KLw(ne4,ne2)=Voigt3*NwexT*((1./reg).*Ne);
  KLw(ne6,ne2)=Voigt3*NwezT*((1./reg).*Ne);
  KLw(ne1,ne3)=Voigt2*NwezT*((1./reg).*Ne);
  KLw(ne2,ne3)=Voigt2*NwezT*((1./reg).*Ne);
  KLw(ne3,ne3)=Voigt1*NwezT*((1./reg).*Ne);
  KLw(ne5,ne3)=Voigt3*NwexT*((1./reg).*Ne);
  KLw(ne6,ne3)=Voigt3*NweyT*((1./reg).*Ne);
end

KQQ(ne1,ne1)=-Me;
KQQ(ne2,ne2)=-Me;
if nsd==3
  KQQ(ne3,ne3)=-Me;
end

KQr(ne1,ne1)=sqrt(kappa)*NwexT*((dTdreg).*Ne);
KQr(ne2,ne1)=sqrt(kappa)*NweyT*((dTdreg).*Ne);
if nsd==3
  KQr(ne3,ne1)=sqrt(kappa)*NwezT*((dTdreg).*Ne);
end

if nsd==2
  KQw(ne1,ne1)=sqrt(kappa)*NwexT*((dTdwxeg).*Ne);
  KQw(ne2,ne1)=sqrt(kappa)*NweyT*((dTdwxeg).*Ne);
  KQw(ne1,ne2)=sqrt(kappa)*NwexT*((dTdwyeg).*Ne);
  KQw(ne2,ne2)=sqrt(kappa)*NweyT*((dTdwyeg).*Ne);
elseif nsd==3
  KQw(ne1,ne1)=sqrt(kappa)*NwexT*((dTdwxeg).*Ne);
  KQw(ne2,ne1)=sqrt(kappa)*NweyT*((dTdwxeg).*Ne);
  KQw(ne3,ne1)=sqrt(kappa)*NwezT*((dTdwxeg).*Ne);
  KQw(ne1,ne2)=sqrt(kappa)*NwexT*((dTdwyeg).*Ne);
  KQw(ne2,ne2)=sqrt(kappa)*NweyT*((dTdwyeg).*Ne);
  KQw(ne3,ne2)=sqrt(kappa)*NwezT*((dTdwyeg).*Ne);
  KQw(ne1,ne3)=sqrt(kappa)*NwexT*((dTdwzeg).*Ne);
  KQw(ne2,ne3)=sqrt(kappa)*NweyT*((dTdwzeg).*Ne);
  KQw(ne3,ne3)=sqrt(kappa)*NwezT*((dTdwzeg).*Ne);
end

KQe(ne1,ne1)=sqrt(kappa)*NwexT*((dTdeeg).*Ne);
KQe(ne2,ne1)=sqrt(kappa)*NweyT*((dTdeeg).*Ne);
if nsd==3
  KQe(ne3,ne1)=sqrt(kappa)*NwezT*((dTdeeg).*Ne);
end

if isTimeDependent
  Krr(ne1,ne1)=alpha(1)/dt*Me;
end

Krw(ne1,ne1)=-Cxe;
Krw(ne1,ne2)=-Cye;
if nsd==3
  Krw(ne1,ne3)=-Cze;
end

if isTimeDependent
  Kww(ne1,ne1)=alpha(1)/dt*Me;
  Kww(ne2,ne2)=alpha(1)/dt*Me;
  if nsd==3
    Kww(ne3,ne3)=alpha(1)/dt*Me;
  end
end

if isConvectiveFlow
  if nsd==2
    Kwr(ne1,ne1)=+NwexT*((wxeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wxeg.*wyeg./reg.^2).*Ne);
    Kwr(ne2,ne1)=+NwexT*((wyeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wyeg.*wyeg./reg.^2).*Ne);
  elseif nsd==3
    Kwr(ne1,ne1)=+NwexT*((wxeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wxeg.*wyeg./reg.^2).*Ne)...
                 +NwezT*((wxeg.*wzeg./reg.^2).*Ne);
    Kwr(ne2,ne1)=+NwexT*((wyeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wyeg.*wyeg./reg.^2).*Ne)...
                 +NwezT*((wyeg.*wzeg./reg.^2).*Ne);
    Kwr(ne3,ne1)=+NwexT*((wzeg.*wxeg./reg.^2).*Ne)...
                 +NweyT*((wzeg.*wyeg./reg.^2).*Ne)...
                 +NwezT*((wzeg.*wzeg./reg.^2).*Ne);
  end
  
  if nsd==2
    Kww(ne1,ne1)=Kww(ne1,ne1)-NwexT*((2*wxeg./reg).*Ne)...
                             -NweyT*((  wyeg./reg).*Ne);
    Kww(ne2,ne1)=Kww(ne2,ne1)-NwexT*((  wyeg./reg).*Ne);
    Kww(ne1,ne2)=Kww(ne1,ne2)-NweyT*((  wxeg./reg).*Ne);
    Kww(ne2,ne2)=Kww(ne2,ne2)-NwexT*((  wxeg./reg).*Ne)...
                             -NweyT*((2*wyeg./reg).*Ne);
  elseif nsd==3
    Kww(ne1,ne1)=Kww(ne1,ne1)-NwexT*((2*wxeg./reg).*Ne)...
                             -NweyT*((  wyeg./reg).*Ne)...
                             -NwezT*((  wzeg./reg).*Ne);
    Kww(ne2,ne1)=Kww(ne2,ne1)-NwexT*((  wyeg./reg).*Ne);
    Kww(ne3,ne1)=Kww(ne3,ne1)-NwezT*((  wzeg./reg).*Ne);
    Kww(ne1,ne2)=Kww(ne1,ne2)-NweyT*((  wxeg./reg).*Ne);
    Kww(ne2,ne2)=Kww(ne2,ne2)-NwexT*((  wxeg./reg).*Ne)...
                             -NweyT*((2*wyeg./reg).*Ne)...
                             -NwezT*((  wzeg./reg).*Ne);
    Kww(ne3,ne2)=Kww(ne3,ne2)-NwezT*((  wzeg./reg).*Ne);
    Kww(ne1,ne3)=Kww(ne1,ne3)-NwezT*((  wxeg./reg).*Ne);
    Kww(ne2,ne3)=Kww(ne2,ne3)-NwezT*((  wyeg./reg).*Ne);
    Kww(ne3,ne3)=Kww(ne3,ne3)-NwexT*((  wxeg./reg).*Ne)...
                             -NweyT*((  wyeg./reg).*Ne)...
                             -NwezT*((2*wzeg./reg).*Ne);
  end
end

if nsd==2
  KwL(ne1,ne1)=Voigt1*CxeT;
  KwL(ne2,ne1)=Voigt2*CyeT;
  KwL(ne1,ne2)=Voigt2*CxeT;
  KwL(ne2,ne2)=Voigt1*CyeT;
  KwL(ne1,ne3)=Voigt3*CyeT;
  KwL(ne2,ne3)=Voigt3*CxeT;
elseif nsd==3
  KwL(ne1,ne1)=Voigt1*CxeT;
  KwL(ne2,ne1)=Voigt2*CyeT;
  KwL(ne3,ne1)=Voigt2*CzeT;
  KwL(ne1,ne2)=Voigt2*CxeT;
  KwL(ne2,ne2)=Voigt1*CyeT;
  KwL(ne3,ne2)=Voigt2*CzeT;
  KwL(ne1,ne3)=Voigt2*CxeT;
  KwL(ne2,ne3)=Voigt2*CyeT;
  KwL(ne3,ne3)=Voigt1*CzeT;
  KwL(ne1,ne4)=Voigt3*CyeT;
  KwL(ne2,ne4)=Voigt3*CxeT;
  KwL(ne1,ne5)=Voigt3*CzeT;
  KwL(ne3,ne5)=Voigt3*CxeT;
  KwL(ne2,ne6)=Voigt3*CzeT;
  KwL(ne3,ne6)=Voigt3*CyeT;
end

Kwr(ne1,ne1)=Kwr(ne1,ne1)+NweT*((ddpdrdxeg).*Ne+(dpdreg).*Nex);
Kwr(ne2,ne1)=Kwr(ne2,ne1)+NweT*((ddpdrdyeg).*Ne+(dpdreg).*Ney);
if nsd==3
  Kwr(ne3,ne1)=Kwr(ne3,ne1)+NweT*((ddpdrdzeg).*Ne+(dpdreg).*Nez);
end

if nsd==2
  Kww(ne1,ne1)=Kww(ne1,ne1)+NweT*((ddpdwxdxeg).*Ne+(dpdwxeg).*Nex);
  Kww(ne2,ne1)=Kww(ne2,ne1)+NweT*((ddpdwxdyeg).*Ne+(dpdwxeg).*Ney);
  Kww(ne1,ne2)=Kww(ne1,ne2)+NweT*((ddpdwydxeg).*Ne+(dpdwyeg).*Nex);
  Kww(ne2,ne2)=Kww(ne2,ne2)+NweT*((ddpdwydyeg).*Ne+(dpdwyeg).*Ney);
elseif nsd==3
  Kww(ne1,ne1)=Kww(ne1,ne1)+NweT*((ddpdwxdxeg).*Ne+(dpdwxeg).*Nex);
  Kww(ne2,ne1)=Kww(ne2,ne1)+NweT*((ddpdwxdyeg).*Ne+(dpdwxeg).*Ney);
  Kww(ne3,ne1)=Kww(ne3,ne1)+NweT*((ddpdwxdzeg).*Ne+(dpdwxeg).*Nez);
  Kww(ne1,ne2)=Kww(ne1,ne2)+NweT*((ddpdwydxeg).*Ne+(dpdwyeg).*Nex);
  Kww(ne2,ne2)=Kww(ne2,ne2)+NweT*((ddpdwydyeg).*Ne+(dpdwyeg).*Ney);
  Kww(ne3,ne2)=Kww(ne3,ne2)+NweT*((ddpdwydzeg).*Ne+(dpdwyeg).*Nez);
  Kww(ne1,ne3)=Kww(ne1,ne3)+NweT*((ddpdwzdxeg).*Ne+(dpdwzeg).*Nex);
  Kww(ne2,ne3)=Kww(ne2,ne3)+NweT*((ddpdwzdyeg).*Ne+(dpdwzeg).*Ney);
  Kww(ne3,ne3)=Kww(ne3,ne3)+NweT*((ddpdwzdzeg).*Ne+(dpdwzeg).*Nez);
end

Kwe(ne1,ne1)=NweT*((ddpdedxeg).*Ne+(dpdeeg).*Nex);
Kwe(ne2,ne1)=NweT*((ddpdedyeg).*Ne+(dpdeeg).*Ney);
if nsd==3
  Kwe(ne3,ne1)=NweT*((ddpdedzeg).*Ne+(dpdeeg).*Nez);
end

if isTimeDependent
  Kee(ne1,ne1)=alpha(1)/dt*Me;
end

Ker(ne1,ne1)=+NwexT*((eeg.*wxeg./reg.^2).*Ne)...
             +NweyT*((eeg.*wyeg./reg.^2).*Ne);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+NwezT*((eeg.*wzeg./reg.^2).*Ne);
end

Kew(ne1,ne1)=-NwexT*((eeg./reg).*Ne);
Kew(ne1,ne2)=-NweyT*((eeg./reg).*Ne);
if nsd==3
  Kew(ne1,ne3)=-NwezT*((eeg./reg).*Ne);
end

Kee(ne1,ne1)=Kee(ne1,ne1)-NwexT*((wxeg./reg).*Ne)...
                         -NweyT*((wyeg./reg).*Ne);
if nsd==3
  Kee(ne1,ne1)=Kee(ne1,ne1)-NwezT*((wzeg./reg).*Ne);
end

if nsd==2
  KeL(ne1,ne1)=+Voigt1*NweT*((dvxdxeg).*Ne+(vxeg).*Nex)...
               +Voigt2*NweT*((dvydyeg).*Ne+(vyeg).*Ney);
  KeL(ne1,ne2)=+Voigt2*NweT*((dvxdxeg).*Ne+(vxeg).*Nex)...
               +Voigt1*NweT*((dvydyeg).*Ne+(vyeg).*Ney);
  KeL(ne1,ne3)=+Voigt3*NweT*((dvxdyeg).*Ne+(vyeg).*Nex)...
               +Voigt3*NweT*((dvydxeg).*Ne+(vxeg).*Ney);
elseif nsd==3
  KeL(ne1,ne1)=+Voigt1*NweT*((dvxdxeg).*Ne+(vxeg).*Nex)...
               +Voigt2*NweT*((dvydyeg).*Ne+(vyeg).*Ney)...
               +Voigt2*NweT*((dvzdzeg).*Ne+(vzeg).*Nez);
  KeL(ne1,ne2)=+Voigt2*NweT*((dvxdxeg).*Ne+(vxeg).*Nex)...
               +Voigt1*NweT*((dvydyeg).*Ne+(vyeg).*Ney)...
               +Voigt2*NweT*((dvzdzeg).*Ne+(vzeg).*Nez);
  KeL(ne1,ne3)=+Voigt2*NweT*((dvxdxeg).*Ne+(vxeg).*Nex)...
               +Voigt2*NweT*((dvydyeg).*Ne+(vyeg).*Ney)...
               +Voigt1*NweT*((dvzdzeg).*Ne+(vzeg).*Nez);
  KeL(ne1,ne4)=+Voigt3*NweT*((dvxdyeg).*Ne+(vyeg).*Nex)...
               +Voigt3*NweT*((dvydxeg).*Ne+(vxeg).*Ney);
  KeL(ne1,ne5)=+Voigt3*NweT*((dvxdzeg).*Ne+(vzeg).*Nex)...
               +Voigt3*NweT*((dvzdxeg).*Ne+(vxeg).*Nez);
  KeL(ne1,ne6)=+Voigt3*NweT*((dvydzeg).*Ne+(vzeg).*Ney)...
               +Voigt3*NweT*((dvzdyeg).*Ne+(vyeg).*Nez);
end

Ker(ne1,ne1)=Ker(ne1,ne1)+NweT*(+(drm2wsdxeg).*Ne+(rm2wsxeg).*Nex...
                                +(drm2wsdyeg).*Ne+(rm2wsyeg).*Ney);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+NweT*((drm2wsdzeg).*Ne+(rm2wszeg).*Nez);
end

if nsd==2
  Kew(ne1,ne1)=Kew(ne1,ne1)-NweT*(+(drm1sxxdxeg+drm1sxydyeg).*Ne...
                                  +(rm1sxxeg).*Nex+(rm1sxyeg).*Ney);
  Kew(ne1,ne2)=Kew(ne1,ne2)-NweT*(+(drm1sxydxeg+drm1syydyeg).*Ne...
                                  +(rm1sxyeg).*Nex+(rm1syyeg).*Ney);
elseif nsd==3
  Kew(ne1,ne1)=Kew(ne1,ne1)-NweT*(+(drm1sxxdxeg+drm1sxydyeg+drm1sxzdzeg).*Ne...
                                  +(rm1sxxeg).*Nex+(rm1sxyeg).*Ney+(rm1sxzeg).*Nez);
  Kew(ne1,ne2)=Kew(ne1,ne2)-NweT*(+(drm1sxydxeg+drm1syydyeg+drm1syzdzeg).*Ne...
                                  +(rm1sxyeg).*Nex+(rm1syyeg).*Ney+(rm1syzeg).*Nez);
  Kew(ne1,ne3)=Kew(ne1,ne3)-NweT*(+(drm1sxzdxeg+drm1syzdyeg+drm1szzdzeg).*Ne...
                                  +(rm1sxzeg).*Nex+(rm1syzeg).*Ney+(rm1szzeg).*Nez);
end

Ker(ne1,ne1)=Ker(ne1,ne1)+NweT*(+(+ddpdrdxeg.*wxeg./reg+ddpdrdyeg.*wyeg./reg...
                                  +dpdreg.*dvxdxeg+dpdreg.*dvydyeg...
                                  -dpdxeg.*wxeg./reg.^2-dpdyeg.*wyeg./reg.^2 ...
                                  -peg.*dwxrm2dxeg-peg.*dwyrm2dyeg).*Ne...
                                +(dpdreg.*wxeg./reg+peg.*wxeg./reg.^2).*Nex...
                                +(dpdreg.*wyeg./reg+peg.*wyeg./reg.^2).*Ney);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+NweT*(+(+ddpdrdzeg.*wzeg./reg...
                                    +dpdreg.*dvzdzeg...
                                    -dpdzeg.*wzeg./reg.^2 ...
                                    -peg.*dwzrm2dzeg).*Ne...
                                  +(dpdreg.*wzeg./reg+peg.*wzeg./reg.^2).*Nez);
end

if nsd==2
  Kew(ne1,ne1)=Kew(ne1,ne1)+NweT*(...
                           +(+ddpdwxdxeg.*wxeg./reg+ddpdwxdyeg.*wyeg./reg...
                             +dpdwxeg.*dvxdxeg+dpdwxeg.*dvydyeg...
                             +dpdxeg./reg...
                             +peg.*drm1dxeg).*Ne...
                           +(dpdwxeg.*wxeg./reg+peg./reg).*Nex...
                           +(dpdwxeg.*wyeg./reg).*Ney);
  Kew(ne1,ne2)=Kew(ne1,ne2)+NweT*(...
                           +(+ddpdwydxeg.*wxeg./reg+ddpdwydyeg.*wyeg./reg...
                             +dpdwyeg.*dvxdxeg+dpdwyeg.*dvydyeg...
                             +dpdyeg./reg...
                             +peg.*drm1dyeg).*Ne...
                           +(dpdwyeg.*wxeg./reg).*Nex...
                           +(dpdwyeg.*wyeg./reg+peg./reg).*Ney);
elseif nsd==3
  Kew(ne1,ne1)=Kew(ne1,ne1)+NweT*(...
                           +(+ddpdwxdxeg.*wxeg./reg+ddpdwxdyeg.*wyeg./reg+ddpdwxdzeg.*wzeg./reg...
                             +dpdwxeg.*dvxdxeg+dpdwxeg.*dvydyeg+dpdwxeg.*dvzdzeg...
                             +dpdxeg./reg...
                             +peg.*drm1dxeg).*Ne...
                           +(dpdwxeg.*wxeg./reg+peg./reg).*Nex...
                           +(dpdwxeg.*wyeg./reg).*Ney...
                           +(dpdwxeg.*wzeg./reg).*Nez);
  Kew(ne1,ne2)=Kew(ne1,ne2)+NweT*(...
                           +(+ddpdwydxeg.*wxeg./reg+ddpdwydyeg.*wyeg./reg+ddpdwydzeg.*wzeg./reg...
                             +dpdwyeg.*dvxdxeg+dpdwyeg.*dvydyeg+dpdwyeg.*dvzdzeg...
                             +dpdyeg./reg...
                             +peg.*drm1dyeg).*Ne...
                           +(dpdwyeg.*wxeg./reg).*Nex...
                           +(dpdwyeg.*wyeg./reg+peg./reg).*Ney...
                           +(dpdwyeg.*wzeg./reg).*Nez);
  Kew(ne1,ne3)=Kew(ne1,ne3)+NweT*(...
                           +(+ddpdwzdxeg.*wxeg./reg+ddpdwzdyeg.*wyeg./reg+ddpdwzdzeg.*wzeg./reg...
                             +dpdwzeg.*dvxdxeg+dpdwzeg.*dvydyeg+dpdwzeg.*dvzdzeg...
                             +dpdzeg./reg...
                             +peg.*drm1dzeg).*Ne...
                           +(dpdwzeg.*wxeg./reg).*Nex...
                           +(dpdwzeg.*wyeg./reg).*Ney...
                           +(dpdwzeg.*wzeg./reg+peg./reg).*Nez);
end

Kee(ne1,ne1)=Kee(ne1,ne1)+NweT*(+(+ddpdedxeg.*wxeg./reg+ddpdedyeg.*wyeg./reg...
                                  +dpdeeg.*dvxdxeg+dpdeeg.*dvydyeg).*Ne...
                                +(dpdeeg.*wxeg./reg).*Nex...
                                +(dpdeeg.*wyeg./reg).*Ney);
if nsd==3
  Kee(ne1,ne1)=Kee(ne1,ne1)+NweT*(+(+ddpdedzeg.*wzeg./reg...
                                    +dpdeeg.*dvzdzeg).*Ne...
                                  +(dpdeeg.*wzeg./reg).*Nez);
end

KeQ(ne1,ne1)=sqrt(kappa)*CxeT;
KeQ(ne1,ne2)=sqrt(kappa)*CyeT;
if nsd==3
  KeQ(ne1,ne3)=sqrt(kappa)*CzeT;
end

if nsd==2
  Ker(ne1,ne1)=Ker(ne1,ne1)+NweT*(((fxeg.*wxeg+fyeg.*wyeg)./reg.^2).*Ne);
elseif nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+NweT*(((fxeg.*wxeg+fyeg.*wyeg+fzeg.*wzeg)./reg.^2).*Ne);
end

Kew(ne1,ne1)=Kew(ne1,ne1)-NweT*((fxeg./reg).*Ne);
Kew(ne1,ne2)=Kew(ne1,ne2)-NweT*((fyeg./reg).*Ne);
if nsd==3
  Kew(ne1,ne3)=Kew(ne1,ne3)-NweT*((fzeg./reg).*Ne);
end

% Compute rhs
if nsd==2
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(Voigt1*wxeg./reg)...
            -NweyT*(Voigt2*wyeg./reg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(Voigt2*wxeg./reg)...
            -NweyT*(Voigt1*wyeg./reg);
  fL(ne3,1)=+NweT*(Lxyeg)...
            -NweyT*(Voigt3*wxeg./reg)...
            -NwexT*(Voigt3*wyeg./reg);
elseif nsd==3
  fL(ne1,1)=+NweT*(Lxxeg)...
            -NwexT*(Voigt1*wxeg./reg)...
            -NweyT*(Voigt2*wyeg./reg)...
            -NwezT*(Voigt2*wzeg./reg);
  fL(ne2,1)=+NweT*(Lyyeg)...
            -NwexT*(Voigt2*wxeg./reg)...
            -NweyT*(Voigt1*wyeg./reg)...
            -NwezT*(Voigt2*wzeg./reg);
  fL(ne3,1)=+NweT*(Lzzeg)...
            -NwexT*(Voigt2*wxeg./reg)...
            -NweyT*(Voigt2*wyeg./reg)...
            -NwezT*(Voigt1*wzeg./reg);
  fL(ne4,1)=+NweT*(Lxyeg)...
            -NwexT*(Voigt3*wyeg./reg)...
            -NweyT*(Voigt3*wxeg./reg);
  fL(ne5,1)=+NweT*(Lxzeg)...
            -NwexT*(Voigt3*wzeg./reg)...
            -NwezT*(Voigt3*wxeg./reg);
  fL(ne6,1)=+NweT*(Lyzeg)...
            -NweyT*(Voigt3*wzeg./reg)...
            -NwezT*(Voigt3*wyeg./reg);
end

fQ(ne1,1)=+NweT*(Qxeg)...
          -NwexT*(sqrt(kappa)*Teg);
fQ(ne2,1)=+NweT*(Qyeg)...
          -NweyT*(sqrt(kappa)*Teg);
if nsd==3
  fQ(ne3,1)=+NweT*(Qzeg)...
            -NwezT*(sqrt(kappa)*Teg);
end

if isTimeDependent
  fr(ne1,1)=-NweT*(1/dt*reg*alpha(1)...
                  +1/dt*roldeg*alpha(2:BDFo+1,1));
end

fr(ne1,1)=fr(ne1,1)+NwexT*(wxeg)...
                   +NweyT*(wyeg);
if nsd==3
  fr(ne1,1)=fr(ne1,1)+NwezT*(wzeg);
end

fr(ne1,1)=fr(ne1,1)+NweT*(Rceg);

if isTimeDependent
  fw(ne1,1)=-NweT*(1/dt*wxeg*alpha(1)...
                  +1/dt*woldxeg*alpha(2:BDFo+1,1));
  fw(ne2,1)=-NweT*(1/dt*wyeg*alpha(1)...
                  +1/dt*woldyeg*alpha(2:BDFo+1,1));
  if nsd==3
    fw(ne3,1)=-NweT*(1/dt*wzeg*alpha(1)...
                    +1/dt*woldzeg*alpha(2:BDFo+1,1));
  end
end

if isConvectiveFlow
  if nsd==2
    fw(ne1,1)=fw(ne1,1)+NwexT*(wxeg.*wxeg./reg)...
                       +NweyT*(wxeg.*wyeg./reg);
    fw(ne2,1)=fw(ne2,1)+NwexT*(wyeg.*wxeg./reg)...
                       +NweyT*(wyeg.*wyeg./reg);
  elseif nsd==3
    fw(ne1,1)=fw(ne1,1)+NwexT*(wxeg.*wxeg./reg)...
                       +NweyT*(wxeg.*wyeg./reg)...
                       +NwezT*(wxeg.*wzeg./reg);
    fw(ne2,1)=fw(ne2,1)+NwexT*(wyeg.*wxeg./reg)...
                       +NweyT*(wyeg.*wyeg./reg)...
                       +NwezT*(wyeg.*wzeg./reg);
    fw(ne3,1)=fw(ne3,1)+NwexT*(wzeg.*wxeg./reg)...
                       +NweyT*(wzeg.*wyeg./reg)...
                       +NwezT*(wzeg.*wzeg./reg);
  end
end

if nsd==2
  fw(ne1,1)=fw(ne1,1)+NweT*((-Voigt1*(Nex*Lxxe)-Voigt2*(Nex*Lyye)...
                             -Voigt3*(Ney*Lxye)...
                             -Nex*pe...
                             +fxeg));
  fw(ne2,1)=fw(ne2,1)+NweT*((-Voigt2*(Ney*Lxxe)-Voigt1*(Ney*Lyye)...
                             -Voigt3*(Nex*Lxye)...
                             -Ney*pe...
                             +fyeg));
elseif nsd==3
  fw(ne1,1)=fw(ne1,1)+NweT*((-Voigt1*(Nex*Lxxe)-Voigt2*(Nex*Lyye)-Voigt2*(Nex*Lzze)...
                             -Voigt3*(Ney*Lxye)-Voigt3*(Nez*Lxze)...
                             -Nex*pe...
                             +fxeg));
  fw(ne2,1)=fw(ne2,1)+NweT*((-Voigt2*(Ney*Lxxe)-Voigt1*(Ney*Lyye)-Voigt2*(Ney*Lzze)...
                             -Voigt3*(Nex*Lxye)-Voigt3*(Nez*Lyze)...
                             -Ney*pe...
                             +fyeg));
  fw(ne3,1)=fw(ne3,1)+NweT*((-Voigt2*(Nez*Lxxe)-Voigt2*(Nez*Lyye)-Voigt1*(Nez*Lzze)...
                             -Voigt3*(Nex*Lxze)-Voigt3*(Ney*Lyze)...
                             -Nez*pe...
                             +fzeg));
end

if isTimeDependent
  fe(ne1,1)=-NweT*(1/dt*eeg*alpha(1)...
                  +1/dt*eoldeg*alpha(2:BDFo+1,1));
end

fe(ne1,1)=fe(ne1,1)+NwexT*(eeg.*wxeg./reg)...
                   +NweyT*(eeg.*wyeg./reg);
if nsd==3
  fe(ne1,1)=fe(ne1,1)+NwezT*(eeg.*wzeg./reg);
end

if nsd==2
  fe(ne1,1)=fe(ne1,1)+NweT*((-Nex*(1./re.*(+wxe.*(Voigt1*Lxxe+Voigt2*Lyye)...
                                           +wye.*(Voigt3*Lxye)))...
                             -Ney*(1./re.*(+wxe.*(Voigt3*Lxye)...
                                           +wye.*(Voigt2*Lxxe+Voigt1*Lyye)))...
                             -Nex*(pe.*wxe./re)...
                             -Ney*(pe.*wye./re)...
                             -Nex*(sqrt(kappa)*Qxe)...
                             -Ney*(sqrt(kappa)*Qye)...
                             +(fxeg.*wxeg+fyeg.*wyeg)./reg...
                             +seg));
elseif nsd==3
  fe(ne1,1)=fe(ne1,1)+NweT*((-Nex*(1./re.*(+wxe.*(Voigt1*Lxxe+Voigt2*Lyye+Voigt2*Lzze)...
                                           +wye.*(Voigt3*Lxye)...
                                           +wze.*(Voigt3*Lxze)))...
                             -Ney*(1./re.*(+wxe.*(Voigt3*Lxye)...
                                           +wye.*(Voigt2*Lxxe+Voigt1*Lyye+Voigt2*Lzze)...
                                           +wze.*(Voigt3*Lyze)))...
                             -Nez*(1./re.*(+wxe.*(Voigt3*Lxze)...
                                           +wye.*(Voigt3*Lyze)...
                                           +wze.*(Voigt2*Lxxe+Voigt2*Lyye+Voigt1*Lzze)))...
                             -Nex*(pe.*wxe./re)...
                             -Ney*(pe.*wye./re)...
                             -Nez*(pe.*wze./re)...
                             -Nex*(sqrt(kappa)*Qxe)...
                             -Ney*(sqrt(kappa)*Qye)...
                             -Nez*(sqrt(kappa)*Qze)...
                             +(fxeg.*wxeg+fyeg.*wyeg+fzeg.*wzeg)./reg...
                             +seg));
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf,nsd);
    
    % Check boundary
    isDirichlet=Faces.Dirichlet(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nf4=nf3+NumElementNodes;
    nf5=nf4+NumElementNodes;
    nf6=nf5+NumElementNodes;
    nefU1=(iFace-1)*(1+nsd+1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
    nefU4=nefU3+NumFaceNodes;
    nefU5=nefU4+NumFaceNodes;
    nefR1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefW1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefW2=nefW1+NumFaceNodes;
    nefW3=nefW2+NumFaceNodes;
    nefE1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    
    % Flip face
    Node2Match1stNode1=Faces.Interior(2,iFace);
    FlipFace=max(Node2Match1stNode1);
    if FlipFace
      order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
      nefU1=nefU1(order);
      nefU2=nefU2(order);
      nefU3=nefU3(order);
      nefU4=nefU4(order);
      nefU5=nefU5(order);
      nefR1=nefR1(order);
      nefW1=nefW1(order);
      nefW2=nefW2(order);
      nefW3=nefW3(order);
      nefE1=nefE1(order);
    end
    
    % Compute variables at nodes
    Lxxf=Lxxe(nf1);
    Lyyf=Lyye(nf1);
    Lxyf=Lxye(nf1);
    if nsd==3
      Lzzf=Lzze(nf1);
      Lxzf=Lxze(nf1);
      Lyzf=Lyze(nf1);
    end
    Qxf=Qxe(nf1);
    Qyf=Qye(nf1);
    if nsd==3
      Qzf=Qze(nf1);
    end
    rf=re(nf1);
    wxf=wxe(nf1);
    wyf=wye(nf1);
    if nsd==3
      wzf=wze(nf1);
    end
    ef=ee(nf1);
    Rf=Ue(nefU1);
    Wxf=Ue(nefU2);
    Wyf=Ue(nefU3);
    if nsd==3
      Wzf=Ue(nefU4);
    end
    if nsd==2
      Ef=Ue(nefU4);
    elseif nsd==3
      Ef=Ue(nefU5);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    Lxxfg=Nf*Lxxf;
    Lyyfg=Nf*Lyyf;
    Lxyfg=Nf*Lxyf;
    if nsd==3
      Lzzfg=Nf*Lzzf;
      Lxzfg=Nf*Lxzf;
      Lyzfg=Nf*Lyzf;
    end
    Qxfg=Nf*Qxf;
    Qyfg=Nf*Qyf;
    if nsd==3
      Qzfg=Nf*Qzf;
    end
    rfg=Nf*rf;
    wxfg=Nf*wxf;
    wyfg=Nf*wyf;
    if nsd==3
      wzfg=Nf*wzf;
    end
    efg=Nf*ef;
    Rfg=Nf*Rf;
    Wxfg=Nf*Wxf;
    Wyfg=Nf*Wyf;
    if nsd==3
      Wzfg=Nf*Wzf;
    end
    Efg=Nf*Ef;
    if nsd==2
      Wxyzfg=[Wxfg,Wyfg];
    elseif nsd==3
      Wxyzfg=[Wxfg,Wyfg,Wzfg];
    end
    pfg=p(Rfg,Wxyzfg,Efg);
    dpdUfg=dpdu(Rfg,Wxyzfg,Efg);
    dpdRfg=dpdUfg(:,1);
    dpdWxfg=dpdUfg(:,2);
    dpdWyfg=dpdUfg(:,3);
    if nsd==3
      dpdWzfg=dpdUfg(:,4);
    end
    if nsd==2
      dpdEfg=dpdUfg(:,4);
    elseif nsd==3
      dpdEfg=dpdUfg(:,5);
    end
    Tfg=T(Rfg,Wxyzfg,Efg);
    dTdUfg=dTdu(Rfg,Wxyzfg,Efg);
    dTdRfg=dTdUfg(:,1);
    dTdWxfg=dTdUfg(:,2);
    dTdWyfg=dTdUfg(:,3);
    if nsd==3
      dTdWzfg=dTdUfg(:,4);
    end
    if nsd==2
      dTdEfg=dTdUfg(:,4);
    elseif nsd==3
      dTdEfg=dTdUfg(:,5);
    end
    if isDirichlet
      rDfg=rD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      wDfg=wD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      wDxfg=wDfg(:,1);
      wDyfg=wDfg(:,2);
      if nsd==3
        wDzfg=wDfg(:,3);
      end
      eDfg=eD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
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
    
    % Compute lhs
    if nsd==2
      KLR(nf1,nefR1)=KLR(nf1,nefR1)+NwfT*((+Voigt1*(Wxfg./Rfg.^2.*nx)...
                                           +Voigt2*(Wyfg./Rfg.^2.*ny)).*Nf);
      KLR(nf2,nefR1)=KLR(nf2,nefR1)+NwfT*((+Voigt2*(Wxfg./Rfg.^2.*nx)...
                                           +Voigt1*(Wyfg./Rfg.^2.*ny)).*Nf);
      KLR(nf3,nefR1)=KLR(nf3,nefR1)+NwfT*((+Voigt3*(Wxfg./Rfg.^2.*ny)...
                                           +Voigt3*(Wyfg./Rfg.^2.*nx)).*Nf);
    elseif nsd==3
      KLR(nf1,nefR1)=KLR(nf1,nefR1)+NwfT*((+Voigt1*(Wxfg./Rfg.^2.*nx)...
                                           +Voigt2*(Wyfg./Rfg.^2.*ny)...
                                           +Voigt2*(Wzfg./Rfg.^2.*nz)).*Nf);
      KLR(nf2,nefR1)=KLR(nf2,nefR1)+NwfT*((+Voigt2*(Wxfg./Rfg.^2.*nx)...
                                           +Voigt1*(Wyfg./Rfg.^2.*ny)...
                                           +Voigt2*(Wzfg./Rfg.^2.*nz)).*Nf);
      KLR(nf3,nefR1)=KLR(nf3,nefR1)+NwfT*((+Voigt2*(Wxfg./Rfg.^2.*nx)...
                                           +Voigt2*(Wyfg./Rfg.^2.*ny)...
                                           +Voigt1*(Wzfg./Rfg.^2.*nz)).*Nf);
      KLR(nf4,nefR1)=KLR(nf4,nefR1)+NwfT*((+Voigt3*(Wxfg./Rfg.^2.*ny)...
                                           +Voigt3*(Wyfg./Rfg.^2.*nx)).*Nf);
      KLR(nf5,nefR1)=KLR(nf5,nefR1)+NwfT*((+Voigt3*(Wxfg./Rfg.^2.*nz)...
                                           +Voigt3*(Wzfg./Rfg.^2.*nx)).*Nf);
      KLR(nf6,nefR1)=KLR(nf6,nefR1)+NwfT*((+Voigt3*(Wyfg./Rfg.^2.*nz)...
                                           +Voigt3*(Wzfg./Rfg.^2.*ny)).*Nf);
    end
      
    if nsd==2
      KLW(nf1,nefW1)=KLW(nf1,nefW1)-Voigt1*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf2,nefW1)=KLW(nf2,nefW1)-Voigt2*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf3,nefW1)=KLW(nf3,nefW1)-Voigt3*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf1,nefW2)=KLW(nf1,nefW2)-Voigt2*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf2,nefW2)=KLW(nf2,nefW2)-Voigt1*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf3,nefW2)=KLW(nf3,nefW2)-Voigt3*NwfT*((1./Rfg.*nx).*Nf);
    elseif nsd==3
      KLW(nf1,nefW1)=KLW(nf1,nefW1)-Voigt1*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf2,nefW1)=KLW(nf2,nefW1)-Voigt2*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf3,nefW1)=KLW(nf3,nefW1)-Voigt2*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf4,nefW1)=KLW(nf4,nefW1)-Voigt3*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf5,nefW1)=KLW(nf5,nefW1)-Voigt3*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf1,nefW2)=KLW(nf1,nefW2)-Voigt2*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf2,nefW2)=KLW(nf2,nefW2)-Voigt1*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf3,nefW2)=KLW(nf3,nefW2)-Voigt2*NwfT*((1./Rfg.*ny).*Nf);
      KLW(nf4,nefW2)=KLW(nf4,nefW2)-Voigt3*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf6,nefW2)=KLW(nf6,nefW2)-Voigt3*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf1,nefW3)=KLW(nf1,nefW3)-Voigt2*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf2,nefW3)=KLW(nf2,nefW3)-Voigt2*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf3,nefW3)=KLW(nf3,nefW3)-Voigt1*NwfT*((1./Rfg.*nz).*Nf);
      KLW(nf5,nefW3)=KLW(nf5,nefW3)-Voigt3*NwfT*((1./Rfg.*nx).*Nf);
      KLW(nf6,nefW3)=KLW(nf6,nefW3)-Voigt3*NwfT*((1./Rfg.*ny).*Nf);
    end
    
    KQR(nf1,nefR1)=KQR(nf1,nefR1)-sqrt(kappa)*NwfT*((dTdRfg.*nx).*Nf);
    KQR(nf2,nefR1)=KQR(nf2,nefR1)-sqrt(kappa)*NwfT*((dTdRfg.*ny).*Nf);
    if nsd==3
      KQR(nf3,nefR1)=KQR(nf3,nefR1)-sqrt(kappa)*NwfT*((dTdRfg.*nz).*Nf);
    end
    
    KQW(nf1,nefW1)=KQW(nf1,nefW1)-sqrt(kappa)*NwfT*((dTdWxfg.*nx).*Nf);
    KQW(nf2,nefW1)=KQW(nf2,nefW1)-sqrt(kappa)*NwfT*((dTdWxfg.*ny).*Nf);
    KQW(nf1,nefW2)=KQW(nf1,nefW2)-sqrt(kappa)*NwfT*((dTdWyfg.*nx).*Nf);
    KQW(nf2,nefW2)=KQW(nf2,nefW2)-sqrt(kappa)*NwfT*((dTdWyfg.*ny).*Nf);
    if nsd==3
      KQW(nf3,nefW1)=KQW(nf3,nefW1)-sqrt(kappa)*NwfT*((dTdWxfg.*nz).*Nf);
      KQW(nf3,nefW2)=KQW(nf3,nefW2)-sqrt(kappa)*NwfT*((dTdWyfg.*nz).*Nf);
      KQW(nf1,nefW3)=KQW(nf1,nefW3)-sqrt(kappa)*NwfT*((dTdWzfg.*nx).*Nf);
      KQW(nf2,nefW3)=KQW(nf2,nefW3)-sqrt(kappa)*NwfT*((dTdWzfg.*ny).*Nf);
      KQW(nf3,nefW3)=KQW(nf3,nefW3)-sqrt(kappa)*NwfT*((dTdWzfg.*nz).*Nf);
    end
    
    KQE(nf1,nefE1)=KQE(nf1,nefE1)-sqrt(kappa)*NwfT*((dTdEfg.*nx).*Nf);
    KQE(nf2,nefE1)=KQE(nf2,nefE1)-sqrt(kappa)*NwfT*((dTdEfg.*ny).*Nf);
    if nsd==3
      KQE(nf3,nefE1)=KQE(nf3,nefE1)-sqrt(kappa)*NwfT*((dTdEfg.*nz).*Nf);
    end
    
    Krr(nf1,nf1)=Krr(nf1,nf1)+tauR*Mf;
    
    KrW(nf1,nefW1)=KrW(nf1,nefW1)+Mfnx;
    KrW(nf1,nefW2)=KrW(nf1,nefW2)+Mfny;
    if nsd==3
      KrW(nf1,nefW3)=KrW(nf1,nefW3)+Mfnz;
    end
    
    KrR(nf1,nefR1)=KrR(nf1,nefR1)-tauR*Mf;
    
    Kww(nf1,nf1)=Kww(nf1,nf1)+tauW*Mf;
    Kww(nf2,nf2)=Kww(nf2,nf2)+tauW*Mf;
    if nsd==3
      Kww(nf3,nf3)=Kww(nf3,nf3)+tauW*Mf;
    end
    
    if isConvectiveFlow
      if nsd==2
        KwR(nf1,nefR1)=KwR(nf1,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Wxfg).*Nf);
        KwR(nf2,nefR1)=KwR(nf2,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Wyfg).*Nf);
      elseif nsd==3
        KwR(nf1,nefR1)=KwR(nf1,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wxfg).*Nf);
        KwR(nf2,nefR1)=KwR(nf2,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wyfg).*Nf);
        KwR(nf3,nefR1)=KwR(nf3,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wzfg).*Nf);
      end

      if nsd==2
        KwW(nf1,nefW1)=KwW(nf1,nefW1)+NwfT*((2*Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Nf);
        KwW(nf2,nefW1)=KwW(nf2,nefW1)+NwfT*((Wyfg./Rfg.*nx).*Nf);
        KwW(nf1,nefW2)=KwW(nf1,nefW2)+NwfT*((Wxfg./Rfg.*ny).*Nf);
        KwW(nf2,nefW2)=KwW(nf2,nefW2)+NwfT*((Wxfg./Rfg.*nx+2*Wyfg./Rfg.*ny).*Nf);
      elseif nsd==3
        KwW(nf1,nefW1)=KwW(nf1,nefW1)+NwfT*((2*Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Nf);
        KwW(nf2,nefW1)=KwW(nf2,nefW1)+NwfT*((Wyfg./Rfg.*nx).*Nf);
        KwW(nf3,nefW1)=KwW(nf3,nefW1)+NwfT*((Wzfg./Rfg.*nx).*Nf);
        KwW(nf1,nefW2)=KwW(nf1,nefW2)+NwfT*((Wxfg./Rfg.*ny).*Nf);
        KwW(nf2,nefW2)=KwW(nf2,nefW2)+NwfT*((Wxfg./Rfg.*nx+2*Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Nf);
        KwW(nf3,nefW2)=KwW(nf3,nefW2)+NwfT*((Wzfg./Rfg.*ny).*Nf);
        KwW(nf1,nefW3)=KwW(nf1,nefW3)+NwfT*((Wxfg./Rfg.*nz).*Nf);
        KwW(nf2,nefW3)=KwW(nf2,nefW3)+NwfT*((Wyfg./Rfg.*nz).*Nf);
        KwW(nf3,nefW3)=KwW(nf3,nefW3)+NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+2*Wzfg./Rfg.*nz).*Nf);
      end
    end
    
    KwW(nf1,nefW1)=KwW(nf1,nefW1)-tauW*Mf;
    KwW(nf2,nefW2)=KwW(nf2,nefW2)-tauW*Mf;
    if nsd==3
      KwW(nf3,nefW3)=KwW(nf3,nefW3)-tauW*Mf;
    end
    
    Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*Mf;
    
    if nsd==2
      KeR(nf1,nefR1)=KeR(nf1,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Efg).*Nf);
    elseif nsd==3
      KeR(nf1,nefR1)=KeR(nf1,nefR1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Efg).*Nf);
    end
    
    KeW(nf1,nefW1)=KeW(nf1,nefW1)+NwfT*((Efg./Rfg.*nx).*Nf);
    KeW(nf1,nefW2)=KeW(nf1,nefW2)+NwfT*((Efg./Rfg.*ny).*Nf);
    if nsd==3
      KeW(nf1,nefW3)=KeW(nf1,nefW3)+NwfT*((Efg./Rfg.*nz).*Nf);
    end
    
    if nsd==2
      KeE(nf1,nefE1)=KeE(nf1,nefE1)+NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg).*Nf);
    elseif nsd==3
      KeE(nf1,nefE1)=KeE(nf1,nefE1)+NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg).*Nf);
    end
    
    KeE(nf1,nefE1)=KeE(nf1,nefE1)-tauE*Mf;
    
    if not(isDirichlet)
      KRW(nefR1,nefW1)=KRW(nefR1,nefW1)-Mfnx;
      KRW(nefR1,nefW2)=KRW(nefR1,nefW2)-Mfny;
      if nsd==3
        KRW(nefR1,nefW3)=KRW(nefR1,nefW3)-Mfnz;
      end
    end
    
    if not(isDirichlet)
      KRr(nefR1,nf1)=KRr(nefR1,nf1)-tauR*Mf;
    end
    
    KRR(nefR1,nefR1)=KRR(nefR1,nefR1)+tauR*Mf;
    
    if not(isDirichlet) && isConvectiveFlow
      if nsd==2
        KWR(nefW1,nefR1)=KWR(nefW1,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Wxfg).*Nf);
        KWR(nefW2,nefR1)=KWR(nefW2,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Wyfg).*Nf);
      elseif nsd==3
        KWR(nefW1,nefR1)=KWR(nefW1,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wxfg).*Nf);
        KWR(nefW2,nefR1)=KWR(nefW2,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wyfg).*Nf);
        KWR(nefW3,nefR1)=KWR(nefW3,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Wzfg).*Nf);
      end
    end
    
    if not(isDirichlet) && isConvectiveFlow
      if nsd==2
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-NwfT*((2*Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-NwfT*((Wyfg./Rfg.*nx).*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-NwfT*((Wxfg./Rfg.*ny).*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-NwfT*((Wxfg./Rfg.*nx+2*Wyfg./Rfg.*ny).*Nf);
      elseif nsd==3
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-NwfT*((2*Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-NwfT*((Wyfg./Rfg.*nx).*Nf);
        KWW(nefW3,nefW1)=KWW(nefW3,nefW1)-NwfT*((Wzfg./Rfg.*nx).*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-NwfT*((Wxfg./Rfg.*ny).*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-NwfT*((Wxfg./Rfg.*nx+2*Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Nf);
        KWW(nefW3,nefW2)=KWW(nefW3,nefW2)-NwfT*((Wzfg./Rfg.*ny).*Nf);
        KWW(nefW1,nefW3)=KWW(nefW1,nefW3)-NwfT*((Wxfg./Rfg.*nz).*Nf);
        KWW(nefW2,nefW3)=KWW(nefW2,nefW3)-NwfT*((Wyfg./Rfg.*nz).*Nf);
        KWW(nefW3,nefW3)=KWW(nefW3,nefW3)-NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+2*Wzfg./Rfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KWL(nefW1,nf1)=KWL(nefW1,nf1)-Voigt1*Mfnx;
        KWL(nefW2,nf1)=KWL(nefW2,nf1)-Voigt2*Mfny;
        KWL(nefW1,nf2)=KWL(nefW1,nf2)-Voigt2*Mfnx;
        KWL(nefW2,nf2)=KWL(nefW2,nf2)-Voigt1*Mfny;
        KWL(nefW1,nf3)=KWL(nefW1,nf3)-Voigt3*Mfny;
        KWL(nefW2,nf3)=KWL(nefW2,nf3)-Voigt3*Mfnx;
      elseif nsd==3
        KWL(nefW1,nf1)=KWL(nefW1,nf1)-Voigt1*Mfnx;
        KWL(nefW2,nf1)=KWL(nefW2,nf1)-Voigt2*Mfny;
        KWL(nefW3,nf1)=KWL(nefW3,nf1)-Voigt2*Mfnz;
        KWL(nefW1,nf2)=KWL(nefW1,nf2)-Voigt2*Mfnx;
        KWL(nefW2,nf2)=KWL(nefW2,nf2)-Voigt1*Mfny;
        KWL(nefW3,nf2)=KWL(nefW3,nf2)-Voigt2*Mfnz;
        KWL(nefW1,nf3)=KWL(nefW1,nf3)-Voigt2*Mfnx;
        KWL(nefW2,nf3)=KWL(nefW2,nf3)-Voigt2*Mfny;
        KWL(nefW3,nf3)=KWL(nefW3,nf3)-Voigt1*Mfnz;
        KWL(nefW1,nf4)=KWL(nefW1,nf4)-Voigt3*Mfny;
        KWL(nefW2,nf4)=KWL(nefW2,nf4)-Voigt3*Mfnx;
        KWL(nefW1,nf5)=KWL(nefW1,nf5)-Voigt3*Mfnz;
        KWL(nefW3,nf5)=KWL(nefW3,nf5)-Voigt3*Mfnx;
        KWL(nefW2,nf6)=KWL(nefW2,nf6)-Voigt3*Mfnz;
        KWL(nefW3,nf6)=KWL(nefW3,nf6)-Voigt3*Mfny;
      end
    end
    
    if not(isDirichlet)
      KWR(nefW1,nefR1)=KWR(nefW1,nefR1)-NwfT*((dpdRfg.*nx).*Nf);
      KWR(nefW2,nefR1)=KWR(nefW2,nefR1)-NwfT*((dpdRfg.*ny).*Nf);
      if nsd==3
        KWR(nefW3,nefR1)=KWR(nefW3,nefR1)-NwfT*((dpdRfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-NwfT*((dpdWxfg.*nx).*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-NwfT*((dpdWxfg.*ny).*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-NwfT*((dpdWyfg.*nx).*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-NwfT*((dpdWyfg.*ny).*Nf);
      elseif nsd==3
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-NwfT*((dpdWxfg.*nx).*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-NwfT*((dpdWxfg.*ny).*Nf);
        KWW(nefW3,nefW1)=KWW(nefW3,nefW1)-NwfT*((dpdWxfg.*nz).*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-NwfT*((dpdWyfg.*nx).*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-NwfT*((dpdWyfg.*ny).*Nf);
        KWW(nefW3,nefW2)=KWW(nefW3,nefW2)-NwfT*((dpdWyfg.*nz).*Nf);
        KWW(nefW1,nefW3)=KWW(nefW1,nefW3)-NwfT*((dpdWzfg.*nx).*Nf);
        KWW(nefW2,nefW3)=KWW(nefW2,nefW3)-NwfT*((dpdWzfg.*ny).*Nf);
        KWW(nefW3,nefW3)=KWW(nefW3,nefW3)-NwfT*((dpdWzfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet)
      KWE(nefW1,nefE1)=KWE(nefW1,nefE1)-NwfT*((dpdEfg.*nx).*Nf);
      KWE(nefW2,nefE1)=KWE(nefW2,nefE1)-NwfT*((dpdEfg.*ny).*Nf);
      if nsd==3
        KWE(nefW3,nefE1)=KWE(nefW3,nefE1)-NwfT*((dpdEfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet)
      KWw(nefW1,nf1)=KWw(nefW1,nf1)-tauW*Mf;
      KWw(nefW2,nf2)=KWw(nefW2,nf2)-tauW*Mf;
      if nsd==3
        KWw(nefW3,nf3)=KWw(nefW3,nf3)-tauW*Mf;
      end
    end
    
    KWW(nefW1,nefW1)=KWW(nefW1,nefW1)+tauW*Mf;
    KWW(nefW2,nefW2)=KWW(nefW2,nefW2)+tauW*Mf;
    if nsd==3
      KWW(nefW3,nefW3)=KWW(nefW3,nefW3)+tauW*Mf;
    end
    
    if not(isDirichlet)
      if nsd==2
        KER(nefE1,nefR1)=KER(nefE1,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg.^2.*Efg).*Nf);
      elseif nsd==3
        KER(nefE1,nefR1)=KER(nefE1,nefR1)+NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg.^2.*Efg).*Nf);
      end
    end
    
    if not(isDirichlet)
      KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-NwfT*((Efg./Rfg.*nx).*Nf);
      KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-NwfT*((Efg./Rfg.*ny).*Nf);
      if nsd==3
        KEW(nefE1,nefW3)=KEW(nefE1,nefW3)-NwfT*((Efg./Rfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KEE(nefE1,nefE1)=KEE(nefE1,nefE1)-NwfT*(((Wxfg.*nx+Wyfg.*ny)./Rfg).*Nf);
      elseif nsd==3
        KEE(nefE1,nefE1)=KEE(nefE1,nefE1)-NwfT*(((Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)./Rfg).*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KEL(nefE1,nf1)=KEL(nefE1,nf1)-NwfT*(((Voigt1*Wxfg.*nx...
                                             +Voigt2*Wyfg.*ny)./Rfg).*Nf);
        KEL(nefE1,nf2)=KEL(nefE1,nf2)-NwfT*(((Voigt2*Wxfg.*nx...
                                             +Voigt1*Wyfg.*ny)./Rfg).*Nf);
        KEL(nefE1,nf3)=KEL(nefE1,nf3)-NwfT*(((Voigt3*Wyfg.*nx...
                                             +Voigt3*Wxfg.*ny)./Rfg).*Nf);
      elseif nsd==3
        KEL(nefE1,nf1)=KEL(nefE1,nf1)-NwfT*(((Voigt1*Wxfg.*nx...
                                             +Voigt2*Wyfg.*ny...
                                             +Voigt2*Wzfg*nz)./Rfg).*Nf);
        KEL(nefE1,nf2)=KEL(nefE1,nf2)-NwfT*(((Voigt2*Wxfg.*nx...
                                             +Voigt1*Wyfg.*ny...
                                             +Voigt2*Wzfg*nz)./Rfg).*Nf);
        KEL(nefE1,nf3)=KEL(nefE1,nf3)-NwfT*(((Voigt2*Wxfg.*nx...
                                             +Voigt2*Wyfg.*ny...
                                             +Voigt1*Wzfg*nz)./Rfg).*Nf);
        KEL(nefE1,nf4)=KEL(nefE1,nf4)-NwfT*(((Voigt3*Wyfg.*nx...
                                             +Voigt3*Wxfg.*ny)./Rfg).*Nf);
        KEL(nefE1,nf5)=KEL(nefE1,nf5)-NwfT*(((Voigt3*Wzfg.*nx...
                                             +Voigt3*Wxfg.*nz)./Rfg).*Nf);
        KEL(nefE1,nf6)=KEL(nefE1,nf6)-NwfT*(((Voigt3*Wzfg.*ny...
                                             +Voigt3*Wyfg.*nz)./Rfg).*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KER(nefE1,nefR1)=KER(nefE1,nefR1)+NwfT*((...
          +(+(Wxfg.*(Voigt1*Lxxfg+Voigt2*Lyyfg)+Voigt3*Wyfg.*Lxyfg).*nx...
            +(Wyfg.*(Voigt2*Lxxfg+Voigt1*Lyyfg)+Voigt3*Wxfg.*Lxyfg).*ny)./Rfg.^2 ...
          -dpdRfg./Rfg.*(Wxfg.*nx+Wyfg.*ny)...
          +pfg./Rfg.^2.*(Wxfg.*nx+Wyfg.*ny)).*Nf);
      elseif nsd==3
        KER(nefE1,nefR1)=KER(nefE1,nefR1)+NwfT*((...
          +(+(Wxfg.*(Voigt1*Lxxfg+Voigt2*Lyyfg+Voigt2*Lzzfg)...
                    +Voigt3*Wyfg.*Lxyfg+Voigt3*Wzfg.*Lxzfg).*nx...
            +(Wyfg.*(Voigt2*Lxxfg+Voigt1*Lyyfg+Voigt2*Lzzfg)...
                    +Voigt3*Wxfg.*Lxyfg+Voigt3*Wzfg.*Lyzfg).*ny...
            +(Wzfg.*(Voigt2*Lxxfg+Voigt2*Lyyfg+Voigt1*Lzzfg)...
                    +Voigt3*Wxfg.*Lxzfg+Voigt3*Wzfg.*Lyzfg).*nz)./Rfg.^2 ...
          -dpdRfg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)...
          +pfg./Rfg.^2.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)).*Nf);
      end
    end

    if not(isDirichlet)
      if nsd==2
        KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-NwfT*((...
          +((Voigt1*Lxxfg+Voigt2*Lyyfg).*nx+Voigt3*Lxyfg.*ny)./Rfg...
          +dpdWxfg./Rfg.*(Wxfg.*nx+Wyfg.*ny)...
          +pfg./Rfg.*nx).*Nf);
        KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-NwfT*((...
          +((Voigt2*Lxxfg+Voigt1*Lyyfg).*ny+Voigt3*Lxyfg.*nx)./Rfg...
          +dpdWyfg./Rfg.*(Wxfg.*nx+Wyfg.*ny)...
          +pfg./Rfg.*ny).*Nf);
      elseif nsd==3
        KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-NwfT*((...
          +((Voigt1*Lxxfg+Voigt2*Lyyfg+Voigt2*Lzzfg).*nx+Voigt3*Lxyfg.*ny+Voigt3*Lxzfg.*nz)./Rfg...
          +dpdWxfg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)...
          +pfg./Rfg.*nx).*Nf);
        KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-NwfT*((...
          +((Voigt2*Lxxfg+Voigt1*Lyyfg+Voigt2*Lzzfg).*ny+Voigt3*Lxyfg.*nx+Voigt3*Lyzfg.*nz)./Rfg...
          +dpdWyfg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)...
          +pfg./Rfg.*ny).*Nf);
        KEW(nefE1,nefW3)=KEW(nefE1,nefW3)-NwfT*((...
          +((Voigt2*Lxxfg+Voigt2*Lyyfg+Voigt1*Lzzfg).*nz+Voigt3*Lxzfg.*nx+Voigt3*Lyzfg.*ny)./Rfg...
          +dpdWzfg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)...
          +pfg./Rfg.*nz).*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KEE(nefE1,nefE1)=KEE(nefE1,nefE1)-NwfT*((dpdEfg./Rfg.*(Wxfg.*nx+Wyfg.*ny)).*Nf);
      elseif nsd==3
        KEE(nefE1,nefE1)=KEE(nefE1,nefE1)-NwfT*((dpdEfg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)).*Nf);
      end
    end
    
    if not(isDirichlet)
      KEQ(nefE1,nf1)=KEQ(nefE1,nf1)-sqrt(kappa)*Mfnx;
      KEQ(nefE1,nf2)=KEQ(nefE1,nf2)-sqrt(kappa)*Mfny;
      if nsd==3
        KEQ(nefE1,nf3)=KEQ(nefE1,nf3)-sqrt(kappa)*Mfnz;
      end
    end
    
    if not(isDirichlet)
      KEe(nefE1,nf1)=KEe(nefE1,nf1)-tauE*Mf;
    end
    
    KEE(nefE1,nefE1)=KEE(nefE1,nefE1)+tauE*Mf;
    
    % Compute rhs
    if nsd==2
      fL(nf1,1)=fL(nf1,1)+NwfT*((+Voigt1*nx.*Wxfg./Rfg...
                                 +Voigt2*ny.*Wyfg./Rfg));
      fL(nf2,1)=fL(nf2,1)+NwfT*((+Voigt2*nx.*Wxfg./Rfg...
                                 +Voigt1*ny.*Wyfg./Rfg));
      fL(nf3,1)=fL(nf3,1)+NwfT*((+Voigt3*ny.*Wxfg./Rfg...
                                 +Voigt3*nx.*Wyfg./Rfg));
    elseif nsd==3
      fL(nf1,1)=fL(nf1,1)+NwfT*((+Voigt1*nx.*Wxfg./Rfg...
                                 +Voigt2*ny.*Wyfg./Rfg...
                                 +Voigt2*nz.*Wzfg./Rfg));
      fL(nf2,1)=fL(nf2,1)+NwfT*((+Voigt2*nx.*Wxfg./Rfg...
                                 +Voigt1*ny.*Wyfg./Rfg...
                                 +Voigt2*nz.*Wzfg./Rfg));
      fL(nf3,1)=fL(nf3,1)+NwfT*((+Voigt2*nx.*Wxfg./Rfg...
                                 +Voigt2*ny.*Wyfg./Rfg...
                                 +Voigt1*nz.*Wzfg./Rfg));
      fL(nf4,1)=fL(nf4,1)+NwfT*((+Voigt3*nx.*Wyfg./Rfg...
                                 +Voigt3*ny.*Wxfg./Rfg));
      fL(nf5,1)=fL(nf5,1)+NwfT*((+Voigt3*nx.*Wzfg./Rfg...
                                 +Voigt3*nz.*Wxfg./Rfg));
      fL(nf6,1)=fL(nf6,1)+NwfT*((+Voigt3*ny.*Wzfg./Rfg...
                                 +Voigt3*nz.*Wyfg./Rfg));
    end
    
    fQ(nf1,1)=fQ(nf1,1)+NwfT*(sqrt(kappa)*nx.*Tfg);
    fQ(nf2,1)=fQ(nf2,1)+NwfT*(sqrt(kappa)*ny.*Tfg);
    if nsd==3
      fQ(nf3,1)=fQ(nf3,1)+NwfT*(sqrt(kappa)*nz.*Tfg);
    end
    
    fr(nf1,1)=fr(nf1,1)-NwfT*((+tauR*rfg...
                               +Wxfg.*nx+Wyfg.*ny...
                               -tauR*Rfg));
    if nsd==3
      fr(nf1,1)=fr(nf1,1)-NwfT*(Wzfg.*nz);
    end
    
    fw(nf1,1)=fw(nf1,1)-NwfT*(tauW*(wxfg-Wxfg));
    fw(nf2,1)=fw(nf2,1)-NwfT*(tauW*(wyfg-Wyfg));
    if nsd==3
      fw(nf3,1)=fw(nf3,1)-NwfT*(tauW*(wzfg-Wzfg));
    end
    
    if isConvectiveFlow
      if nsd==2
        fw(nf1,1)=fw(nf1,1)-NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wxfg);
        fw(nf2,1)=fw(nf2,1)-NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wyfg);
      elseif nsd==3
        fw(nf1,1)=fw(nf1,1)-NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wxfg);
        fw(nf2,1)=fw(nf2,1)-NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wyfg);
        fw(nf3,1)=fw(nf3,1)-NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wzfg);
      end
    end
    
    fe(nf1,1)=fe(nf1,1)-NwfT*((+tauE*(efg-Efg)...
                               +Efg./Rfg.*(Wxfg.*nx+Wyfg.*ny)));
    if nsd==3
      fe(nf1,1)=fe(nf1,1)-NwfT*(Efg./Rfg.*(Wzfg.*nz));
    end
    
    if not(isDirichlet)
      fR(nefR1,1)=fR(nefR1,1)+NwfT*((+Wxfg.*nx+Wyfg.*ny...
                                     +tauR*rfg));
      if nsd==3
        fR(nefR1,1)=fR(nefR1,1)+NwfT*(Wzfg.*nz);
      end
    end
    
    fR(nefR1,1)=fR(nefR1,1)-NwfT*(tauR*Rfg);
    
    if isDirichlet
      fR(nefR1,1)=fR(nefR1,1)+NwfT*(tauR*rDfg);
    end
    
    if not(isDirichlet) && isConvectiveFlow
      if nsd==2
        fW(nefW1,1)=fW(nefW1,1)+NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wxfg);
        fW(nefW2,1)=fW(nefW2,1)+NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wyfg);
      elseif nsd==3
        fW(nefW1,1)=fW(nefW1,1)+NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wxfg);
        fW(nefW2,1)=fW(nefW2,1)+NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wyfg);
        fW(nefW3,1)=fW(nefW3,1)+NwfT*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wzfg);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        fW(nefW1,1)=fW(nefW1,1)+NwfT*((+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg...
                                       +Voigt3*ny.*Lxyfg...
                                       +pfg.*nx...
                                       +tauW*wxfg));
        fW(nefW2,1)=fW(nefW2,1)+NwfT*((+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg...
                                       +Voigt3*nx.*Lxyfg...
                                       +pfg.*ny...
                                       +tauW*wyfg));
      elseif nsd==3
        fW(nefW1,1)=fW(nefW1,1)+NwfT*((+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg+Voigt2*nx.*Lzzfg...
                                       +Voigt3*ny.*Lxyfg+Voigt3*nz.*Lxzfg...
                                       +pfg.*nx...
                                       +tauW*wxfg));
        fW(nefW2,1)=fW(nefW2,1)+NwfT*((+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg+Voigt2*ny.*Lzzfg...
                                       +Voigt3*nx.*Lxyfg+Voigt3*nz.*Lyzfg...
                                       +pfg.*ny...
                                       +tauW*wyfg));
        fW(nefW3,1)=fW(nefW3,1)+NwfT*((+Voigt2*nz.*Lxxfg+Voigt2*nz.*Lyyfg+Voigt1*nz.*Lzzfg...
                                       +Voigt3*nx.*Lxzfg+Voigt3*ny.*Lyzfg...
                                       +pfg.*nz...
                                       +tauW*wzfg));
      end
    end
    
    fW(nefW1,1)=fW(nefW1,1)-NwfT*(tauW*Wxfg);
    fW(nefW2,1)=fW(nefW2,1)-NwfT*(tauW*Wyfg);
    if nsd==3
      fW(nefW3,1)=fW(nefW3,1)-NwfT*(tauW*Wzfg);
    end
    
    if isDirichlet
      fW(nefW1,1)=fW(nefW1,1)+NwfT*(tauW*wDxfg);
      fW(nefW2,1)=fW(nefW2,1)+NwfT*(tauW*wDyfg);
      if nsd==3
        fW(nefW3,1)=fW(nefW3,1)+NwfT*(tauW*wDzfg);
      end
    end
    
    if not(isDirichlet)
      fE(nefE1,1)=fE(nefE1,1)+NwfT*((+Efg./Rfg.*(Wxfg.*nx+Wyfg.*ny)...
                                     +(Wxfg./Rfg.*(Voigt1*Lxxfg+Voigt2*Lyyfg)...
                                      +Wyfg./Rfg.*Voigt3.*Lxyfg).*nx...
                                     +(Wyfg./Rfg.*(Voigt2*Lxxfg+Voigt1*Lyyfg)...
                                      +Wxfg./Rfg.*Voigt3.*Lxyfg).*ny...
                                     +pfg./Rfg.*(Wxfg.*nx+Wyfg.*ny)...
                                     +sqrt(kappa)*(Qxfg.*nx+Qyfg.*ny)...
                                     +tauE*efg));
      if nsd==3
        fE(nefE1,1)=fE(nefE1,1)+NwfT*((+Efg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)...
                                       +(Wxfg./Rfg.*(Voigt1*Lxxfg+Voigt2*Lyyfg+Voigt2*Lzzfg)...
                                        +Wyfg./Rfg.*Voigt3.*Lxyfg+Wzfg./Rfg.*Voigt3.*Lxzfg).*nx...
                                       +(Wyfg./Rfg.*(Voigt2*Lxxfg+Voigt1*Lyyfg+Voigt2*Lzzfg)...
                                        +Wxfg./Rfg.*Voigt3.*Lxyfg+Wzfg./Rfg.*Voigt3.*Lyzfg).*ny...
                                       +(Wzfg./Rfg.*(Voigt2*Lxxfg+Voigt2*Lyyfg+Voigt1*Lzzfg)...
                                        +Wxfg./Rfg.*Voigt3.*Lxzfg+Wyfg./Rfg.*Voigt3.*Lyzfg).*nz...
                                       +pfg./Rfg.*(Wxfg.*nx+Wyfg.*ny+Wzfg.*nz)...
                                       +sqrt(kappa)*(Qxfg.*nx+Qyfg.*ny+Qzfg.*nz)...
                                       +tauE*efg));
      end
    end
    
    fE(nefE1,1)=fE(nefE1,1)-NwfT*(tauE*Efg);
    
    if isDirichlet
      fE(nefE1,1)=fE(nefE1,1)+NwfT*(tauE*eDfg);
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
iL=1:msd*NumElementNodes;
iQ=iL(end)+(1:nsd*NumElementNodes);
ir=iQ(end)+(1:NumElementNodes);
iw=ir(end)+(1:nsd*NumElementNodes);
ie=iw(end)+(1:NumElementNodes);
iR=reshape((0:NumElementFaces-1)*(1+nsd+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iW=reshape((0:NumElementFaces-1)*(1+nsd+1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[])+NumFaceNodes;
iE=reshape((0:NumElementFaces-1)*(1+nsd+1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+(1+nsd)*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((msd+nsd+1+nsd+1)*NumElementNodes,(msd+nsd+1+nsd+1)*NumElementNodes);
LhsLG=zeros((msd+nsd+1+nsd+1)*NumElementNodes,(1+nsd+1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((1+nsd+1)*NumElementFaces*NumFaceNodes,(msd+nsd+1+nsd+1)*NumElementNodes);
LhsGG=zeros((1+nsd+1)*NumElementFaces*NumFaceNodes,(1+nsd+1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((msd+nsd+1+nsd+1)*NumElementNodes,1);
RhsG=zeros((1+nsd+1)*NumElementFaces*NumFaceNodes,1);

% Remove energy equation
if not(SolveEnergyEquation)
  KeL=zeros(size(KeL));
  KeQ=zeros(size(KeQ));
  Ker=zeros(size(Ker));
  Kew=zeros(size(Kew));
  Kee=eye(size(Kee));
  KeR=zeros(size(KeR));
  KeW=zeros(size(KeW));
  KeE=zeros(size(KeE));
  KEL=zeros(size(KEL));
  KEQ=zeros(size(KEQ));
  KEe=zeros(size(KEe));
  KER=zeros(size(KER));
  KEW=zeros(size(KEW));
  KEE=eye(size(KEE));
  fe=zeros(size(fe));
  fE=zeros(size(fE));
end

if matchField(Parameters,'SolveEnergyEquationOnly','yes')
  KLL=eye(size(KLL));
  KLr=zeros(size(KLr));
  KLw=zeros(size(KLw));
  KLR=zeros(size(KLR));
  KLW=zeros(size(KLW));
  KQr=zeros(size(KQr));
  KQw=zeros(size(KQw));
  KQR=zeros(size(KQR));
  KQW=zeros(size(KQW));
  Krr=eye(size(Krr));
  Krw=zeros(size(Krw));
  KrR=zeros(size(KrR));
  KrW=zeros(size(KrW));
  KwL=zeros(size(KwL));
  Kwr=zeros(size(Kwr));
  Kww=eye(size(Kww));
  Kwe=zeros(size(Kwe));
  KwR=zeros(size(KwR));
  KwW=zeros(size(KwW));
  KeL=zeros(size(KeL));
  Ker=zeros(size(Ker));
  Kew=zeros(size(Kew));
  KeR=zeros(size(KeR));
  KeW=zeros(size(KeW));
  KRr=zeros(size(KRr));
  KRR=eye(size(KRR));
  KRW=zeros(size(KRW));
  KWL=zeros(size(KWL));
  KWw=zeros(size(KWw));
  KWR=zeros(size(KWR));
  KWW=eye(size(KWW));
  KWE=zeros(size(KWE));
  KEL=zeros(size(KEL));
  KER=zeros(size(KER));
  KEW=zeros(size(KEW));
  fL=zeros(size(fL));
  fr=zeros(size(fr));
  fw=zeros(size(fw));
  fR=zeros(size(fR));
  fW=zeros(size(fW));
end

% Lhs local-local
LhsLL(iL,iL)=KLL;
LhsLL(iL,ir)=KLr;
LhsLL(iL,iw)=KLw;
LhsLL(iQ,iQ)=KQQ;
LhsLL(iQ,ir)=KQr;
LhsLL(iQ,iw)=KQw;
LhsLL(iQ,ie)=KQe;
LhsLL(ir,ir)=Krr;
LhsLL(ir,iw)=Krw;
LhsLL(iw,iL)=KwL;
LhsLL(iw,ir)=Kwr;
LhsLL(iw,iw)=Kww;
LhsLL(iw,ie)=Kwe;
LhsLL(ie,iL)=KeL;
LhsLL(ie,iQ)=KeQ;
LhsLL(ie,ir)=Ker;
LhsLL(ie,iw)=Kew;
LhsLL(ie,ie)=Kee;

% Lhs local-global
LhsLG(iL,iR)=KLR;
LhsLG(iL,iW)=KLW;
LhsLG(iQ,iR)=KQR;
LhsLG(iQ,iW)=KQW;
LhsLG(iQ,iE)=KQE;
LhsLG(ir,iR)=KrR;
LhsLG(ir,iW)=KrW;
LhsLG(iw,iR)=KwR;
LhsLG(iw,iW)=KwW;
LhsLG(ie,iR)=KeR;
LhsLG(ie,iW)=KeW;
LhsLG(ie,iE)=KeE;

% Rhs local
RhsL(iL,1)=fL;
RhsL(iQ,1)=fQ;
RhsL(ir,1)=fr;
RhsL(iw,1)=fw;
RhsL(ie,1)=fe;

% Lhs global-local
LhsGL(iR,ir)=KRr;
LhsGL(iW,iL)=KWL;
LhsGL(iW,iw)=KWw;
LhsGL(iE,iL)=KEL;
LhsGL(iE,iQ)=KEQ;
LhsGL(iE,ie)=KEe;

% Lhs global-global
LhsGG(iR,iR)=KRR;
LhsGG(iR,iW)=KRW;
LhsGG(iW,iR)=KWR;
LhsGG(iW,iW)=KWW;
LhsGG(iW,iE)=KWE;
LhsGG(iE,iR)=KER;
LhsGG(iE,iW)=KEW;
LhsGG(iE,iE)=KEE;

% Rhs global
RhsG(iR,1)=fR;
RhsG(iW,1)=fW;
RhsG(iE,1)=fE;

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
  Nodes,Faces,SolutionGlobal,SolutionLocal,Parameters,RefElement,Sizes)

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
kappa=Parameters.ThermalConductivity;
T=Parameters.ComputeTemperature;
Xe=Nodes';

% Get solution
Le=reshape(SolutionLocal(:,1:msd),[],1);
Qe=reshape(SolutionLocal(:,msd+(1:nsd)),[],1);
re=reshape(SolutionLocal(:,msd+nsd+1),[],1);
we=reshape(SolutionLocal(:,msd+nsd+1+(1:nsd)),[],1);
ee=reshape(SolutionLocal(:,msd+nsd+1+nsd+1),[],1);
Ue=SolutionGlobal;

% Initialize lhs
KVpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
KVtp=zeros(nsd,nsd*NumElementNodesPost);
KVrp=zeros(qsd,nsd*NumElementNodesPost);
KTpp=zeros(NumElementNodesPost,NumElementNodesPost);
KTtp=zeros(1,NumElementNodesPost);

% Initialize rhs
fVp=zeros(nsd*NumElementNodesPost,1);
fVt=zeros(nsd,1);
fVr=zeros(qsd,1);
fTp=zeros(NumElementNodesPost,1);
fTt=zeros(1,1);

% Get reference data
FaceNodes=RefElement.FaceNodesElem;

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions(1,RefElement.PostLow,RefElement.Post,Xe,nsd);
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
Qxe=Qe(nle1);
Qye=Qe(nle2);
if nsd==3
  Qze=Qe(nle3);
end
wxe=we(nle1);
wye=we(nle2);
if nsd==3
  wze=we(nle3);
end

% Compute variables at Gauss points
Lxxeg=Nle*Lxxe;
Lyyeg=Nle*Lyye;
Lxyeg=Nle*Lxye;
if nsd==3
  Lzzeg=Nle*Lzze;
  Lxzeg=Nle*Lxze;
  Lyzeg=Nle*Lyze;
end
Qxeg=Nle*Qxe;
Qyeg=Nle*Qye;
if nsd==3
  Qzeg=Nle*Qze;
end
reg=Nle*re;
wxeg=Nle*wxe;
wyeg=Nle*wye;
if nsd==3
  wzeg=Nle*wze;
end
eeg=Nle*ee;
if nsd==2
  wxyzeg=[wxeg,wyeg];
elseif nsd==3
  wxyzeg=[wxeg,wyeg,wzeg];
end
Teg=T(reg,wxyzeg,eeg);

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
  KVpp(ne1,ne1)=-Voigt1*Kxxe-Voigt3*Kyye;
  KVpp(ne1,ne2)=-Voigt2*Kxye-Voigt3*Kyxe;
  KVpp(ne2,ne1)=-Voigt2*Kyxe-Voigt3*Kxye;
  KVpp(ne2,ne2)=-Voigt1*Kyye-Voigt3*Kxxe;
elseif nsd==3
  KVpp(ne1,ne1)=-Voigt1*Kxxe-Voigt3*Kyye-Voigt3*Kzze;
  KVpp(ne1,ne2)=-Voigt2*Kxye-Voigt3*Kyxe;
  KVpp(ne1,ne3)=-Voigt2*Kxze-Voigt3*Kzxe;
  KVpp(ne2,ne1)=-Voigt2*Kyxe-Voigt3*Kxye;
  KVpp(ne2,ne2)=-Voigt1*Kyye-Voigt3*Kxxe-Voigt3*Kzze;
  KVpp(ne2,ne3)=-Voigt2*Kyze-Voigt3*Kzye;
  KVpp(ne3,ne1)=-Voigt2*Kzxe-Voigt3*Kxze;
  KVpp(ne3,ne2)=-Voigt2*Kzye-Voigt3*Kyze;
  KVpp(ne3,ne3)=-Voigt1*Kzze-Voigt3*Kxxe-Voigt3*Kyye;
end

KVtp(1,ne1)=Nw1eT*Ne;
KVtp(2,ne2)=Nw1eT*Ne;
if nsd==3
  KVtp(3,ne3)=Nw1eT*Ne;
end

if nsd==2
  KVrp(1,ne1)=-Nw1eT*Ney;
  KVrp(1,ne2)=+Nw1eT*Nex;
elseif nsd==3
  KVrp(1,ne2)=-Nw1eT*Nez;
  KVrp(1,ne3)=+Nw1eT*Ney;
  KVrp(2,ne1)=+Nw1eT*Nez;
  KVrp(2,ne3)=-Nw1eT*Nex;
  KVrp(3,ne1)=-Nw1eT*Ney;
  KVrp(3,ne2)=+Nw1eT*Nex;
end

KTpp(ne1,ne1)=-sqrt(kappa)*Kxxe...
              -sqrt(kappa)*Kyye;
if nsd==3
    KTpp(ne1,ne1)=KTpp(ne1,ne1)-sqrt(kappa)*Kzze;
end

KTtp(1,ne1)=Nw1eT*Ne;

% Compute rhs
if nsd==2
  fVp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg);
  fVp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg);
elseif nsd==3
  fVp(ne1,1)=NwexT*(Lxxeg)+NweyT*(Lxyeg)+NwezT*(Lxzeg);
  fVp(ne2,1)=NwexT*(Lxyeg)+NweyT*(Lyyeg)+NwezT*(Lyzeg);
  fVp(ne3,1)=NwexT*(Lxzeg)+NweyT*(Lyzeg)+NwezT*(Lzzeg);
end

fVt(1,1)=Nw1eT*(wxeg./reg);
fVt(2,1)=Nw1eT*(wyeg./reg);
if nsd==3
  fVt(3,1)=Nw1eT*(wzeg./reg);
end

fTp(ne1,1)=NwexT*(Qxeg)...
          +NweyT*(Qyeg);
if nsd==3
    fTp(ne1,1)=fTp(ne1,1)+NwezT*(Qzeg);
end

fTt(1,1)=Nw1eT*(Teg);

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions(0,RefElement.PostLow,RefElement.Post,Xf,nsd);
  N1f=ones(length(wfg),1);
  
  % Indices
  nlefU1=(iFace-1)*(1+nsd+1)*NumFaceNodes+(1:NumFaceNodes);
  nlefU2=nlefU1+NumFaceNodes;
  nlefU3=nlefU2+NumFaceNodes;
  nlefU4=nlefU3+NumFaceNodes;
  
  % Flip face
  Node2Match1stNode1=Faces.Interior(2,iFace);
  FlipFace=max(Node2Match1stNode1);
  if FlipFace
    order=flipFace(nsd,Parameters.Degree,Node2Match1stNode1);
    nlefU1=nlefU1(order);
    nlefU2=nlefU2(order);
    nlefU3=nlefU3(order);
    nlefU4=nlefU4(order);
  end
  
  % Compute variables at nodes
  Rf=Ue(nlefU1);
  Wxf=Ue(nlefU2);
  Wyf=Ue(nlefU3);
  if nsd==3
    Wzf=Ue(nlefU4);
  end
  
  % Compute variables at Gauss points
  Rfg=Nlf*Rf;
  Wxfg=Nlf*Wxf;
  Wyfg=Nlf*Wyf;
  if nsd==3
    Wzfg=Nlf*Wzf;
  end
  
  % Compute basic matrices
  Nw1fT=(wfg.*N1f)';
  
  % Compute rhs
  if nsd==2
    fVr(1)=fVr(1)+Nw1fT*(-Wxfg./Rfg.*ny+Wyfg./Rfg.*nx);
  elseif nsd==3
    fVr(1)=fVr(1)+Nw1fT*(-Wyfg./Rfg.*nz+Wzfg./Rfg.*ny);
    fVr(2)=fVr(2)+Nw1fT*(+Wxfg./Rfg.*nz-Wzfg./Rfg.*nx);
    fVr(3)=fVr(3)+Nw1fT*(-Wxfg./Rfg.*ny+Wyfg./Rfg.*nx);
  end
end

% Indices
iVp=1:nsd*NumElementNodesPost;
iVt=iVp(end)+(1:nsd);
iVr=iVt(end)+(1:qsd);
iTp=iVr(end)+(1:NumElementNodesPost); iTp2=iVp(end)+(1:NumElementNodesPost);
iTt=iTp(end)+1;

% Initialization of lhs and rhs
LhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd+NumElementNodesPost+1,...
              nsd*NumElementNodesPost+NumElementNodesPost);
RhsPost=zeros(nsd*NumElementNodesPost+nsd+qsd+NumElementNodesPost+1,1);

% Lhs for post-processing
LhsPost(iVp,iVp)=KVpp;
LhsPost(iVt,iVp)=KVtp;
LhsPost(iVr,iVp)=KVrp;
LhsPost(iTp,iTp2)=KTpp;
LhsPost(iTt,iTp2)=KTtp;

% Rhs for post-processing
RhsPost(iVp,1)=fVp;
RhsPost(iVt,1)=fVt;
RhsPost(iVr,1)=fVr;
RhsPost(iTp,1)=fTp;
RhsPost(iTt,1)=fTt;

end