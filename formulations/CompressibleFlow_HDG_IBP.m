classdef CompressibleFlow_HDG_IBP < Formulation  
  
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
KwE=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
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

% Get Gauss data
gwe=RefElement.GaussWeigthsElem;
gwf=RefElement.GaussWeigthsFace;
nge=(1:length(gwe))';
ngf=(1:length(gwf))';

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,~,pinvNe]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);
WEG=sparse(nge,nge,weg);

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
if nsd==2
  peg=p(reg,[wxeg,wyeg],eeg);
  dpdueg=dpdu(reg,[wxeg,wyeg],eeg);
  dpdreg=dpdueg(:,1);
  dpdwxeg=dpdueg(:,2);
  dpdwyeg=dpdueg(:,3);
  dpdeeg=dpdueg(:,4);
  Teg=T(reg,[wxeg,wyeg],eeg);
  dTdueg=dTdu(reg,[wxeg,wyeg],eeg);
  dTdreg=dTdueg(:,1);
  dTdwxeg=dTdueg(:,2);
  dTdwyeg=dTdueg(:,3);
  dTdeeg=dTdueg(:,4);
elseif nsd==3
  peg=p(reg,[wxeg,wyeg,wzeg],eeg);
  dpdueg=dpdu(reg,[wxeg,wyeg,wzeg],eeg);
  dpdreg=dpdueg(:,1);
  dpdwxeg=dpdueg(:,2);
  dpdwyeg=dpdueg(:,3);
  dpdwzeg=dpdueg(:,4);
  dpdeeg=dpdueg(:,5);
  Teg=T(reg,[wxeg,wyeg,wzeg],eeg);
  dTdueg=dTdu(reg,[wxeg,wyeg,wzeg],eeg);
  dTdreg=dTdueg(:,1);
  dTdwxeg=dTdueg(:,2);
  dTdwyeg=dTdueg(:,3);
  dTdwzeg=dTdueg(:,4);
  dTdeeg=dTdueg(:,5);
end
dvxdxeg=Nex*(pinvNe*(wxeg./reg));
dvxdyeg=Ney*(pinvNe*(wxeg./reg));
dvydxeg=Nex*(pinvNe*(wyeg./reg));
dvydyeg=Ney*(pinvNe*(wyeg./reg));
if nsd==3
  dvxdzeg=Nez*(pinvNe*(wxeg./reg));
  dvydzeg=Nez*(pinvNe*(wyeg./reg));
  dvzdxeg=Nex*(pinvNe*(wzeg./reg));
  dvzdyeg=Ney*(pinvNe*(wzeg./reg));
  dvzdzeg=Nez*(pinvNe*(wzeg./reg));
end
vxeg=wxeg./reg;
vyeg=wyeg./reg;
if nsd==3
  vzeg=wzeg./reg;
end
if nsd==2
  sxxeg=-Voigt1*Lxxeg-Voigt2*Lyyeg;
  syyeg=-Voigt1*Lyyeg-Voigt2*Lxxeg;
  sxyeg=-Voigt3*Lxyeg;
elseif nsd==3
  sxxeg=-Voigt1*Lxxeg-Voigt2*Lyyeg-Voigt2*Lzzeg;
  syyeg=-Voigt1*Lyyeg-Voigt2*Lxxeg-Voigt2*Lzzeg;
  szzeg=-Voigt1*Lzzeg-Voigt2*Lxxeg-Voigt2*Lyyeg;
  sxyeg=-Voigt3*Lxyeg;
  sxzeg=-Voigt3*Lxzeg;
  syzeg=-Voigt3*Lyzeg;
end
if nsd==2
  drm2wsdxeg=Nex*(pinvNe*(1./reg.^2.*(wxeg.*sxxeg+wyeg.*sxyeg)));
  drm2wsdyeg=Ney*(pinvNe*(1./reg.^2.*(wxeg.*sxyeg+wyeg.*syyeg)));
elseif nsd==3
  drm2wsdxeg=Nex*(pinvNe*(1./reg.^2.*(wxeg.*sxxeg+wyeg.*sxyeg+wzeg.*sxzeg)));
  drm2wsdyeg=Ney*(pinvNe*(1./reg.^2.*(wxeg.*sxyeg+wyeg.*syyeg+wzeg.*syzeg)));
  drm2wsdzeg=Nez*(pinvNe*(1./reg.^2.*(wxeg.*sxzeg+wyeg.*syzeg+wzeg.*szzeg)));
end
if nsd==2
  rm2wsxeg=1./reg.^2.*(wxeg.*sxxeg+wyeg.*sxyeg);
  rm2wsyeg=1./reg.^2.*(wxeg.*sxyeg+wyeg.*syyeg);
elseif nsd==3
  rm2wsxeg=1./reg.^2.*(wxeg.*sxxeg+wyeg.*sxyeg+wzeg.*sxzeg);
  rm2wsyeg=1./reg.^2.*(wxeg.*sxyeg+wyeg.*syyeg+wzeg.*syzeg);
  rm2wszeg=1./reg.^2.*(wxeg.*sxzeg+wyeg.*syzeg+wzeg.*szzeg);
end
drm1sxxdxeg=Nex*(pinvNe*(1./reg.*sxxeg));
drm1sxydxeg=Nex*(pinvNe*(1./reg.*sxyeg));
drm1sxydyeg=Ney*(pinvNe*(1./reg.*sxyeg));
drm1syydyeg=Ney*(pinvNe*(1./reg.*syyeg));
if nsd==3
  drm1sxzdxeg=Nex*(pinvNe*(1./reg.*sxzeg));
  drm1sxzdzeg=Nez*(pinvNe*(1./reg.*sxzeg));
  drm1syzdyeg=Ney*(pinvNe*(1./reg.*syzeg));
  drm1syzdzeg=Nez*(pinvNe*(1./reg.*syzeg));
  drm1szzdzeg=Nez*(pinvNe*(1./reg.*szzeg));
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
Me=Ne'*(WEG*Ne);
Me=(Me+Me')/2;
Cxe=Nex'*(WEG*Ne);
Cye=Ney'*(WEG*Ne);
if nsd==3
  Cze=Nez'*(WEG*Ne);
end

% Compute common terms
INVDENS=sparse(nge,nge,1./reg);
MOMX=sparse(nge,nge,wxeg);
MOMY=sparse(nge,nge,wyeg);
if nsd==3
  MOMZ=sparse(nge,nge,wzeg);
end
ENER=sparse(nge,nge,eeg);
PRES=sparse(nge,nge,peg);
FORCEX=sparse(nge,nge,fxeg);
FORCEY=sparse(nge,nge,fyeg);
if nsd==3
  FORCEZ=sparse(nge,nge,fzeg);
end
DPRESDDENS=sparse(nge,nge,dpdreg);
DPRESDMOMX=sparse(nge,nge,dpdwxeg);
DPRESDMOMY=sparse(nge,nge,dpdwyeg);
if nsd==3
  DPRESDMOMZ=sparse(nge,nge,dpdwzeg);
end
DPRESDENER=sparse(nge,nge,dpdeeg);
DTEMPDDENS=sparse(nge,nge,dTdreg);
DTEMPDMOMX=sparse(nge,nge,dTdwxeg);
DTEMPDMOMY=sparse(nge,nge,dTdwyeg);
if nsd==3
  DTEMPDMOMZ=sparse(nge,nge,dTdwzeg);
end
DTEMPDENER=sparse(nge,nge,dTdeeg);
DVELXDX=sparse(nge,nge,dvxdxeg);
DVELXDY=sparse(nge,nge,dvxdyeg);
DVELYDX=sparse(nge,nge,dvydxeg);
DVELYDY=sparse(nge,nge,dvydyeg);
if nsd==3
  DVELXDZ=sparse(nge,nge,dvxdzeg);
  DVELYDZ=sparse(nge,nge,dvydzeg);
  DVELZDX=sparse(nge,nge,dvzdxeg);
  DVELZDY=sparse(nge,nge,dvzdyeg);
  DVELZDZ=sparse(nge,nge,dvzdzeg);
end
VELX=sparse(nge,nge,vxeg);
VELY=sparse(nge,nge,vyeg);
if nsd==3
  VELZ=sparse(nge,nge,vzeg);
end
DINVDENS2MOMSTRESSDX=sparse(nge,nge,drm2wsdxeg);
DINVDENS2MOMSTRESSDY=sparse(nge,nge,drm2wsdyeg);
if nsd==3
  DINVDENS2MOMSTRESSDZ=sparse(nge,nge,drm2wsdzeg);
end
INVDENS2MOMSTRESSX=sparse(nge,nge,rm2wsxeg);
INVDENS2MOMSTRESSY=sparse(nge,nge,rm2wsyeg);
if nsd==3
  INVDENS2MOMSTRESSZ=sparse(nge,nge,rm2wszeg);
end
DINVDENSSTRESSXXDX=sparse(nge,nge,drm1sxxdxeg);
DINVDENSSTRESSXYDX=sparse(nge,nge,drm1sxydxeg);
DINVDENSSTRESSXYDY=sparse(nge,nge,drm1sxydyeg);
DINVDENSSTRESSYYDY=sparse(nge,nge,drm1syydyeg);
if nsd==3
  DINVDENSSTRESSXZDX=sparse(nge,nge,drm1sxzdxeg);
  DINVDENSSTRESSXZDZ=sparse(nge,nge,drm1sxzdzeg);
  DINVDENSSTRESSYZDY=sparse(nge,nge,drm1syzdyeg);
  DINVDENSSTRESSYZDZ=sparse(nge,nge,drm1syzdzeg);
  DINVDENSSTRESSZZDZ=sparse(nge,nge,drm1szzdzeg);
end
INVDENSSTRESSXX=sparse(nge,nge,rm1sxxeg);
INVDENSSTRESSXY=sparse(nge,nge,rm1sxyeg);
INVDENSSTRESSYY=sparse(nge,nge,rm1syyeg);
if nsd==3
  INVDENSSTRESSXZ=sparse(nge,nge,rm1sxzeg);
  INVDENSSTRESSYZ=sparse(nge,nge,rm1syzeg);
  INVDENSSTRESSZZ=sparse(nge,nge,rm1szzeg);
end
MomInvdens2_x=MOMX*(INVDENS.^2)*WEG;
MomInvdens2_y=MOMY*(INVDENS.^2)*WEG;
if nsd==3
  MomInvdens2_z=MOMZ*(INVDENS.^2)*WEG;
end
Invdens=INVDENS*WEG;
dTempdDens=DTEMPDDENS*WEG;
dTempdMomx=DTEMPDMOMX*WEG;
dTempdMomy=DTEMPDMOMY*WEG;
if nsd==3
  dTempdMomz=DTEMPDMOMZ*WEG;
end
dTempdEner=DTEMPDENER*WEG;
if nsd==2
  MomInvdens2_xx=MOMX*MOMX*(INVDENS.^2)*WEG;
  MomInvdens2_xy=MOMX*MOMY*(INVDENS.^2)*WEG;
  MomInvdens2_yx=MOMY*MOMX*(INVDENS.^2)*WEG;
  MomInvdens2_yy=MOMY*MOMY*(INVDENS.^2)*WEG;
elseif nsd==3
  MomInvdens2_xx=MOMX*MOMX*(INVDENS.^2)*WEG;
  MomInvdens2_xy=MOMX*MOMY*(INVDENS.^2)*WEG;
  MomInvdens2_xz=MOMX*MOMZ*(INVDENS.^2)*WEG;
  MomInvdens2_yx=MOMY*MOMX*(INVDENS.^2)*WEG;
  MomInvdens2_yy=MOMY*MOMY*(INVDENS.^2)*WEG;
  MomInvdens2_yz=MOMY*MOMZ*(INVDENS.^2)*WEG;
  MomInvdens2_zx=MOMZ*MOMX*(INVDENS.^2)*WEG;
  MomInvdens2_zy=MOMZ*MOMY*(INVDENS.^2)*WEG;
  MomInvdens2_zz=MOMZ*MOMZ*(INVDENS.^2)*WEG;
end
if nsd==2
  Conv_xx1=(2*MOMX*INVDENS)*WEG;
  Conv_xx2=(  MOMY*INVDENS)*WEG;
  Conv_xy =(  MOMX*INVDENS)*WEG;
  Conv_yx =(  MOMY*INVDENS)*WEG;
  Conv_yy1=(  MOMX*INVDENS)*WEG;
  Conv_yy2=(2*MOMY*INVDENS)*WEG;
elseif nsd==3
  Conv_xx1=(2*MOMX*INVDENS)*WEG;
  Conv_xx2=(  MOMY*INVDENS)*WEG;
  Conv_xx3=(  MOMZ*INVDENS)*WEG;
  Conv_xy =(  MOMX*INVDENS)*WEG;
  Conv_xz =(  MOMX*INVDENS)*WEG;
  Conv_yx =(  MOMY*INVDENS)*WEG;
  Conv_yy1=(  MOMX*INVDENS)*WEG;
  Conv_yy2=(2*MOMY*INVDENS)*WEG;
  Conv_yy3=(  MOMZ*INVDENS)*WEG;
  Conv_yz =(  MOMY*INVDENS)*WEG;
  Conv_zx =(  MOMZ*INVDENS)*WEG;
  Conv_zy =(  MOMZ*INVDENS)*WEG;
  Conv_zz1=(  MOMX*INVDENS)*WEG;
  Conv_zz2=(  MOMY*INVDENS)*WEG;
  Conv_zz3=(2*MOMZ*INVDENS)*WEG;
end
dPresdDens=DPRESDDENS*WEG;
dPresdMomx=DPRESDMOMX*WEG;
dPresdMomy=DPRESDMOMY*WEG;
if nsd==3
  dPresdMomz=DPRESDMOMZ*WEG;
end
dPresdEner=DPRESDENER*WEG;
EnerMomxInvdens2=ENER*MOMX*(INVDENS.^2)*WEG;
EnerMomyInvdens2=ENER*MOMY*(INVDENS.^2)*WEG;
if nsd==3
  EnerMomzInvdens2=ENER*MOMZ*(INVDENS.^2)*WEG;
end
EnerInvdens=ENER*INVDENS*WEG;
MomxInvdens=MOMX*INVDENS*WEG;
MomyInvdens=MOMY*INVDENS*WEG;
if nsd==3
  MomzInvdens=MOMZ*INVDENS*WEG;
end
dVelxdx=DVELXDX*WEG;
dVelxdy=DVELXDY*WEG;
dVelydx=DVELYDX*WEG;
dVelydy=DVELYDY*WEG;
if nsd==3
  dVelxdz=DVELXDZ*WEG;
  dVelydz=DVELYDZ*WEG;
  dVelzdx=DVELZDX*WEG;
  dVelzdy=DVELZDY*WEG;
  dVelzdz=DVELZDZ*WEG;
end
Velx=VELX*WEG;
Vely=VELY*WEG;
if nsd==3
  Velz=VELZ*WEG;
end
dInvdens2MomStressdx=DINVDENS2MOMSTRESSDX*WEG;
dInvdens2MomStressdy=DINVDENS2MOMSTRESSDY*WEG;
if nsd==3
  dInvdens2MomStressdz=DINVDENS2MOMSTRESSDZ*WEG;
end
Invdens2MomStress_x=INVDENS2MOMSTRESSX*WEG;
Invdens2MomStress_y=INVDENS2MOMSTRESSY*WEG;
if nsd==3
  Invdens2MomStress_z=INVDENS2MOMSTRESSZ*WEG;
end
dInvdensStressxxdx=DINVDENSSTRESSXXDX*WEG;
dInvdensStressxydx=DINVDENSSTRESSXYDX*WEG;
dInvdensStressxydy=DINVDENSSTRESSXYDY*WEG;
dInvdensStressyydy=DINVDENSSTRESSYYDY*WEG;
if nsd==3
  dInvdensStressxzdx=DINVDENSSTRESSXZDX*WEG;
  dInvdensStressyzdy=DINVDENSSTRESSYZDY*WEG;
  dInvdensStressxzdz=DINVDENSSTRESSXZDZ*WEG;
  dInvdensStressyzdz=DINVDENSSTRESSYZDZ*WEG;
  dInvdensStresszzdz=DINVDENSSTRESSZZDZ*WEG;
end
InvdensStressxx=INVDENSSTRESSXX*WEG;
InvdensStressxy=INVDENSSTRESSXY*WEG;
InvdensStressyy=INVDENSSTRESSYY*WEG;
if nsd==3
  InvdensStressxz=INVDENSSTRESSXZ*WEG;
  InvdensStressyz=INVDENSSTRESSYZ*WEG;
  InvdensStresszz=INVDENSSTRESSZZ*WEG;
end
dPresdDensMomxInvdens=DPRESDDENS*MOMX*INVDENS*WEG;
dPresdDensMomyInvdens=DPRESDDENS*MOMY*INVDENS*WEG;
if nsd==3
  dPresdDensMomzInvdens=DPRESDDENS*MOMZ*INVDENS*WEG;
end
PresMomxInvdens2=PRES*MOMX*(INVDENS.^2)*WEG;
PresMomyInvdens2=PRES*MOMY*(INVDENS.^2)*WEG;
if nsd==3
  PresMomzInvdens2=PRES*MOMZ*(INVDENS.^2)*WEG;
end
dPresdMomxMomxInvdens=DPRESDMOMX*MOMX*INVDENS*WEG;
dPresdMomxMomyInvdens=DPRESDMOMX*MOMY*INVDENS*WEG;
dPresdMomyMomxInvdens=DPRESDMOMY*MOMX*INVDENS*WEG;
dPresdMomyMomyInvdens=DPRESDMOMY*MOMY*INVDENS*WEG;
if nsd==3
  dPresdMomxMomzInvdens=DPRESDMOMX*MOMZ*INVDENS*WEG;
  dPresdMomyMomzInvdens=DPRESDMOMY*MOMZ*INVDENS*WEG;
  dPresdMomzMomxInvdens=DPRESDMOMZ*MOMX*INVDENS*WEG;
  dPresdMomzMomyInvdens=DPRESDMOMZ*MOMY*INVDENS*WEG;
  dPresdMomzMomzInvdens=DPRESDMOMZ*MOMZ*INVDENS*WEG;
end
PresInvdens=PRES*INVDENS*WEG;
dPresdEnerMomxInvdens=DPRESDENER*MOMX*INVDENS*WEG;
dPresdEnerMomyInvdens=DPRESDENER*MOMY*INVDENS*WEG;
if nsd==3
  dPresdEnerMomzInvdens=DPRESDENER*MOMZ*INVDENS*WEG;
end
if nsd==2
  ForceMomInvdens2=(FORCEX*MOMX+FORCEY*MOMY)*(INVDENS.^2)*WEG;
elseif nsd==3
  ForceMomInvdens2=(FORCEX*MOMX+FORCEY*MOMY+FORCEZ*MOMZ)*(INVDENS.^2)*WEG;
end
ForcexInvdens=FORCEX*INVDENS*WEG;
ForceyInvdens=FORCEY*INVDENS*WEG;
if nsd==3
  ForcezInvdens=FORCEZ*INVDENS*WEG;
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
  KLr(ne1,ne1)=-Voigt1*Nex'*(MomInvdens2_x*Ne)...
               -Voigt2*Ney'*(MomInvdens2_y*Ne);
  KLr(ne2,ne1)=-Voigt2*Nex'*(MomInvdens2_x*Ne)...
               -Voigt1*Ney'*(MomInvdens2_y*Ne);
  KLr(ne3,ne1)=-Voigt3*Ney'*(MomInvdens2_x*Ne)...
               -Voigt3*Nex'*(MomInvdens2_y*Ne);
elseif nsd==3
  KLr(ne1,ne1)=-Voigt1*Nex'*(MomInvdens2_x*Ne)...
               -Voigt2*Ney'*(MomInvdens2_y*Ne)...
               -Voigt2*Nez'*(MomInvdens2_z*Ne);
  KLr(ne2,ne1)=-Voigt2*Nex'*(MomInvdens2_x*Ne)...
               -Voigt1*Ney'*(MomInvdens2_y*Ne)...
               -Voigt2*Nez'*(MomInvdens2_z*Ne);
  KLr(ne3,ne1)=-Voigt2*Nex'*(MomInvdens2_x*Ne)...
               -Voigt2*Ney'*(MomInvdens2_y*Ne)...
               -Voigt1*Nez'*(MomInvdens2_z*Ne);
  KLr(ne4,ne1)=-Voigt3*Ney'*(MomInvdens2_x*Ne)...
               -Voigt3*Nex'*(MomInvdens2_y*Ne);
  KLr(ne5,ne1)=-Voigt3*Nez'*(MomInvdens2_x*Ne)...
               -Voigt3*Nex'*(MomInvdens2_z*Ne);
  KLr(ne6,ne1)=-Voigt3*Nez'*(MomInvdens2_y*Ne)...
               -Voigt3*Ney'*(MomInvdens2_z*Ne);
end

if nsd==2
  KLw(ne1,ne1)=Voigt1*Nex'*(Invdens*Ne);
  KLw(ne2,ne1)=Voigt2*Nex'*(Invdens*Ne);
  KLw(ne3,ne1)=Voigt3*Ney'*(Invdens*Ne);
  KLw(ne1,ne2)=Voigt2*Ney'*(Invdens*Ne);
  KLw(ne2,ne2)=Voigt1*Ney'*(Invdens*Ne);
  KLw(ne3,ne2)=Voigt3*Nex'*(Invdens*Ne);
elseif nsd==3
  KLw(ne1,ne1)=Voigt1*Nex'*(Invdens*Ne);
  KLw(ne2,ne1)=Voigt2*Nex'*(Invdens*Ne);
  KLw(ne3,ne1)=Voigt2*Nex'*(Invdens*Ne);
  KLw(ne4,ne1)=Voigt3*Ney'*(Invdens*Ne);
  KLw(ne5,ne1)=Voigt3*Nez'*(Invdens*Ne);
  KLw(ne1,ne2)=Voigt2*Ney'*(Invdens*Ne);
  KLw(ne2,ne2)=Voigt1*Ney'*(Invdens*Ne);
  KLw(ne3,ne2)=Voigt2*Ney'*(Invdens*Ne);
  KLw(ne4,ne2)=Voigt3*Nex'*(Invdens*Ne);
  KLw(ne6,ne2)=Voigt3*Nez'*(Invdens*Ne);
  KLw(ne1,ne3)=Voigt2*Nez'*(Invdens*Ne);
  KLw(ne2,ne3)=Voigt2*Nez'*(Invdens*Ne);
  KLw(ne3,ne3)=Voigt1*Nez'*(Invdens*Ne);
  KLw(ne5,ne3)=Voigt3*Nex'*(Invdens*Ne);
  KLw(ne6,ne3)=Voigt3*Ney'*(Invdens*Ne);
end

KQQ(ne1,ne1)=-Me;
KQQ(ne2,ne2)=-Me;
if nsd==3
  KQQ(ne3,ne3)=-Me;
end

KQr(ne1,ne1)=sqrt(kappa)*Nex'*(dTempdDens*Ne);
KQr(ne2,ne1)=sqrt(kappa)*Ney'*(dTempdDens*Ne);
if nsd==3
  KQr(ne3,ne1)=sqrt(kappa)*Nez'*(dTempdDens*Ne);
end

if nsd==2
  KQw(ne1,ne1)=sqrt(kappa)*Nex'*(dTempdMomx*Ne);
  KQw(ne2,ne1)=sqrt(kappa)*Ney'*(dTempdMomx*Ne);
  KQw(ne1,ne2)=sqrt(kappa)*Nex'*(dTempdMomy*Ne);
  KQw(ne2,ne2)=sqrt(kappa)*Ney'*(dTempdMomy*Ne);
elseif nsd==3
  KQw(ne1,ne1)=sqrt(kappa)*Nex'*(dTempdMomx*Ne);
  KQw(ne2,ne1)=sqrt(kappa)*Ney'*(dTempdMomx*Ne);
  KQw(ne3,ne1)=sqrt(kappa)*Nez'*(dTempdMomx*Ne);
  KQw(ne1,ne2)=sqrt(kappa)*Nex'*(dTempdMomy*Ne);
  KQw(ne2,ne2)=sqrt(kappa)*Ney'*(dTempdMomy*Ne);
  KQw(ne3,ne2)=sqrt(kappa)*Nez'*(dTempdMomy*Ne);
  KQw(ne1,ne3)=sqrt(kappa)*Nex'*(dTempdMomz*Ne);
  KQw(ne2,ne3)=sqrt(kappa)*Ney'*(dTempdMomz*Ne);
  KQw(ne3,ne3)=sqrt(kappa)*Nez'*(dTempdMomz*Ne);
end

KQe(ne1,ne1)=sqrt(kappa)*Nex'*(dTempdEner*Ne);
KQe(ne2,ne1)=sqrt(kappa)*Ney'*(dTempdEner*Ne);
if nsd==3
  KQe(ne3,ne1)=sqrt(kappa)*Nez'*(dTempdEner*Ne);
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
    Kwr(ne1,ne1)=+Nex'*(MomInvdens2_xx*Ne)...
                 +Ney'*(MomInvdens2_xy*Ne);
    Kwr(ne2,ne1)=+Nex'*(MomInvdens2_yx*Ne)...
                 +Ney'*(MomInvdens2_yy*Ne);
  elseif nsd==3
    Kwr(ne1,ne1)=+Nex'*(MomInvdens2_xx*Ne)...
                 +Ney'*(MomInvdens2_xy*Ne)...
                 +Nez'*(MomInvdens2_xz*Ne);
    Kwr(ne2,ne1)=+Nex'*(MomInvdens2_yx*Ne)...
                 +Ney'*(MomInvdens2_yy*Ne)...
                 +Nez'*(MomInvdens2_yz*Ne);
    Kwr(ne3,ne1)=+Nex'*(MomInvdens2_zx*Ne)...
                 +Ney'*(MomInvdens2_zy*Ne)...
                 +Nez'*(MomInvdens2_zz*Ne);
  end
  
  if nsd==2
    Kww(ne1,ne1)=Kww(ne1,ne1)-Nex'*(Conv_xx1*Ne)...
                             -Ney'*(Conv_xx2*Ne);
    Kww(ne2,ne1)=Kww(ne2,ne1)-Nex'*(Conv_yx *Ne);
    Kww(ne1,ne2)=Kww(ne1,ne2)-Ney'*(Conv_xy *Ne);
    Kww(ne2,ne2)=Kww(ne2,ne2)-Nex'*(Conv_yy1*Ne)...
                             -Ney'*(Conv_yy2*Ne);
  elseif nsd==3
    Kww(ne1,ne1)=Kww(ne1,ne1)-Nex'*(Conv_xx1*Ne)...
                             -Ney'*(Conv_xx2*Ne)...
                             -Nez'*(Conv_xx3*Ne);
    Kww(ne2,ne1)=Kww(ne2,ne1)-Nex'*(Conv_yx *Ne);
    Kww(ne3,ne1)=Kww(ne3,ne1)-Nez'*(Conv_zx *Ne);
    Kww(ne1,ne2)=Kww(ne1,ne2)-Ney'*(Conv_xy *Ne);
    Kww(ne2,ne2)=Kww(ne2,ne2)-Nex'*(Conv_yy1*Ne)...
                             -Ney'*(Conv_yy2*Ne)...
                             -Nez'*(Conv_yy3*Ne);
    Kww(ne3,ne2)=Kww(ne3,ne2)-Nez'*(Conv_zy *Ne);
    Kww(ne1,ne3)=Kww(ne1,ne3)-Nez'*(Conv_xz *Ne);
    Kww(ne2,ne3)=Kww(ne2,ne3)-Nez'*(Conv_yz *Ne);
    Kww(ne3,ne3)=Kww(ne3,ne3)-Nex'*(Conv_zz1*Ne)...
                             -Ney'*(Conv_zz2*Ne)...
                             -Nez'*(Conv_zz3*Ne);
  end
end

if nsd==2
  KwL(ne1,ne1)=Voigt1*Cxe';
  KwL(ne2,ne1)=Voigt2*Cye';
  KwL(ne1,ne2)=Voigt2*Cxe';
  KwL(ne2,ne2)=Voigt1*Cye';
  KwL(ne1,ne3)=Voigt3*Cye';
  KwL(ne2,ne3)=Voigt3*Cxe';
elseif nsd==3
  KwL(ne1,ne1)=Voigt1*Cxe';
  KwL(ne2,ne1)=Voigt2*Cye';
  KwL(ne3,ne1)=Voigt2*Cze';
  KwL(ne1,ne2)=Voigt2*Cxe';
  KwL(ne2,ne2)=Voigt1*Cye';
  KwL(ne3,ne2)=Voigt2*Cze';
  KwL(ne1,ne3)=Voigt2*Cxe';
  KwL(ne2,ne3)=Voigt2*Cye';
  KwL(ne3,ne3)=Voigt1*Cze';
  KwL(ne1,ne4)=Voigt3*Cye';
  KwL(ne2,ne4)=Voigt3*Cxe';
  KwL(ne1,ne5)=Voigt3*Cze';
  KwL(ne3,ne5)=Voigt3*Cxe';
  KwL(ne2,ne6)=Voigt3*Cze';
  KwL(ne3,ne6)=Voigt3*Cye';
end

Kwr(ne1,ne1)=Kwr(ne1,ne1)-Nex'*(dPresdDens*Ne);
Kwr(ne2,ne1)=Kwr(ne2,ne1)-Ney'*(dPresdDens*Ne);
if nsd==3
  Kwr(ne3,ne1)=Kwr(ne3,ne1)-Nez'*(dPresdDens*Ne);
end

if nsd==2
  Kww(ne1,ne1)=Kww(ne1,ne1)-Nex'*(dPresdMomx*Ne);
  Kww(ne2,ne1)=Kww(ne2,ne1)-Ney'*(dPresdMomx*Ne);
  Kww(ne1,ne2)=Kww(ne1,ne2)-Nex'*(dPresdMomy*Ne);
  Kww(ne2,ne2)=Kww(ne2,ne2)-Ney'*(dPresdMomy*Ne);
elseif nsd==3
  Kww(ne1,ne1)=Kww(ne1,ne1)-Nex'*(dPresdMomx*Ne);
  Kww(ne2,ne1)=Kww(ne2,ne1)-Ney'*(dPresdMomx*Ne);
  Kww(ne3,ne1)=Kww(ne3,ne1)-Nez'*(dPresdMomx*Ne);
  Kww(ne1,ne2)=Kww(ne1,ne2)-Nex'*(dPresdMomy*Ne);
  Kww(ne2,ne2)=Kww(ne2,ne2)-Ney'*(dPresdMomy*Ne);
  Kww(ne3,ne2)=Kww(ne3,ne2)-Nez'*(dPresdMomy*Ne);
  Kww(ne1,ne3)=Kww(ne1,ne3)-Nex'*(dPresdMomz*Ne);
  Kww(ne2,ne3)=Kww(ne2,ne3)-Ney'*(dPresdMomz*Ne);
  Kww(ne3,ne3)=Kww(ne3,ne3)-Nez'*(dPresdMomz*Ne);
end

Kwe(ne1,ne1)=-Nex'*(dPresdEner*Ne);
Kwe(ne2,ne1)=-Ney'*(dPresdEner*Ne);
if nsd==3
  Kwe(ne3,ne1)=-Nez'*(dPresdEner*Ne);
end

if isTimeDependent
  Kee(ne1,ne1)=alpha(1)/dt*Me;
end

Ker(ne1,ne1)=+Nex'*(EnerMomxInvdens2*Ne)...
             +Ney'*(EnerMomyInvdens2*Ne);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+Nez'*(EnerMomzInvdens2*Ne);
end

Kew(ne1,ne1)=-Nex'*(EnerInvdens*Ne);
Kew(ne1,ne2)=-Ney'*(EnerInvdens*Ne);
if nsd==3
  Kew(ne1,ne3)=-Nez'*(EnerInvdens*Ne);
end

Kee(ne1,ne1)=Kee(ne1,ne1)-Nex'*(MomxInvdens*Ne)...
                         -Ney'*(MomyInvdens*Ne);
if nsd==3
  Kee(ne1,ne1)=Kee(ne1,ne1)-Nez'*(MomzInvdens*Ne);
end

if nsd==2
  KeL(ne1,ne1)=+Voigt1*Ne'*(dVelxdx*Ne)...
               +Voigt2*Ne'*(dVelydy*Ne);
  KeL(ne1,ne2)=+Voigt2*Ne'*(dVelxdx*Ne)...
               +Voigt1*Ne'*(dVelydy*Ne);
  KeL(ne1,ne3)=+Voigt3*Ne'*(dVelxdy*Ne)...
               +Voigt3*Ne'*(dVelydx*Ne);
elseif nsd==3
  KeL(ne1,ne1)=+Voigt1*Ne'*(dVelxdx*Ne)...
               +Voigt2*Ne'*(dVelydy*Ne)...
               +Voigt2*Ne'*(dVelzdz*Ne);
  KeL(ne1,ne2)=+Voigt2*Ne'*(dVelxdx*Ne)...
               +Voigt1*Ne'*(dVelydy*Ne)...
               +Voigt2*Ne'*(dVelzdz*Ne);
  KeL(ne1,ne3)=+Voigt2*Ne'*(dVelxdx*Ne)...
               +Voigt2*Ne'*(dVelydy*Ne)...
               +Voigt1*Ne'*(dVelzdz*Ne);
  KeL(ne1,ne4)=+Voigt3*Ne'*(dVelxdy*Ne)...
               +Voigt3*Ne'*(dVelydx*Ne);
  KeL(ne1,ne5)=+Voigt3*Ne'*(dVelxdz*Ne)...
               +Voigt3*Ne'*(dVelzdx*Ne);
  KeL(ne1,ne6)=+Voigt3*Ne'*(dVelydz*Ne)...
               +Voigt3*Ne'*(dVelzdy*Ne);
end

if nsd==2
  KeL(ne1,ne1)=KeL(ne1,ne1)+Voigt1*Ne'*(Velx*Nex)...
                           +Voigt2*Ne'*(Vely*Ney);
  KeL(ne1,ne2)=KeL(ne1,ne2)+Voigt2*Ne'*(Velx*Nex)...
                           +Voigt1*Ne'*(Vely*Ney);
  KeL(ne1,ne3)=KeL(ne1,ne3)+Voigt3*Ne'*(Vely*Nex)...
                           +Voigt3*Ne'*(Velx*Ney);
elseif nsd==3
  KeL(ne1,ne1)=KeL(ne1,ne1)+Voigt1*Ne'*(Velx*Nex)...
                           +Voigt2*Ne'*(Vely*Ney)...
                           +Voigt2*Ne'*(Velz*Nez);
  KeL(ne1,ne2)=KeL(ne1,ne2)+Voigt2*Ne'*(Velx*Nex)...
                           +Voigt1*Ne'*(Vely*Ney)...
                           +Voigt2*Ne'*(Velz*Nez);
  KeL(ne1,ne3)=KeL(ne1,ne3)+Voigt2*Ne'*(Velx*Nex)...
                           +Voigt2*Ne'*(Vely*Ney)...
                           +Voigt1*Ne'*(Velz*Nez);
  KeL(ne1,ne4)=KeL(ne1,ne4)+Voigt3*Ne'*(Vely*Nex)...
                           +Voigt3*Ne'*(Velx*Ney);
  KeL(ne1,ne5)=KeL(ne1,ne5)+Voigt3*Ne'*(Velz*Nex)...
                           +Voigt3*Ne'*(Velx*Nez);
  KeL(ne1,ne6)=KeL(ne1,ne6)+Voigt3*Ne'*(Velz*Ney)...
                           +Voigt3*Ne'*(Vely*Nez);
end

Ker(ne1,ne1)=Ker(ne1,ne1)+Ne'*(dInvdens2MomStressdx*Ne)...
                         +Ne'*(dInvdens2MomStressdy*Ne);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+Ne'*(dInvdens2MomStressdz*Ne);
end

Ker(ne1,ne1)=Ker(ne1,ne1)+Ne'*(Invdens2MomStress_x*Nex)...
                         +Ne'*(Invdens2MomStress_y*Ney);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+Ne'*(Invdens2MomStress_z*Nez);
end

if nsd==2
  Kew(ne1,ne1)=Kew(ne1,ne1)-Ne'*(dInvdensStressxxdx*Ne)...
                           -Ne'*(dInvdensStressxydy*Ne);
  Kew(ne1,ne2)=Kew(ne1,ne2)-Ne'*(dInvdensStressxydx*Ne)...
                           -Ne'*(dInvdensStressyydy*Ne);
elseif nsd==3
  Kew(ne1,ne1)=Kew(ne1,ne1)-Ne'*(dInvdensStressxxdx*Ne)...
                           -Ne'*(dInvdensStressxydy*Ne)...
                           -Ne'*(dInvdensStressxzdz*Ne);
  Kew(ne1,ne2)=Kew(ne1,ne2)-Ne'*(dInvdensStressxydx*Ne)...
                           -Ne'*(dInvdensStressyydy*Ne)...
                           -Ne'*(dInvdensStressyzdz*Ne);
  Kew(ne1,ne3)=Kew(ne1,ne3)-Ne'*(dInvdensStressxzdx*Ne)...
                           -Ne'*(dInvdensStressyzdy*Ne)...
                           -Ne'*(dInvdensStresszzdz*Ne);
end

if nsd==2
  Kew(ne1,ne1)=Kew(ne1,ne1)-Ne'*(InvdensStressxx*Nex)...
                           -Ne'*(InvdensStressxy*Ney);
  Kew(ne1,ne2)=Kew(ne1,ne2)-Ne'*(InvdensStressxy*Nex)...
                           -Ne'*(InvdensStressyy*Ney);
elseif nsd==3
  Kew(ne1,ne1)=Kew(ne1,ne1)-Ne'*(InvdensStressxx*Nex)...
                           -Ne'*(InvdensStressxy*Ney)...
                           -Ne'*(InvdensStressxz*Nez);
  Kew(ne1,ne2)=Kew(ne1,ne2)-Ne'*(InvdensStressxy*Nex)...
                           -Ne'*(InvdensStressyy*Ney)...
                           -Ne'*(InvdensStressyz*Nez);
  Kew(ne1,ne3)=Kew(ne1,ne3)-Ne'*(InvdensStressxz*Nex)...
                           -Ne'*(InvdensStressyz*Ney)...
                           -Ne'*(InvdensStresszz*Nez);
end

Ker(ne1,ne1)=Ker(ne1,ne1)-Nex'*(dPresdDensMomxInvdens*Ne)...
                         -Ney'*(dPresdDensMomyInvdens*Ne);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)-Nez'*(dPresdDensMomzInvdens*Ne);
end

Ker(ne1,ne1)=Ker(ne1,ne1)+Nex'*(PresMomxInvdens2*Ne)...
                         +Ney'*(PresMomyInvdens2*Ne);
if nsd==3
  Ker(ne1,ne1)=Ker(ne1,ne1)+Nez'*(PresMomzInvdens2*Ne);
end

if nsd==2
  Kew(ne1,ne1)=Kew(ne1,ne1)-Nex'*(dPresdMomxMomxInvdens*Ne)...
                           -Ney'*(dPresdMomxMomyInvdens*Ne);
  Kew(ne1,ne2)=Kew(ne1,ne2)-Nex'*(dPresdMomyMomxInvdens*Ne)...
                           -Ney'*(dPresdMomyMomyInvdens*Ne);  
elseif nsd==3
  Kew(ne1,ne1)=Kew(ne1,ne1)-Nex'*(dPresdMomxMomxInvdens*Ne)...
                           -Ney'*(dPresdMomxMomyInvdens*Ne)...
                           -Nez'*(dPresdMomxMomzInvdens*Ne);
  Kew(ne1,ne2)=Kew(ne1,ne2)-Nex'*(dPresdMomyMomxInvdens*Ne)...
                           -Ney'*(dPresdMomyMomyInvdens*Ne)...
                           -Nez'*(dPresdMomyMomzInvdens*Ne);
  Kew(ne1,ne3)=Kew(ne1,ne3)-Nex'*(dPresdMomzMomxInvdens*Ne)...
                           -Ney'*(dPresdMomzMomyInvdens*Ne)...
                           -Nez'*(dPresdMomzMomzInvdens*Ne);
end

Kew(ne1,ne1)=Kew(ne1,ne1)-Nex'*(PresInvdens*Ne);
Kew(ne1,ne2)=Kew(ne1,ne2)-Ney'*(PresInvdens*Ne);
if nsd==3
  Kew(ne1,ne3)=Kew(ne1,ne3)-Nez'*(PresInvdens*Ne);
end

Kee(ne1,ne1)=Kee(ne1,ne1)-Nex'*(dPresdEnerMomxInvdens*Ne)...
                         -Ney'*(dPresdEnerMomyInvdens*Ne);
if nsd==3
  Kee(ne1,ne1)=Kee(ne1,ne1)-Nez'*(dPresdEnerMomzInvdens*Ne);
end

KeQ(ne1,ne1)=sqrt(kappa)*Cxe';
KeQ(ne1,ne2)=sqrt(kappa)*Cye';
if nsd==3
  KeQ(ne1,ne3)=sqrt(kappa)*Cze';
end

Ker(ne1,ne1)=Ker(ne1,ne1)+Ne'*(ForceMomInvdens2*Ne);

Kew(ne1,ne1)=Kew(ne1,ne1)-Ne'*(ForcexInvdens*Ne);
Kew(ne1,ne2)=Kew(ne1,ne2)-Ne'*(ForceyInvdens*Ne);
if nsd==3
  Kew(ne1,ne3)=Kew(ne1,ne3)-Ne'*(ForcezInvdens*Ne);
end

% Compute rhs
if nsd==2
  fL(ne1,1)=+Ne'*(Lxxeg.*weg)...
            -Nex'*(Voigt1*wxeg./reg.*weg)...
            -Ney'*(Voigt2*wyeg./reg.*weg);
  fL(ne2,1)=+Ne'*(Lyyeg.*weg)...
            -Nex'*(Voigt2*wxeg./reg.*weg)...
            -Ney'*(Voigt1*wyeg./reg.*weg);
  fL(ne3,1)=+Ne'*(Lxyeg.*weg)...
            -Ney'*(Voigt3*wxeg./reg.*weg)...
            -Nex'*(Voigt3*wyeg./reg.*weg);
elseif nsd==3
  fL(ne1,1)=+Ne'*(Lxxeg.*weg)...
            -Nex'*(Voigt1*wxeg./reg.*weg)...
            -Ney'*(Voigt2*wyeg./reg.*weg)...
            -Nez'*(Voigt2*wzeg./reg.*weg);
  fL(ne2,1)=+Ne'*(Lyyeg.*weg)...
            -Nex'*(Voigt2*wxeg./reg.*weg)...
            -Ney'*(Voigt1*wyeg./reg.*weg)...
            -Nez'*(Voigt2*wzeg./reg.*weg);
  fL(ne3,1)=+Ne'*(Lzzeg.*weg)...
            -Nex'*(Voigt2*wxeg./reg.*weg)...
            -Ney'*(Voigt2*wyeg./reg.*weg)...
            -Nez'*(Voigt1*wzeg./reg.*weg);
  fL(ne4,1)=+Ne'*(Lxyeg.*weg)...
            -Nex'*(Voigt3*wyeg./reg.*weg)...
            -Ney'*(Voigt3*wxeg./reg.*weg);
  fL(ne5,1)=+Ne'*(Lxzeg.*weg)...
            -Nex'*(Voigt3*wzeg./reg.*weg)...
            -Nez'*(Voigt3*wxeg./reg.*weg);
  fL(ne6,1)=+Ne'*(Lyzeg.*weg)...
            -Ney'*(Voigt3*wzeg./reg.*weg)...
            -Nez'*(Voigt3*wyeg./reg.*weg);
end

fQ(ne1,1)=+Ne'*(Qxeg.*weg)...
          -Nex'*(sqrt(kappa)*Teg.*weg);
fQ(ne2,1)=+Ne'*(Qyeg.*weg)...
          -Ney'*(sqrt(kappa)*Teg.*weg);
if nsd==3
  fQ(ne3,1)=+Ne'*(Qzeg.*weg)...
            -Nez'*(sqrt(kappa)*Teg.*weg);
end

if isTimeDependent
  fr(ne1,1)=-Ne'*((1/dt*reg*alpha(1)...
                  +1/dt*roldeg*alpha(2:BDFo+1,1)).*weg);
end

fr(ne1,1)=fr(ne1,1)+Nex'*(wxeg.*weg)...
                   +Ney'*(wyeg.*weg);
if nsd==3
  fr(ne1,1)=fr(ne1,1)+Nez'*(wzeg.*weg);
end

fr(ne1,1)=fr(ne1,1)+Ne'*(Rceg.*weg);

if isTimeDependent
  fw(ne1,1)=-Ne'*((1/dt*wxeg*alpha(1)...
                  +1/dt*woldxeg*alpha(2:BDFo+1,1)).*weg);
  fw(ne2,1)=-Ne'*((1/dt*wyeg*alpha(1)...
                  +1/dt*woldyeg*alpha(2:BDFo+1,1)).*weg);
  if nsd==3
    fw(ne3,1)=-Ne'*((1/dt*wzeg*alpha(1)...
                    +1/dt*woldzeg*alpha(2:BDFo+1,1)).*weg);
  end
end

if isConvectiveFlow
  if nsd==2
    fw(ne1,1)=fw(ne1,1)+Nex'*(wxeg.*wxeg./reg.*weg)...
                       +Ney'*(wxeg.*wyeg./reg.*weg);
    fw(ne2,1)=fw(ne2,1)+Nex'*(wyeg.*wxeg./reg.*weg)...
                       +Ney'*(wyeg.*wyeg./reg.*weg);
  elseif nsd==3
    fw(ne1,1)=fw(ne1,1)+Nex'*(wxeg.*wxeg./reg.*weg)...
                       +Ney'*(wxeg.*wyeg./reg.*weg)...
                       +Nez'*(wxeg.*wzeg./reg.*weg);
    fw(ne2,1)=fw(ne2,1)+Nex'*(wyeg.*wxeg./reg.*weg)...
                       +Ney'*(wyeg.*wyeg./reg.*weg)...
                       +Nez'*(wyeg.*wzeg./reg.*weg);
    fw(ne3,1)=fw(ne3,1)+Nex'*(wzeg.*wxeg./reg.*weg)...
                       +Ney'*(wzeg.*wyeg./reg.*weg)...
                       +Nez'*(wzeg.*wzeg./reg.*weg);
  end
end

if nsd==2
  fw(ne1,1)=fw(ne1,1)-Ne'*((+Voigt1*(Nex*Lxxe)+Voigt2*(Nex*Lyye)...
                            +Voigt3*(Ney*Lxye)).*weg);
  fw(ne2,1)=fw(ne2,1)-Ne'*((+Voigt2*(Ney*Lxxe)+Voigt1*(Ney*Lyye)...
                            +Voigt3*(Nex*Lxye)).*weg);
elseif nsd==3
  fw(ne1,1)=fw(ne1,1)-Ne'*((+Voigt1*(Nex*Lxxe)+Voigt2*(Nex*Lyye)+Voigt2*(Nex*Lzze)...
                            +Voigt3*(Ney*Lxye)+Voigt3*(Nez*Lxze)).*weg);
  fw(ne2,1)=fw(ne2,1)-Ne'*((+Voigt2*(Ney*Lxxe)+Voigt1*(Ney*Lyye)+Voigt2*(Ney*Lzze)...
                            +Voigt3*(Nex*Lxye)+Voigt3*(Nez*Lyze)).*weg);
  fw(ne3,1)=fw(ne3,1)-Ne'*((+Voigt2*(Nez*Lxxe)+Voigt2*(Nez*Lyye)+Voigt1*(Nez*Lzze)...
                            +Voigt3*(Nex*Lxze)+Voigt3*(Ney*Lyze)).*weg);
end

fw(ne1,1)=fw(ne1,1)+Nex'*(peg.*weg);
fw(ne2,1)=fw(ne2,1)+Ney'*(peg.*weg);
if nsd==3
  fw(ne3,1)=fw(ne3,1)+Nez'*(peg.*weg);
end

fw(ne1,1)=fw(ne1,1)+Ne'*(fxeg.*weg);
fw(ne2,1)=fw(ne2,1)+Ne'*(fyeg.*weg);
if nsd==3
  fw(ne3,1)=fw(ne3,1)+Ne'*(fzeg.*weg);
end

if isTimeDependent
  fe(ne1,1)=-Ne'*((1/dt*eeg*alpha(1)...
                  +1/dt*eoldeg*alpha(2:BDFo+1,1)).*weg);
end

fe(ne1,1)=fe(ne1,1)+Nex'*(eeg.*wxeg./reg.*weg)...
                   +Ney'*(eeg.*wyeg./reg.*weg);
if nsd==3
  fe(ne1,1)=fe(ne1,1)+Nez'*(eeg.*wzeg./reg.*weg);
end

if nsd==2
  fe(ne1,1)=fe(ne1,1)-Ne'*((+Nex*(1./re.*(+wxe.*(Voigt1*Lxxe+Voigt2*Lyye)...
                                          +wye.*(Voigt3*Lxye)))...
                            +Ney*(1./re.*(+wxe.*(Voigt3*Lxye)...
                                          +wye.*(Voigt2*Lxxe+Voigt1*Lyye)))...
                            +Nex*(sqrt(kappa)*Qxe)...
                            +Ney*(sqrt(kappa)*Qye)...
                            -(fxeg.*wxeg+fyeg.*wyeg)./reg).*weg);
elseif nsd==3
  fe(ne1,1)=fe(ne1,1)-Ne'*((+Nex*(1./re.*(+wxe.*(Voigt1*Lxxe+Voigt2*Lyye+Voigt2*Lzze)...
                                          +wye.*(Voigt3*Lxye)...
                                          +wze.*(Voigt3*Lxze)))...
                            +Ney*(1./re.*(+wxe.*(Voigt3*Lxye)...
                                          +wye.*(Voigt2*Lxxe+Voigt1*Lyye+Voigt2*Lzze)...
                                          +wze.*(Voigt3*Lyze)))...
                            +Nez*(1./re.*(+wxe.*(Voigt3*Lxze)...
                                          +wye.*(Voigt3*Lyze)...
                                          +wze.*(Voigt2*Lxxe+Voigt2*Lyye+Voigt1*Lzze)))...
                            +Nex*(sqrt(kappa)*Qxe)...
                            +Ney*(sqrt(kappa)*Qye)...
                            +Nez*(sqrt(kappa)*Qze)...
                            -(fxeg.*wxeg+fyeg.*wyeg+fzeg.*wzeg)./reg).*weg);
end

fe(ne1,1)=fe(ne1,1)+Nex'*(peg.*wxeg./reg.*weg)...
                   +Ney'*(peg.*wyeg./reg.*weg);  
if nsd==3
  fe(ne1,1)=fe(ne1,1)+Nez'*(peg.*wzeg./reg.*weg);
end

fe(ne1,1)=fe(ne1,1)+Ne'*(seg.*weg);

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
    NORMX=sparse(ngf,ngf,nx);
    NORMY=sparse(ngf,ngf,ny);
    if nsd==3
      NORMZ=sparse(ngf,ngf,nz);
    end
    WFG=sparse(ngf,ngf,wfg);
    
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
      pfg=p(Rfg,[Wxfg,Wyfg],Efg);
      dpdUfg=dpdu(Rfg,[Wxfg,Wyfg],Efg);
      dpdRfg=dpdUfg(:,1);
      dpdWxfg=dpdUfg(:,2);
      dpdWyfg=dpdUfg(:,3);
      dpdEfg=dpdUfg(:,4);
      Tfg=T(Rfg,[Wxfg,Wyfg],Efg);
      dTdUfg=dTdu(Rfg,[Wxfg,Wyfg],Efg);
      dTdRfg=dTdUfg(:,1);
      dTdWxfg=dTdUfg(:,2);
      dTdWyfg=dTdUfg(:,3);
      dTdEfg=dTdUfg(:,4);
    elseif nsd==3
      pfg=p(Rfg,[Wxfg,Wyfg,Wzfg],Efg);
      dpdUfg=dpdu(Rfg,[Wxfg,Wyfg,Wzfg],Efg);
      dpdRfg=dpdUfg(:,1);
      dpdWxfg=dpdUfg(:,2);
      dpdWyfg=dpdUfg(:,3);
      dpdWzfg=dpdUfg(:,4);
      dpdEfg=dpdUfg(:,5);
      Tfg=T(Rfg,[Wxfg,Wyfg,Wzfg],Efg);
      dTdUfg=dTdu(Rfg,[Wxfg,Wyfg,Wzfg],Efg);
      dTdRfg=dTdUfg(:,1);
      dTdWxfg=dTdUfg(:,2);
      dTdWyfg=dTdUfg(:,3);
      dTdWzfg=dTdUfg(:,4);
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
    Mf=Nf'*(WFG*Nf);
    Mf=(Mf+Mf')/2;
    Mfnx=Nf'*(NORMX*WFG*Nf);
    Mfny=Nf'*(NORMY*WFG*Nf);
    if nsd==3
      Mfnz=Nf'*(NORMZ*WFG*Nf);
    end
    
    % Compute common terms
    SSTRXX=sparse(ngf,ngf,Lxxfg);
    SSTRYY=sparse(ngf,ngf,Lyyfg);
    SSTRXY=sparse(ngf,ngf,Lxyfg);
    if nsd==3
      SSTRZZ=sparse(ngf,ngf,Lzzfg);
      SSTRXZ=sparse(ngf,ngf,Lxzfg);
      SSTRYZ=sparse(ngf,ngf,Lyzfg);
    end
    INVDENS=sparse(ngf,ngf,1./Rfg);
    MOMX=sparse(ngf,ngf,Wxfg);
    MOMY=sparse(ngf,ngf,Wyfg);
    if nsd==3
      MOMZ=sparse(ngf,ngf,Wzfg);
    end
    ENER=sparse(ngf,ngf,Efg);
    PRES=sparse(ngf,ngf,pfg);
    DTEMPDDENS=sparse(ngf,ngf,dTdRfg);
    DTEMPDMOMX=sparse(ngf,ngf,dTdWxfg);
    DTEMPDMOMY=sparse(ngf,ngf,dTdWyfg);
    if nsd==3
      DTEMPDMOMZ=sparse(ngf,ngf,dTdWzfg);
    end
    DTEMPDENER=sparse(ngf,ngf,dTdEfg);
    DPRESDDENS=sparse(ngf,ngf,dpdRfg);
    DPRESDMOMX=sparse(ngf,ngf,dpdWxfg);
    DPRESDMOMY=sparse(ngf,ngf,dpdWyfg);
    if nsd==3
      DPRESDMOMZ=sparse(ngf,ngf,dpdWzfg);
    end
    DPRESDENER=sparse(ngf,ngf,dpdEfg);
    Invdens_nx=INVDENS*NORMX*WFG;
    Invdens_ny=INVDENS*NORMY*WFG;
    if nsd==3
      Invdens_nz=INVDENS*NORMZ*WFG;
    end
    if nsd==2
      MomInvdens2_xnx=MOMX*(INVDENS.^2)*NORMX*WFG;
      MomInvdens2_xny=MOMX*(INVDENS.^2)*NORMY*WFG;
      MomInvdens2_ynx=MOMY*(INVDENS.^2)*NORMX*WFG;
      MomInvdens2_yny=MOMY*(INVDENS.^2)*NORMY*WFG;
    elseif nsd==3
      MomInvdens2_xnx=MOMX*(INVDENS.^2)*NORMX*WFG;
      MomInvdens2_xny=MOMX*(INVDENS.^2)*NORMY*WFG;
      MomInvdens2_xnz=MOMX*(INVDENS.^2)*NORMZ*WFG;
      MomInvdens2_ynx=MOMY*(INVDENS.^2)*NORMX*WFG;
      MomInvdens2_yny=MOMY*(INVDENS.^2)*NORMY*WFG;
      MomInvdens2_ynz=MOMY*(INVDENS.^2)*NORMZ*WFG;
      MomInvdens2_znx=MOMZ*(INVDENS.^2)*NORMX*WFG;
      MomInvdens2_zny=MOMZ*(INVDENS.^2)*NORMY*WFG;
      MomInvdens2_znz=MOMZ*(INVDENS.^2)*NORMZ*WFG;
    end
    dTempdDens_nx=DTEMPDDENS*NORMX*WFG;
    dTempdDens_ny=DTEMPDDENS*NORMY*WFG;
    if nsd==3
      dTempdDens_nz=DTEMPDDENS*NORMZ*WFG;
    end
    dTempdMomx_nx=DTEMPDMOMX*NORMX*WFG;
    dTempdMomx_ny=DTEMPDMOMX*NORMY*WFG;
    dTempdMomy_nx=DTEMPDMOMY*NORMX*WFG;
    dTempdMomy_ny=DTEMPDMOMY*NORMY*WFG;
    if nsd==3
      dTempdMomx_nz=DTEMPDMOMX*NORMZ*WFG;
      dTempdMomy_nz=DTEMPDMOMY*NORMZ*WFG;
      dTempdMomz_nx=DTEMPDMOMZ*NORMX*WFG;
      dTempdMomz_ny=DTEMPDMOMZ*NORMY*WFG;
      dTempdMomz_nz=DTEMPDMOMZ*NORMZ*WFG;
    end
    dTempdEner_nx=DTEMPDENER*NORMX*WFG;
    dTempdEner_ny=DTEMPDENER*NORMY*WFG;
    if nsd==3
      dTempdEner_nz=DTEMPDENER*NORMZ*WFG;
    end
    if nsd==2
      Conv_xnxyny_x=(MOMX*NORMX+MOMY*NORMY)*(INVDENS.^2)*MOMX*WFG;
      Conv_xnxyny_y=(MOMX*NORMX+MOMY*NORMY)*(INVDENS.^2)*MOMY*WFG;
    elseif nsd==3
      Conv_xnxynyznz_x=(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*(INVDENS.^2)*MOMX*WFG;
      Conv_xnxynyznz_y=(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*(INVDENS.^2)*MOMY*WFG;
      Conv_xnxynyznz_z=(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*(INVDENS.^2)*MOMZ*WFG;
    end
    if nsd==2
      Conv_xx1=(2*MOMX*INVDENS)*NORMX*WFG;
      Conv_xx2=(  MOMY*INVDENS)*NORMY*WFG;
      Conv_xy =(  MOMX*INVDENS)*NORMY*WFG;
      Conv_yx =(  MOMY*INVDENS)*NORMX*WFG;
      Conv_yy1=(  MOMX*INVDENS)*NORMX*WFG;
      Conv_yy2=(2*MOMY*INVDENS)*NORMY*WFG;
    elseif nsd==3
      Conv_xx1=(2*MOMX*INVDENS)*NORMX*WFG;
      Conv_xx2=(  MOMY*INVDENS)*NORMY*WFG;
      Conv_xx3=(  MOMZ*INVDENS)*NORMZ*WFG;
      Conv_xy =(  MOMX*INVDENS)*NORMY*WFG;
      Conv_xz =(  MOMX*INVDENS)*NORMZ*WFG;
      Conv_yx =(  MOMY*INVDENS)*NORMX*WFG;
      Conv_yy1=(  MOMX*INVDENS)*NORMX*WFG;
      Conv_yy2=(2*MOMY*INVDENS)*NORMY*WFG;
      Conv_yy3=(  MOMZ*INVDENS)*NORMZ*WFG;
      Conv_yz =(  MOMY*INVDENS)*NORMZ*WFG;
      Conv_zx =(  MOMZ*INVDENS)*NORMX*WFG;
      Conv_zy =(  MOMZ*INVDENS)*NORMY*WFG;
      Conv_zz1=(  MOMX*INVDENS)*NORMX*WFG;
      Conv_zz2=(  MOMY*INVDENS)*NORMY*WFG;
      Conv_zz3=(2*MOMZ*INVDENS)*NORMZ*WFG;
    end
    if nsd==2
      Invdens2MomnEner=(INVDENS.^2)*(MOMX*NORMX+MOMY*NORMY)*ENER*WFG;
    elseif nsd==3
      Invdens2MomnEner=(INVDENS.^2)*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*ENER*WFG;
    end
    InvdensEner_nx=INVDENS*ENER*NORMX*WFG;
    InvdensEner_ny=INVDENS*ENER*NORMY*WFG;
    if nsd==3
      InvdensEner_nz=INVDENS*ENER*NORMZ*WFG;
    end
    if nsd==2
      InvdensMomn=INVDENS*(MOMX*NORMX+MOMY*NORMY)*WFG;
    elseif nsd==3
      InvdensMomn=INVDENS*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
    end
    if nsd==2
      Stress_KELxx=INVDENS*(MOMX*Voigt1*NORMX+MOMY*Voigt2*NORMY)*WFG;
      Stress_KELyy=INVDENS*(MOMX*Voigt2*NORMX+MOMY*Voigt1*NORMY)*WFG;
      Stress_KELxy=INVDENS*(MOMY*Voigt3*NORMX+MOMX*Voigt3*NORMY)*WFG;
    elseif nsd==3
      Stress_KELxx=INVDENS*(MOMX*Voigt1*NORMX+MOMY*Voigt2*NORMY+MOMZ*Voigt2*NORMZ)*WFG;
      Stress_KELyy=INVDENS*(MOMX*Voigt2*NORMX+MOMY*Voigt1*NORMY+MOMZ*Voigt2*NORMZ)*WFG;
      Stress_KELzz=INVDENS*(MOMX*Voigt2*NORMX+MOMY*Voigt2*NORMY+MOMZ*Voigt1*NORMZ)*WFG;
      Stress_KELxy=INVDENS*(MOMY*Voigt3*NORMX+MOMX*Voigt3*NORMY)*WFG;
      Stress_KELxz=INVDENS*(MOMZ*Voigt3*NORMX+MOMX*Voigt3*NORMZ)*WFG;
      Stress_KELyz=INVDENS*(MOMZ*Voigt3*NORMY+MOMY*Voigt3*NORMZ)*WFG;
    end
    if nsd==2
      Stress_KER=(INVDENS.^2)*((MOMX*(Voigt1*SSTRXX+Voigt2*SSTRYY)...
                               +MOMY*Voigt3*SSTRXY)*NORMX...
                              +(MOMY*(Voigt2*SSTRXX+Voigt1*SSTRYY)...
                               +MOMX*Voigt3*SSTRXY)*NORMY)*WFG;
    elseif nsd==3
      Stress_KER=(INVDENS.^2)*((MOMX*(Voigt1*SSTRXX+Voigt2*SSTRYY+Voigt2*SSTRZZ)...
                               +MOMY*Voigt3*SSTRXY+MOMZ*Voigt3*SSTRXZ)*NORMX...
                              +(MOMY*(Voigt2*SSTRXX+Voigt1*SSTRYY+Voigt2*SSTRZZ)...
                               +MOMX*Voigt3*SSTRXY+MOMZ*Voigt3*SSTRYZ)*NORMY...
                              +(MOMZ*(Voigt2*SSTRXX+Voigt2*SSTRYY+Voigt1*SSTRZZ)...
                               +MOMX*Voigt3*SSTRXZ+MOMY*Voigt3*SSTRYZ)*NORMZ)*WFG;
    end
    if nsd==2
      Stress_KEWx=INVDENS*((Voigt1*SSTRXX+Voigt2*SSTRYY)*NORMX...
                                        +(Voigt3*SSTRXY)*NORMY)*WFG;
      Stress_KEWy=INVDENS*((Voigt2*SSTRXX+Voigt1*SSTRYY)*NORMY...
                                        +(Voigt3*SSTRXY)*NORMX)*WFG;
    elseif nsd==3
      Stress_KEWx=INVDENS*((Voigt1*SSTRXX+Voigt2*SSTRYY+Voigt2*SSTRZZ)*NORMX...
                                                      +(Voigt3*SSTRXY)*NORMY...
                                                      +(Voigt3*SSTRXZ)*NORMZ)*WFG;
      Stress_KEWy=INVDENS*((Voigt2*SSTRXX+Voigt1*SSTRYY+Voigt2*SSTRZZ)*NORMY...
                                                      +(Voigt3*SSTRXY)*NORMX...
                                                      +(Voigt3*SSTRYZ)*NORMZ)*WFG;
      Stress_KEWz=INVDENS*((Voigt2*SSTRXX+Voigt2*SSTRYY+Voigt1*SSTRZZ)*NORMZ...
                                                      +(Voigt3*SSTRXZ)*NORMX...
                                                      +(Voigt3*SSTRYZ)*NORMY)*WFG;
    end
    dPresdDens_nx=DPRESDDENS*NORMX*WFG;
    dPresdDens_ny=DPRESDDENS*NORMY*WFG;
    if nsd==3
      dPresdDens_nz=DPRESDDENS*NORMZ*WFG;
    end
    dPresdMomx_nx=DPRESDMOMX*NORMX*WFG;
    dPresdMomx_ny=DPRESDMOMX*NORMY*WFG;
    dPresdMomy_nx=DPRESDMOMY*NORMX*WFG;
    dPresdMomy_ny=DPRESDMOMY*NORMY*WFG;
    if nsd==3
      dPresdMomx_nz=DPRESDMOMX*NORMZ*WFG;
      dPresdMomy_nz=DPRESDMOMY*NORMZ*WFG;
      dPresdMomz_nx=DPRESDMOMZ*NORMX*WFG;
      dPresdMomz_ny=DPRESDMOMZ*NORMY*WFG;
      dPresdMomz_nz=DPRESDMOMZ*NORMZ*WFG;
    end
    dPresdEner_nx=DPRESDENER*NORMX*WFG;
    dPresdEner_ny=DPRESDENER*NORMY*WFG;
    if nsd==3
      dPresdEner_nz=DPRESDENER*NORMZ*WFG;
    end
    if nsd==2
      dPresdDensVeln=DPRESDDENS*INVDENS*(MOMX*NORMX+MOMY*NORMY)*WFG;
    elseif nsd==3
      dPresdDensVeln=DPRESDDENS*INVDENS*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
    end
    if nsd==2
      PresInvdens2Momn=PRES*(INVDENS.^2)*(MOMX*NORMX+MOMY*NORMY)*WFG;
    elseif nsd==3
      PresInvdens2Momn=PRES*(INVDENS.^2)*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
    end
    if nsd==2
      dPresdMomxVeln=DPRESDMOMX*INVDENS*(MOMX*NORMX+MOMY*NORMY)*WFG;
      dPresdMomyVeln=DPRESDMOMY*INVDENS*(MOMX*NORMX+MOMY*NORMY)*WFG;
    elseif nsd==3
      dPresdMomxVeln=DPRESDMOMX*INVDENS*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
      dPresdMomyVeln=DPRESDMOMY*INVDENS*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
      dPresdMomzVeln=DPRESDMOMZ*INVDENS*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
    end
    PresInvdens_nx=PRES*INVDENS*NORMX*WFG;
    PresInvdens_ny=PRES*INVDENS*NORMY*WFG;
    if nsd==3
      PresInvdens_nz=PRES*INVDENS*NORMZ*WFG;
    end
    if nsd==2
      dPresdEnerVeln=DPRESDENER*INVDENS*(MOMX*NORMX+MOMY*NORMY)*WFG;
    elseif nsd==3
      dPresdEnerVeln=DPRESDENER*INVDENS*(MOMX*NORMX+MOMY*NORMY+MOMZ*NORMZ)*WFG;
    end
    
    % Compute lhs
    if nsd==2
      KLR(nf1,nefR1)=KLR(nf1,nefR1)+Voigt1*Nf'*(MomInvdens2_xnx*Nf)...
                                   +Voigt2*Nf'*(MomInvdens2_yny*Nf);
      KLR(nf2,nefR1)=KLR(nf2,nefR1)+Voigt2*Nf'*(MomInvdens2_xnx*Nf)...
                                   +Voigt1*Nf'*(MomInvdens2_yny*Nf);
      KLR(nf3,nefR1)=KLR(nf3,nefR1)+Voigt3*Nf'*(MomInvdens2_xny*Nf)...
                                   +Voigt3*Nf'*(MomInvdens2_ynx*Nf);
    elseif nsd==3
      KLR(nf1,nefR1)=KLR(nf1,nefR1)+Voigt1*Nf'*(MomInvdens2_xnx*Nf)...
                                   +Voigt2*Nf'*(MomInvdens2_yny*Nf)...
                                   +Voigt2*Nf'*(MomInvdens2_znz*Nf);
      KLR(nf2,nefR1)=KLR(nf2,nefR1)+Voigt2*Nf'*(MomInvdens2_xnx*Nf)...
                                   +Voigt1*Nf'*(MomInvdens2_yny*Nf)...
                                   +Voigt2*Nf'*(MomInvdens2_znz*Nf);
      KLR(nf3,nefR1)=KLR(nf3,nefR1)+Voigt2*Nf'*(MomInvdens2_xnx*Nf)...
                                   +Voigt2*Nf'*(MomInvdens2_yny*Nf)...
                                   +Voigt1*Nf'*(MomInvdens2_znz*Nf);
      KLR(nf4,nefR1)=KLR(nf4,nefR1)+Voigt3*Nf'*(MomInvdens2_xny*Nf)...
                                   +Voigt3*Nf'*(MomInvdens2_ynx*Nf);
      KLR(nf5,nefR1)=KLR(nf5,nefR1)+Voigt3*Nf'*(MomInvdens2_xnz*Nf)...
                                   +Voigt3*Nf'*(MomInvdens2_znx*Nf);
      KLR(nf6,nefR1)=KLR(nf6,nefR1)+Voigt3*Nf'*(MomInvdens2_ynz*Nf)...
                                   +Voigt3*Nf'*(MomInvdens2_zny*Nf);
    end
      
    if nsd==2
      KLW(nf1,nefW1)=KLW(nf1,nefW1)-Voigt1*Nf'*(Invdens_nx*Nf);
      KLW(nf2,nefW1)=KLW(nf2,nefW1)-Voigt2*Nf'*(Invdens_nx*Nf);
      KLW(nf3,nefW1)=KLW(nf3,nefW1)-Voigt3*Nf'*(Invdens_ny*Nf);
      KLW(nf1,nefW2)=KLW(nf1,nefW2)-Voigt2*Nf'*(Invdens_ny*Nf);
      KLW(nf2,nefW2)=KLW(nf2,nefW2)-Voigt1*Nf'*(Invdens_ny*Nf);
      KLW(nf3,nefW2)=KLW(nf3,nefW2)-Voigt3*Nf'*(Invdens_nx*Nf);
    elseif nsd==3
      KLW(nf1,nefW1)=KLW(nf1,nefW1)-Voigt1*Nf'*(Invdens_nx*Nf);
      KLW(nf2,nefW1)=KLW(nf2,nefW1)-Voigt2*Nf'*(Invdens_nx*Nf);
      KLW(nf3,nefW1)=KLW(nf3,nefW1)-Voigt2*Nf'*(Invdens_nx*Nf);
      KLW(nf4,nefW1)=KLW(nf4,nefW1)-Voigt3*Nf'*(Invdens_ny*Nf);
      KLW(nf5,nefW1)=KLW(nf5,nefW1)-Voigt3*Nf'*(Invdens_nz*Nf);
      KLW(nf1,nefW2)=KLW(nf1,nefW2)-Voigt2*Nf'*(Invdens_ny*Nf);
      KLW(nf2,nefW2)=KLW(nf2,nefW2)-Voigt1*Nf'*(Invdens_ny*Nf);
      KLW(nf3,nefW2)=KLW(nf3,nefW2)-Voigt2*Nf'*(Invdens_ny*Nf);
      KLW(nf4,nefW2)=KLW(nf4,nefW2)-Voigt3*Nf'*(Invdens_nx*Nf);
      KLW(nf6,nefW2)=KLW(nf6,nefW2)-Voigt3*Nf'*(Invdens_nz*Nf);
      KLW(nf1,nefW3)=KLW(nf1,nefW3)-Voigt2*Nf'*(Invdens_nz*Nf);
      KLW(nf2,nefW3)=KLW(nf2,nefW3)-Voigt2*Nf'*(Invdens_nz*Nf);
      KLW(nf3,nefW3)=KLW(nf3,nefW3)-Voigt1*Nf'*(Invdens_nz*Nf);
      KLW(nf5,nefW3)=KLW(nf5,nefW3)-Voigt3*Nf'*(Invdens_nx*Nf);
      KLW(nf6,nefW3)=KLW(nf6,nefW3)-Voigt3*Nf'*(Invdens_ny*Nf);
    end
    
    KQR(nf1,nefR1)=KQR(nf1,nefR1)-sqrt(kappa)*Nf'*(dTempdDens_nx*Nf);
    KQR(nf2,nefR1)=KQR(nf2,nefR1)-sqrt(kappa)*Nf'*(dTempdDens_ny*Nf);
    if nsd==3
      KQR(nf3,nefR1)=KQR(nf3,nefR1)-sqrt(kappa)*Nf'*(dTempdDens_nz*Nf);
    end
    
    if nsd==2
      KQW(nf1,nefW1)=KQW(nf1,nefW1)-sqrt(kappa)*Nf'*(dTempdMomx_nx*Nf);
      KQW(nf2,nefW1)=KQW(nf2,nefW1)-sqrt(kappa)*Nf'*(dTempdMomx_ny*Nf);
      KQW(nf1,nefW2)=KQW(nf1,nefW2)-sqrt(kappa)*Nf'*(dTempdMomy_nx*Nf);
      KQW(nf2,nefW2)=KQW(nf2,nefW2)-sqrt(kappa)*Nf'*(dTempdMomy_ny*Nf);
    elseif nsd==3
      KQW(nf1,nefW1)=KQW(nf1,nefW1)-sqrt(kappa)*Nf'*(dTempdMomx_nx*Nf);
      KQW(nf2,nefW1)=KQW(nf2,nefW1)-sqrt(kappa)*Nf'*(dTempdMomx_ny*Nf);
      KQW(nf3,nefW1)=KQW(nf3,nefW1)-sqrt(kappa)*Nf'*(dTempdMomx_nz*Nf);
      KQW(nf1,nefW2)=KQW(nf1,nefW2)-sqrt(kappa)*Nf'*(dTempdMomy_nx*Nf);
      KQW(nf2,nefW2)=KQW(nf2,nefW2)-sqrt(kappa)*Nf'*(dTempdMomy_ny*Nf);
      KQW(nf3,nefW2)=KQW(nf3,nefW2)-sqrt(kappa)*Nf'*(dTempdMomy_nz*Nf);
      KQW(nf1,nefW3)=KQW(nf1,nefW3)-sqrt(kappa)*Nf'*(dTempdMomz_nx*Nf);
      KQW(nf2,nefW3)=KQW(nf2,nefW3)-sqrt(kappa)*Nf'*(dTempdMomz_ny*Nf);
      KQW(nf3,nefW3)=KQW(nf3,nefW3)-sqrt(kappa)*Nf'*(dTempdMomz_nz*Nf);
    end
    
    KQE(nf1,nefE1)=KQE(nf1,nefE1)-sqrt(kappa)*Nf'*(dTempdEner_nx*Nf);
    KQE(nf2,nefE1)=KQE(nf2,nefE1)-sqrt(kappa)*Nf'*(dTempdEner_ny*Nf);
    if nsd==3
      KQE(nf3,nefE1)=KQE(nf3,nefE1)-sqrt(kappa)*Nf'*(dTempdEner_nz*Nf);
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
        KwR(nf1,nefR1)=KwR(nf1,nefR1)-Nf'*(Conv_xnxyny_x*Nf);
        KwR(nf2,nefR1)=KwR(nf2,nefR1)-Nf'*(Conv_xnxyny_y*Nf);
      elseif nsd==3
        KwR(nf1,nefR1)=KwR(nf1,nefR1)-Nf'*(Conv_xnxynyznz_x*Nf);
        KwR(nf2,nefR1)=KwR(nf2,nefR1)-Nf'*(Conv_xnxynyznz_y*Nf);
        KwR(nf3,nefR1)=KwR(nf3,nefR1)-Nf'*(Conv_xnxynyznz_z*Nf);
      end
    end
    
    if isConvectiveFlow
      if nsd==2
        KwW(nf1,nefW1)=KwW(nf1,nefW1)+Nf'*((Conv_xx1+Conv_xx2)*Nf);
        KwW(nf2,nefW1)=KwW(nf2,nefW1)+Nf'* (Conv_yx*Nf);
        KwW(nf1,nefW2)=KwW(nf1,nefW2)+Nf'* (Conv_xy*Nf);
        KwW(nf2,nefW2)=KwW(nf2,nefW2)+Nf'*((Conv_yy1+Conv_yy2)*Nf);
      elseif nsd==3
        KwW(nf1,nefW1)=KwW(nf1,nefW1)+Nf'*((Conv_xx1+Conv_xx2+Conv_xx3)*Nf);
        KwW(nf2,nefW1)=KwW(nf2,nefW1)+Nf'* (Conv_yx*Nf);
        KwW(nf3,nefW1)=KwW(nf3,nefW1)+Nf'* (Conv_zx*Nf);
        KwW(nf1,nefW2)=KwW(nf1,nefW2)+Nf'* (Conv_xy*Nf);
        KwW(nf2,nefW2)=KwW(nf2,nefW2)+Nf'*((Conv_yy1+Conv_yy2+Conv_yy3)*Nf);
        KwW(nf3,nefW2)=KwW(nf3,nefW2)+Nf'* (Conv_zy*Nf);
        KwW(nf1,nefW3)=KwW(nf1,nefW3)+Nf'* (Conv_xz*Nf);
        KwW(nf2,nefW3)=KwW(nf2,nefW3)+Nf'* (Conv_yz*Nf);
        KwW(nf3,nefW3)=KwW(nf3,nefW3)+Nf'*((Conv_zz1+Conv_zz2+Conv_zz3)*Nf);
      end
    end
    
    KwR(nf1,nefR1)=KwR(nf1,nefR1)+Nf'*(dPresdDens_nx*Nf);
    KwR(nf2,nefR1)=KwR(nf2,nefR1)+Nf'*(dPresdDens_ny*Nf);
    if nsd==3
      KwR(nf3,nefR1)=KwR(nf3,nefR1)+Nf'*(dPresdDens_nz*Nf);
    end
    
    if nsd==2
      KwW(nf1,nefW1)=KwW(nf1,nefW1)+Nf'*(dPresdMomx_nx*Nf);
      KwW(nf2,nefW1)=KwW(nf2,nefW1)+Nf'*(dPresdMomx_ny*Nf);
      KwW(nf1,nefW2)=KwW(nf1,nefW2)+Nf'*(dPresdMomy_nx*Nf);
      KwW(nf2,nefW2)=KwW(nf2,nefW2)+Nf'*(dPresdMomy_ny*Nf);
    elseif nsd==3
      KwW(nf1,nefW1)=KwW(nf1,nefW1)+Nf'*(dPresdMomx_nx*Nf);
      KwW(nf2,nefW1)=KwW(nf2,nefW1)+Nf'*(dPresdMomx_ny*Nf);
      KwW(nf3,nefW1)=KwW(nf3,nefW1)+Nf'*(dPresdMomx_nz*Nf);
      KwW(nf1,nefW2)=KwW(nf1,nefW2)+Nf'*(dPresdMomy_nx*Nf);
      KwW(nf2,nefW2)=KwW(nf2,nefW2)+Nf'*(dPresdMomy_ny*Nf);
      KwW(nf3,nefW2)=KwW(nf3,nefW2)+Nf'*(dPresdMomy_nz*Nf);
      KwW(nf1,nefW3)=KwW(nf1,nefW3)+Nf'*(dPresdMomz_nx*Nf);
      KwW(nf2,nefW3)=KwW(nf2,nefW3)+Nf'*(dPresdMomz_ny*Nf);
      KwW(nf3,nefW3)=KwW(nf3,nefW3)+Nf'*(dPresdMomz_nz*Nf);
    end
    
    KwE(nf1,nefE1)=KwE(nf1,nefE1)+Nf'*(dPresdEner_nx*Nf);
    KwE(nf2,nefE1)=KwE(nf2,nefE1)+Nf'*(dPresdEner_ny*Nf);
    if nsd==3
      KwE(nf3,nefE1)=KwE(nf3,nefE1)+Nf'*(dPresdEner_nz*Nf);
    end
    
    KwW(nf1,nefW1)=KwW(nf1,nefW1)-tauW*Mf;
    KwW(nf2,nefW2)=KwW(nf2,nefW2)-tauW*Mf;
    if nsd==3
      KwW(nf3,nefW3)=KwW(nf3,nefW3)-tauW*Mf;
    end
    
    Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*Mf;
    
    KeR(nf1,nefR1)=KeR(nf1,nefR1)-Nf'*(Invdens2MomnEner*Nf);
    
    KeW(nf1,nefW1)=KeW(nf1,nefW1)+Nf'*(InvdensEner_nx*Nf);
    KeW(nf1,nefW2)=KeW(nf1,nefW2)+Nf'*(InvdensEner_ny*Nf);
    if nsd==3
      KeW(nf1,nefW3)=KeW(nf1,nefW3)+Nf'*(InvdensEner_nz*Nf);
    end
    
    KeE(nf1,nefE1)=KeE(nf1,nefE1)+Nf'*(InvdensMomn*Nf);
    
    KeR(nf1,nefR1)=KeR(nf1,nefR1)+Nf'*(dPresdDensVeln*Nf);
    
    KeR(nf1,nefR1)=KeR(nf1,nefR1)-Nf'*(PresInvdens2Momn*Nf);
    
    KeW(nf1,nefW1)=KeW(nf1,nefW1)+Nf'*(dPresdMomxVeln*Nf);
    KeW(nf1,nefW2)=KeW(nf1,nefW2)+Nf'*(dPresdMomyVeln*Nf);
    if nsd==3
      KeW(nf1,nefW3)=KeW(nf1,nefW3)+Nf'*(dPresdMomzVeln*Nf);
    end
    
    KeW(nf1,nefW1)=KeW(nf1,nefW1)+Nf'*(PresInvdens_nx*Nf);
    KeW(nf1,nefW2)=KeW(nf1,nefW2)+Nf'*(PresInvdens_ny*Nf);
    if nsd==3
      KeW(nf1,nefW3)=KeW(nf1,nefW3)+Nf'*(PresInvdens_nz*Nf);
    end
    
    KeE(nf1,nefE1)=KeE(nf1,nefE1)+Nf'*(dPresdEnerVeln*Nf);
    
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
        KWR(nefW1,nefR1)=KWR(nefW1,nefR1)+Nf'*(Conv_xnxyny_x*Nf);
        KWR(nefW2,nefR1)=KWR(nefW2,nefR1)+Nf'*(Conv_xnxyny_y*Nf);
      elseif nsd==3
        KWR(nefW1,nefR1)=KWR(nefW1,nefR1)+Nf'*(Conv_xnxynyznz_x*Nf);
        KWR(nefW2,nefR1)=KWR(nefW2,nefR1)+Nf'*(Conv_xnxynyznz_y*Nf);
        KWR(nefW3,nefR1)=KWR(nefW3,nefR1)+Nf'*(Conv_xnxynyznz_z*Nf);
      end
    end
    
    if not(isDirichlet) && isConvectiveFlow
      if nsd==2
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-Nf'*((Conv_xx1+Conv_xx2)*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-Nf'* (Conv_yx*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-Nf'* (Conv_xy*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-Nf'*((Conv_yy1+Conv_yy2)*Nf);
      elseif nsd==3
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-Nf'*((Conv_xx1+Conv_xx2+Conv_xx3)*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-Nf'* (Conv_yx*Nf);
        KWW(nefW3,nefW1)=KWW(nefW3,nefW1)-Nf'* (Conv_zx*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-Nf'* (Conv_xy*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-Nf'*((Conv_yy1+Conv_yy2+Conv_yy3)*Nf);
        KWW(nefW3,nefW2)=KWW(nefW3,nefW2)-Nf'* (Conv_zy*Nf);
        KWW(nefW1,nefW3)=KWW(nefW1,nefW3)-Nf'* (Conv_xz*Nf);
        KWW(nefW2,nefW3)=KWW(nefW2,nefW3)-Nf'* (Conv_yz*Nf);
        KWW(nefW3,nefW3)=KWW(nefW3,nefW3)-Nf'*((Conv_zz1+Conv_zz2+Conv_zz3)*Nf);
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
      KWR(nefW1,nefR1)=KWR(nefW1,nefR1)-Nf'*(dPresdDens_nx*Nf);
      KWR(nefW2,nefR1)=KWR(nefW2,nefR1)-Nf'*(dPresdDens_ny*Nf);
      if nsd==3
        KWR(nefW3,nefR1)=KWR(nefW3,nefR1)-Nf'*(dPresdDens_nz*Nf);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-Nf'*(dPresdMomx_nx*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-Nf'*(dPresdMomx_ny*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-Nf'*(dPresdMomy_nx*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-Nf'*(dPresdMomy_ny*Nf);
      elseif nsd==3
        KWW(nefW1,nefW1)=KWW(nefW1,nefW1)-Nf'*(dPresdMomx_nx*Nf);
        KWW(nefW2,nefW1)=KWW(nefW2,nefW1)-Nf'*(dPresdMomx_ny*Nf);
        KWW(nefW3,nefW1)=KWW(nefW3,nefW1)-Nf'*(dPresdMomx_nz*Nf);
        KWW(nefW1,nefW2)=KWW(nefW1,nefW2)-Nf'*(dPresdMomy_nx*Nf);
        KWW(nefW2,nefW2)=KWW(nefW2,nefW2)-Nf'*(dPresdMomy_ny*Nf);
        KWW(nefW3,nefW2)=KWW(nefW3,nefW2)-Nf'*(dPresdMomy_nz*Nf);
        KWW(nefW1,nefW3)=KWW(nefW1,nefW3)-Nf'*(dPresdMomz_nx*Nf);
        KWW(nefW2,nefW3)=KWW(nefW2,nefW3)-Nf'*(dPresdMomz_ny*Nf);
        KWW(nefW3,nefW3)=KWW(nefW3,nefW3)-Nf'*(dPresdMomz_nz*Nf);
      end
    end
    
    if not(isDirichlet)
      KWE(nefW1,nefE1)=KWE(nefW1,nefE1)-Nf'*(dPresdEner_nx*Nf);
      KWE(nefW2,nefE1)=KWE(nefW2,nefE1)-Nf'*(dPresdEner_ny*Nf);
      if nsd==3
        KWE(nefW3,nefE1)=KWE(nefW3,nefE1)-Nf'*(dPresdEner_nz*Nf);
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
      KER(nefE1,nefR1)=KER(nefE1,nefR1)+Nf'*(Invdens2MomnEner*Nf);
    end
    
    if not(isDirichlet)
      KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-Nf'*(InvdensEner_nx*Nf);
      KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-Nf'*(InvdensEner_ny*Nf);
      if nsd==3
        KEW(nefE1,nefW3)=KEW(nefE1,nefW3)-Nf'*(InvdensEner_nz*Nf);
      end
    end
    
    if not(isDirichlet)
      KEE(nefE1,nefE1)=KEE(nefE1,nefE1)-Nf'*(InvdensMomn*Nf);
    end
    
    if not(isDirichlet)
      if nsd==2
        KEL(nefE1,nf1)=KEL(nefE1,nf1)-Nf'*(Stress_KELxx*Nf);
        KEL(nefE1,nf2)=KEL(nefE1,nf2)-Nf'*(Stress_KELyy*Nf);
        KEL(nefE1,nf3)=KEL(nefE1,nf3)-Nf'*(Stress_KELxy*Nf);
      elseif nsd==3
        KEL(nefE1,nf1)=KEL(nefE1,nf1)-Nf'*(Stress_KELxx*Nf);
        KEL(nefE1,nf2)=KEL(nefE1,nf2)-Nf'*(Stress_KELyy*Nf);
        KEL(nefE1,nf3)=KEL(nefE1,nf3)-Nf'*(Stress_KELzz*Nf);
        KEL(nefE1,nf4)=KEL(nefE1,nf4)-Nf'*(Stress_KELxy*Nf);
        KEL(nefE1,nf5)=KEL(nefE1,nf5)-Nf'*(Stress_KELxz*Nf);
        KEL(nefE1,nf6)=KEL(nefE1,nf6)-Nf'*(Stress_KELyz*Nf);
      end
    end
    
    if not(isDirichlet)
      KER(nefE1,nefR1)=KER(nefE1,nefR1)+Nf'*(Stress_KER*Nf);
    end
    
    if not(isDirichlet)
      KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-Nf'*(Stress_KEWx*Nf);
      KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-Nf'*(Stress_KEWy*Nf);
      if nsd==3
        KEW(nefE1,nefW3)=KEW(nefE1,nefW3)-Nf'*(Stress_KEWz*Nf);
      end
    end
    
    if not(isDirichlet)
      KER(nefE1,nefR1)=KER(nefE1,nefR1)-Nf'*(dPresdDensVeln*Nf);
    end
    
    if not(isDirichlet)
      KER(nefE1,nefR1)=KER(nefE1,nefR1)+Nf'*(PresInvdens2Momn*Nf);
    end
    
    if not(isDirichlet)
      KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-Nf'*(dPresdMomxVeln*Nf);
      KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-Nf'*(dPresdMomyVeln*Nf);
      if nsd==3
        KEW(nefE1,nefW3)=KEW(nefE1,nefW3)-Nf'*(dPresdMomzVeln*Nf);
      end
    end
    
    if not(isDirichlet)
      KEW(nefE1,nefW1)=KEW(nefE1,nefW1)-Nf'*(PresInvdens_nx*Nf);
      KEW(nefE1,nefW2)=KEW(nefE1,nefW2)-Nf'*(PresInvdens_ny*Nf);
      if nsd==3
        KEW(nefE1,nefW3)=KEW(nefE1,nefW3)-Nf'*(PresInvdens_nz*Nf);
      end
    end
    
    if not(isDirichlet)
      KEE(nefE1,nefE1)=KEE(nefE1,nefE1)-Nf'*(dPresdEnerVeln*Nf);
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
      fL(nf1,1)=fL(nf1,1)+Nf'*((+Voigt1*nx.*Wxfg./Rfg...
                                +Voigt2*ny.*Wyfg./Rfg).*wfg);
      fL(nf2,1)=fL(nf2,1)+Nf'*((+Voigt2*nx.*Wxfg./Rfg...
                                +Voigt1*ny.*Wyfg./Rfg).*wfg);
      fL(nf3,1)=fL(nf3,1)+Nf'*((+Voigt3*ny.*Wxfg./Rfg...
                                +Voigt3*nx.*Wyfg./Rfg).*wfg);
    elseif nsd==3
      fL(nf1,1)=fL(nf1,1)+Nf'*((+Voigt1*nx.*Wxfg./Rfg...
                                +Voigt2*ny.*Wyfg./Rfg...
                                +Voigt2*nz.*Wzfg./Rfg).*wfg);
      fL(nf2,1)=fL(nf2,1)+Nf'*((+Voigt2*nx.*Wxfg./Rfg...
                                +Voigt1*ny.*Wyfg./Rfg...
                                +Voigt2*nz.*Wzfg./Rfg).*wfg);
      fL(nf3,1)=fL(nf3,1)+Nf'*((+Voigt2*nx.*Wxfg./Rfg...
                                +Voigt2*ny.*Wyfg./Rfg...
                                +Voigt1*nz.*Wzfg./Rfg).*wfg);
      fL(nf4,1)=fL(nf4,1)+Nf'*((+Voigt3*nx.*Wyfg./Rfg...
                                +Voigt3*ny.*Wxfg./Rfg).*wfg);
      fL(nf5,1)=fL(nf5,1)+Nf'*((+Voigt3*nx.*Wzfg./Rfg...
                                +Voigt3*nz.*Wxfg./Rfg).*wfg);
      fL(nf6,1)=fL(nf6,1)+Nf'*((+Voigt3*ny.*Wzfg./Rfg...
                                +Voigt3*nz.*Wyfg./Rfg).*wfg);
    end
    
    fQ(nf1,1)=fQ(nf1,1)+Nf'*(sqrt(kappa)*nx.*Tfg.*wfg);
    fQ(nf2,1)=fQ(nf2,1)+Nf'*(sqrt(kappa)*ny.*Tfg.*wfg);
    if nsd==3
      fQ(nf3,1)=fQ(nf3,1)+Nf'*(sqrt(kappa)*nz.*Tfg.*wfg);
    end
    
    fr(nf1,1)=fr(nf1,1)-Nf'*(tauR*rfg.*wfg);
    
    fr(nf1,1)=fr(nf1,1)-Nf'*((Wxfg.*nx+Wyfg.*ny).*wfg);
    if nsd==3
      fr(nf1,1)=fr(nf1,1)-Nf'*(Wzfg.*nz.*wfg);
    end
    
    fr(nf1,1)=fr(nf1,1)+Nf'*(tauR*Rfg.*wfg);
    
    fw(nf1,1)=fw(nf1,1)-Nf'*(tauW*wxfg.*wfg);
    fw(nf2,1)=fw(nf2,1)-Nf'*(tauW*wyfg.*wfg);
    if nsd==3
      fw(nf3,1)=fw(nf3,1)-Nf'*(tauW*wzfg.*wfg);
    end
    
    if isConvectiveFlow
      if nsd==2
        fw(nf1,1)=fw(nf1,1)-Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wxfg.*wfg);
        fw(nf2,1)=fw(nf2,1)-Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wyfg.*wfg);
      elseif nsd==3
        fw(nf1,1)=fw(nf1,1)-Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wxfg.*wfg);
        fw(nf2,1)=fw(nf2,1)-Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wyfg.*wfg);
        fw(nf3,1)=fw(nf3,1)-Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wzfg.*wfg);
      end
    end
    
    fw(nf1,1)=fw(nf1,1)-Nf'*(pfg.*nx.*wfg);
    fw(nf2,1)=fw(nf2,1)-Nf'*(pfg.*ny.*wfg);
    if nsd==3
      fw(nf3,1)=fw(nf3,1)-Nf'*(pfg.*nz.*wfg);
    end
    
    fw(nf1,1)=fw(nf1,1)+Nf'*(tauW*Wxfg.*wfg);
    fw(nf2,1)=fw(nf2,1)+Nf'*(tauW*Wyfg.*wfg);
    if nsd==3
      fw(nf3,1)=fw(nf3,1)+Nf'*(tauW*Wzfg.*wfg);
    end
    
    fe(nf1,1)=fe(nf1,1)-Nf'*(tauE*efg.*wfg);
    
    fe(nf1,1)=fe(nf1,1)-Nf'*(Efg./Rfg.*(Wxfg.*nx+Wyfg.*ny).*wfg);
    if nsd==3
      fe(nf1,1)=fe(nf1,1)-Nf'*(Efg./Rfg.*(Wzfg.*nz).*wfg);
    end
    
    fe(nf1,1)=fe(nf1,1)-Nf'*(pfg./Rfg.*(Wxfg.*nx+Wyfg.*ny).*wfg);
    if nsd==3
      fe(nf1,1)=fe(nf1,1)-Nf'*(pfg./Rfg.*Wzfg.*nz.*wfg);
    end
    
    fe(nf1,1)=fe(nf1,1)+Nf'*(tauE*Efg.*wfg);
    
    if not(isDirichlet)
      fR(nefR1,1)=fR(nefR1,1)+Nf'*((Wxfg.*nx+Wyfg.*ny).*wfg);
      if nsd==3
        fR(nefR1,1)=fR(nefR1,1)+Nf'*(Wzfg.*nz.*wfg);
      end
    end
    
    if not(isDirichlet)
      fR(nefR1,1)=fR(nefR1,1)+Nf'*(tauR*rfg.*wfg);
    end
    
    fR(nefR1,1)=fR(nefR1,1)-Nf'*(tauR*Rfg.*wfg);
    
    if isDirichlet
      fR(nefR1,1)=fR(nefR1,1)+Nf'*(tauR*rDfg.*wfg);
    end
    
    if not(isDirichlet) && isConvectiveFlow
      if nsd==2
        fW(nefW1,1)=fW(nefW1,1)+Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wxfg.*wfg);
        fW(nefW2,1)=fW(nefW2,1)+Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny).*Wyfg.*wfg);
      elseif nsd==3
        fW(nefW1,1)=fW(nefW1,1)+Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wxfg.*wfg);
        fW(nefW2,1)=fW(nefW2,1)+Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wyfg.*wfg);
        fW(nefW3,1)=fW(nefW3,1)+Nf'*((Wxfg./Rfg.*nx+Wyfg./Rfg.*ny+Wzfg./Rfg.*nz).*Wzfg.*wfg);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        fW(nefW1,1)=fW(nefW1,1)+Nf'*((+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg...
                                      +Voigt3*ny.*Lxyfg).*wfg);
        fW(nefW2,1)=fW(nefW2,1)+Nf'*((+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg...
                                      +Voigt3*nx.*Lxyfg).*wfg);
      elseif nsd==3
        fW(nefW1,1)=fW(nefW1,1)+Nf'*((+Voigt1*nx.*Lxxfg+Voigt2*nx.*Lyyfg+Voigt2*nx.*Lzzfg...
                                      +Voigt3*ny.*Lxyfg+Voigt3*nz.*Lxzfg).*wfg);
        fW(nefW2,1)=fW(nefW2,1)+Nf'*((+Voigt2*ny.*Lxxfg+Voigt1*ny.*Lyyfg+Voigt2*ny.*Lzzfg...
                                      +Voigt3*nx.*Lxyfg+Voigt3*nz.*Lyzfg).*wfg);
        fW(nefW3,1)=fW(nefW3,1)+Nf'*((+Voigt2*nz.*Lxxfg+Voigt2*nz.*Lyyfg+Voigt1*nz.*Lzzfg...
                                      +Voigt3*nx.*Lxzfg+Voigt3*ny.*Lyzfg).*wfg);
      end
    end
    
    if not(isDirichlet)
      fW(nefW1,1)=fW(nefW1,1)+Nf'*(pfg.*nx.*wfg);
      fW(nefW2,1)=fW(nefW2,1)+Nf'*(pfg.*ny.*wfg);
      if nsd==3
        fW(nefW3,1)=fW(nefW3,1)+Nf'*(pfg.*nz.*wfg);
      end
    end
    
    if not(isDirichlet)
      fW(nefW1,1)=fW(nefW1,1)+Nf'*(tauW*wxfg.*wfg);
      fW(nefW2,1)=fW(nefW2,1)+Nf'*(tauW*wyfg.*wfg);
      if nsd==3
        fW(nefW3,1)=fW(nefW3,1)+Nf'*(tauW*wzfg.*wfg);
      end
    end
    
    fW(nefW1,1)=fW(nefW1,1)-Nf'*(tauW*Wxfg.*wfg);
    fW(nefW2,1)=fW(nefW2,1)-Nf'*(tauW*Wyfg.*wfg);
    if nsd==3
      fW(nefW3,1)=fW(nefW3,1)-Nf'*(tauW*Wzfg.*wfg);
    end
    
    if isDirichlet
      fW(nefW1,1)=fW(nefW1,1)+Nf'*(tauW*wDxfg.*wfg);
      fW(nefW2,1)=fW(nefW2,1)+Nf'*(tauW*wDyfg.*wfg);
      if nsd==3
        fW(nefW3,1)=fW(nefW3,1)+Nf'*(tauW*wDzfg.*wfg);
      end
    end
    
    if not(isDirichlet)
      fE(nefE1,1)=fE(nefE1,1)+Nf'*(Efg./Rfg.*(Wxfg.*nx+Wyfg.*ny).*wfg);
      if nsd==3
        fE(nefE1,1)=fE(nefE1,1)+Nf'*(Efg./Rfg.*Wzfg.*nz.*wfg);
      end
    end
    
    if not(isDirichlet)
      if nsd==2
        fE(nefE1,1)=fE(nefE1,1)...
                    +Nf'*((+(Wxfg./Rfg.*(Voigt1*Lxxfg+Voigt2*Lyyfg)...
                            +Wyfg./Rfg.*Voigt3.*Lxyfg).*nx...
                           +(Wyfg./Rfg.*(Voigt2*Lxxfg+Voigt1*Lyyfg)...
                            +Wxfg./Rfg.*Voigt3.*Lxyfg).*ny).*wfg);
      elseif nsd==3
        fE(nefE1,1)=fE(nefE1,1)...
                    +Nf'*((+(Wxfg./Rfg.*(Voigt1*Lxxfg+Voigt2*Lyyfg+Voigt2*Lzzfg)...
                            +Wyfg./Rfg.*Voigt3.*Lxyfg+Wzfg./Rfg.*Voigt3.*Lxzfg).*nx...
                           +(Wyfg./Rfg.*(Voigt2*Lxxfg+Voigt1*Lyyfg+Voigt2*Lzzfg)...
                            +Wxfg./Rfg.*Voigt3.*Lxyfg+Wzfg./Rfg.*Voigt3.*Lyzfg).*ny...
                           +(Wzfg./Rfg.*(Voigt2*Lxxfg+Voigt2*Lyyfg+Voigt1*Lzzfg)...
                            +Wxfg./Rfg.*Voigt3.*Lxzfg+Wyfg./Rfg.*Voigt3.*Lyzfg).*nz).*wfg);
      end
    end
    
    if not(isDirichlet)
      fE(nefE1,1)=fE(nefE1,1)+Nf'*(pfg./Rfg.*(Wxfg.*nx+Wyfg.*ny).*wfg);
      if nsd==3
        fE(nefE1,1)=fE(nefE1,1)+Nf'*(pfg./Rfg.*Wzfg.*nz.*wfg);
      end
    end
    
    if not(isDirichlet)
      fE(nefE1,1)=fE(nefE1,1)+Nf'*(sqrt(kappa)*(Qxfg.*nx+Qyfg.*ny).*wfg);
      if nsd==3
        fE(nefE1,1)=fE(nefE1,1)+Nf'*(sqrt(kappa)*(Qzfg.*nz).*wfg);
      end
    end
    
    if not(isDirichlet)
      fE(nefE1,1)=fE(nefE1,1)+Nf'*(tauE*efg.*wfg);
    end
    
    fE(nefE1,1)=fE(nefE1,1)-Nf'*(tauE*Efg.*wfg);
    
    if isDirichlet
      fE(nefE1,1)=fE(nefE1,1)+Nf'*(tauE*eDfg.*wfg);
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

% **************************************************************************************************
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
% **************************************************************************************************

% Lhs local-local
LhsLL(iL,[iL,   ir,iw   ])=[KLL,    KLr,KLw    ];
LhsLL(iQ,[   iQ,ir,iw,ie])=[    KQQ,KQr,KQw,KQe];
LhsLL(ir,[      ir,iw   ])=[        Krr,Krw    ];
LhsLL(iw,[iL,   ir,iw,ie])=[KwL,    Kwr,Kww,Kwe];
LhsLL(ie,[iL,iQ,ir,iw,ie])=[KeL,KeQ,Ker,Kew,Kee];

% Lhs local-global
LhsLG(iL,[iR,iW   ])=[KLR,KLW    ];
LhsLG(iQ,[iR,iW,iE])=[KQR,KQW,KQE];
LhsLG(ir,[iR,iW   ])=[KrR,KrW    ];
LhsLG(iw,[iR,iW,iE])=[KwR,KwW,KwE];
LhsLG(ie,[iR,iW,iE])=[KeR,KeW,KeE];

% Rhs local
RhsL(iL,1)=fL;
RhsL(iQ,1)=fQ;
RhsL(ir,1)=fr;
RhsL(iw,1)=fw;
RhsL(ie,1)=fe;

% Lhs global-local
LhsGL(iR,       ir       )=         KRr         ;
LhsGL(iW,[iL,      iw   ])=[KWL,        KWw    ];
LhsGL(iE,[iL,iQ,      ie])=[KEL,KEQ,        KEe];

% Lhs global-global
LhsGG(iR,[iR,iW   ])=[KRR,KRW    ];
LhsGG(iW,[iR,iW,iE])=[KWR,KWW,KWE];
LhsGG(iE,[iR,iW,iE])=[KER,KEW,KEE];

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

% Get Gauss data
gwe=RefElement.Post.GaussWeigthsElem;
gwf=RefElement.Post.GaussWeigthsFace;
nge=(1:length(gwe))';
ngf=(1:length(gwf))';

% Get reference data
FaceNodes=RefElement.FaceNodesElem;
N1e=ones(length(nge),1);
N1f=ones(length(ngf),1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions(1,RefElement.PostLow,RefElement.Post,Xe,nsd);
WEG=sparse(nge,nge,weg);

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
  Teg=T(reg,[wxeg,wyeg],eeg);
elseif nsd==3
  Teg=T(reg,[wxeg,wyeg,wzeg],eeg);
end

% Compute basic matrices
if nsd==2
  Kxxe=Nex'*(WEG*Nex); Kxye=Nex'*(WEG*Ney);
  Kyxe=Ney'*(WEG*Nex); Kyye=Ney'*(WEG*Ney);
elseif nsd==3
  Kxxe=Nex'*(WEG*Nex); Kxye=Nex'*(WEG*Ney); Kxze=Nex'*(WEG*Nez);
  Kyxe=Ney'*(WEG*Nex); Kyye=Ney'*(WEG*Ney); Kyze=Ney'*(WEG*Nez);
  Kzxe=Nez'*(WEG*Nex); Kzye=Nez'*(WEG*Ney); Kzze=Nez'*(WEG*Nez);
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

KVtp(1,ne1)=N1e'*(WEG*Ne);
KVtp(2,ne2)=N1e'*(WEG*Ne);
if nsd==3
  KVtp(3,ne3)=N1e'*(WEG*Ne);
end

if nsd==2
  KVrp(1,ne1)=-N1e'*(WEG*Ney);
  KVrp(1,ne2)=+N1e'*(WEG*Nex);
elseif nsd==3
  KVrp(1,ne2)=-N1e'*(WEG*Nez);
  KVrp(1,ne3)=+N1e'*(WEG*Ney);
  KVrp(2,ne1)=+N1e'*(WEG*Nez);
  KVrp(2,ne3)=-N1e'*(WEG*Nex);
  KVrp(3,ne1)=-N1e'*(WEG*Ney);
  KVrp(3,ne2)=+N1e'*(WEG*Nex);
end

KTpp(ne1,ne1)=-sqrt(kappa)*Kxxe...
              -sqrt(kappa)*Kyye;
if nsd==3
    KTpp(ne1,ne1)=KTpp(ne1,ne1)-sqrt(kappa)*Kzze;
end

KTtp(1,ne1)=N1e'*(WEG*Ne);

% Compute rhs
if nsd==2
  fVp(ne1,1)=Nex'*(Lxxeg.*weg)+Ney'*(Lxyeg.*weg);
  fVp(ne2,1)=Nex'*(Lxyeg.*weg)+Ney'*(Lyyeg.*weg);
elseif nsd==3
  fVp(ne1,1)=Nex'*(Lxxeg.*weg)+Ney'*(Lxyeg.*weg)+Nez'*(Lxzeg.*weg);
  fVp(ne2,1)=Nex'*(Lxyeg.*weg)+Ney'*(Lyyeg.*weg)+Nez'*(Lyzeg.*weg);
  fVp(ne3,1)=Nex'*(Lxzeg.*weg)+Ney'*(Lyzeg.*weg)+Nez'*(Lzzeg.*weg);
end

fVt(1,1)=N1e'*(wxeg./reg.*weg);
fVt(2,1)=N1e'*(wyeg./reg.*weg);
if nsd==3
  fVt(3,1)=N1e'*(wzeg./reg.*weg);
end

fTp(ne1,1)=Nex'*(Qxeg.*weg)...
          +Ney'*(Qyeg.*weg);
if nsd==3
    fTp(ne1,1)=fTp(ne1,1)+Nez'*(Qzeg.*weg);
end

fTt(1,1)=N1e'*(Teg.*weg);

% Faces loop
for iFace=1:NumElementFaces
  % Compute weights at Gauss points
  Xf=Xe(FaceNodes(iFace,:),:);
  [~,nx,ny,nz,wfg,Nlf]=mapShapeFunctions(0,RefElement.PostLow,RefElement.Post,Xf,nsd);
  
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
  
  % Compute rhs
  if nsd==2
    fVr(1)=fVr(1)+N1f'*((-Wxfg./Rfg.*ny+Wyfg./Rfg.*nx).*wfg);
  elseif nsd==3
    fVr(1)=fVr(1)+N1f'*((-Wyfg./Rfg.*nz+Wzfg./Rfg.*ny).*wfg);
    fVr(2)=fVr(2)+N1f'*((+Wxfg./Rfg.*nz-Wzfg./Rfg.*nx).*wfg);
    fVr(3)=fVr(3)+N1f'*((-Wxfg./Rfg.*ny+Wyfg./Rfg.*nx).*wfg);
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
LhsPost(iVp,iVp     )=KVpp     ;
LhsPost(iVt,iVp     )=KVtp     ;
LhsPost(iVr,iVp     )=KVrp     ;
LhsPost(iTp,    iTp2)=     KTpp;
LhsPost(iTt,    iTp2)=     KTtp;

% Rhs for post-processing
RhsPost(iVp,1)=fVp;
RhsPost(iVt,1)=fVt;
RhsPost(iVr,1)=fVr;
RhsPost(iTp,1)=fTp;
RhsPost(iTt,1)=fTt;

end