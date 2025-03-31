classdef Plasma1FluidElectromagnetic_HDG < Formulation  
  
  properties
    
    % Number of global components
    NumGlobalComp=@(NumSpaceDim) 1+NumSpaceDim+1+NumSpaceDim-1;
    
    % Number of Voigt components
    NumVoigtComp=@(NumSpaceDim) 0;
    
    % Number of local components
    NumLocalComp=@(NumSpaceDim) 1+NumSpaceDim+1+...
                                NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim+...
                                NumSpaceDim*(NumSpaceDim-1)/2+NumSpaceDim;
    
    % Number of post-processed components
    NumPostComp=@(NumSpaceDim) NumSpaceDim;

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
      
      if matchField(Time,'SteadySolution','yes')
        % Consider the values [1,0.1,...,0.1,1] for the fluid variables to avoid singularities
        Block(iD,iD).SolutionLocal(:,1:1+Sizes(iD).NumSpaceDim+1)=...
          [1,0.1*ones(1,Sizes(iD).NumSpaceDim),1].*ones(Sizes(iD).NumLocalNodes,1);
        for iFace=1:Sizes(iD).NumFaces
          Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp*Sizes(iD).NumFaceNodes+...
            (1:(1+Sizes(iD).NumSpaceDim+1)*Sizes(iD).NumFaceNodes))=...
            repelem([1,0.1*ones(1,Sizes(iD).NumSpaceDim),1],Sizes(iD).NumFaceNodes)';
        end
      else
        % Consider the values [1,0,...,0,1] for the fluid variables
        Block(iD,iD).SolutionLocal(:,1:1+Sizes(iD).NumSpaceDim+1)=...
          [1,0*ones(1,Sizes(iD).NumSpaceDim),1].*ones(Sizes(iD).NumLocalNodes,1);
        for iFace=1:Sizes(iD).NumFaces
          Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp*Sizes(iD).NumFaceNodes+...
            (1:(1+Sizes(iD).NumSpaceDim+1)*Sizes(iD).NumFaceNodes))=...
            repelem([1,0*ones(1,Sizes(iD).NumSpaceDim),1],Sizes(iD).NumFaceNodes)';
        end
      end
    end
    
    %% Compute initial conditions
    function [Block]=computeInitialConditions(~,iD,Block,Parameters,Mesh,Faces,Time,RefElement,...
        Sizes)
      if strcmp(Time.TimeDependent,'yes')
        if not(matchField(Time,'SteadySolution','yes'))
          Ne=RefElement(iD,iD).ShapeFunctionsElem;
          for iElem=1:Sizes(iD).NumElements
            Ce=Mesh(iD).Elements(:,iElem)';
            Xe=Mesh(iD).Nodes(:,Ce)';
            Xeg=Ne*Xe;
            Block(iD,iD).SolutionLocal((iElem-1)*Sizes(iD).NumElementNodes+...
              (1:Sizes(iD).NumElementNodes),:)=Ne\[...
              Parameters(iD).Density(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime),...
              Parameters(iD).Momentum(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime),...
              Parameters(iD).Energy(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime),...
              Parameters(iD).MagneticField(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime),...
              Parameters(iD).ElectricField(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime),...
              Parameters(iD).MagneticFieldAux(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime),...
              Parameters(iD).ElectricFieldAux(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.InitialTime)];
          end
          
          % Impose initial conditions on the fluid traces for a correct spectral decomposition
          for iFace=1:Sizes(iD).NumFaces
            [iElem,iElemFace]=find(Faces(iD,iD).Connectivity==iFace);
            [iElem,Order]=min(iElem);
            iElemFace=iElemFace(Order);
            FaceNodes=RefElement(iD,iD).FaceNodesElem;
            Ce=Mesh(iD).Elements(:,iElem)';
            Xe=Mesh(iD).Nodes(:,Ce)';
            Xf=Xe(FaceNodes(iElemFace,:),:);
            Block(iD,iD).SolutionGlobal((iFace-1)*Sizes(iD).NumGlobalComp*Sizes(iD).NumFaceNodes+...
            (1:(1+Sizes(iD).NumSpaceDim+1)*Sizes(iD).NumFaceNodes))=reshape(...
              [Parameters(iD).Density( Xf(:,1),Xf(:,2),Xf(:,3),Time.InitialTime),...
               Parameters(iD).Momentum(Xf(:,1),Xf(:,2),Xf(:,3),Time.InitialTime),...
               Parameters(iD).Energy(  Xf(:,1),Xf(:,2),Xf(:,3),Time.InitialTime)],[],1);
          end
        end
        Block(iD,iD).SolutionOld(:,:,1)=Block(iD,iD).SolutionLocal;
        
        % Get ghost solutions for BDF order > 2
        if Time.BDFOrder>2
          for iBDF=1:Time.BDFOrder-1
            Time.TimeOld=Time.Time-iBDF*Time.TimeStepSize;
            Ne=RefElement(iD,iD).ShapeFunctionsElem;
            for iElem=1:Sizes(iD).NumElements
              Ce=Mesh(iD).Elements(:,iElem)';
              Xe=Mesh(iD).Nodes(:,Ce)';
              Xeg=Ne*Xe;
              Block(iD,iD).SolutionOld((iElem-1)*Sizes(iD).NumElementNodes+...
                (1:Sizes(iD).NumElementNodes),:,1+iBDF)=Ne\[...
                Parameters(iD).Density(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld),...
                Parameters(iD).Momentum(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld),...
                Parameters(iD).Energy(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld),...
                Parameters(iD).MagneticField(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld),...
                Parameters(iD).ElectricField(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld),...
                Parameters(iD).MagneticFieldAux(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld),...
                Parameters(iD).ElectricFieldAux(Xeg(:,1),Xeg(:,2),Xeg(:,3),Time.TimeOld)];
            end
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
      Block(iD,iD).LhsGlobal=fsparse(Block(iD,iD).LhsRowIndices,....
                                     Block(iD,iD).LhsColIndices,LhsCoef(:));
      Block(iD,iD).RhsGlobal=fsparse(Block(iD,iD).RhsRowIndices,1,RhsCoef(:));
      Elements(iD).MatLocal=MatLocal;
      Elements(iD).VecLocal=VecLocal;
    end
    
    %% Do post-process
    function [Elements]=doPostProcess(~,Elements,~,Parameters,~,Time,RefElement,Sizes)
      RefElement.PostLowHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree  ,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+1,Parameters.Degree+2,Parameters.NodesDistribution);
      RefElement.PostHighHigh=createReferenceElement(Sizes.NumSpaceDim,...
                              Parameters.Degree+2,Parameters.Degree+2,Parameters.NodesDistribution);
      NodesElem=Elements.Nodes;
      SolutionLocalElem=Elements.SolutionLocal;
      SolutionOldElem=Elements.SolutionOld;
      LhsPost=cell(Sizes.NumElements,1);
      RhsPost=cell(Sizes.NumElements,1);
      parfor iElem=1:Sizes.NumElements
        [LhsPostElem,RhsPostElem]=...
          doPostProcessElement(NodesElem{iElem},...
          SolutionLocalElem{iElem},SolutionOldElem{iElem},...
          Parameters,Time,RefElement,Sizes);
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
        Results(iD).Density=[];
        Results(iD).Momentum=[];
        Results(iD).Energy=[];
        Results(iD).MagneticField=[];
        Results(iD).ElectricField=[];
        Results(iD).MagneticFieldAux=[];
        Results(iD).ElectricFieldAux=[];
      end
      Results(iD).Time(iST)=Time.Time;
      Results(iD).Density(:,:,iST)=Block(iD,iD).SolutionLocal(:,1);
      Results(iD).Momentum(:,:,iST)=Block(iD,iD).SolutionLocal(:,1+(1:Sizes(iD).NumSpaceDim));
      Results(iD).Energy(:,:,iST)=Block(iD,iD).SolutionLocal(:,1+Sizes(iD).NumSpaceDim+1);
      Results(iD).MagneticField(:,:,iST)=Block(iD,iD).SolutionLocal(:,1+Sizes(iD).NumSpaceDim+...
        1+(1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2));
      Results(iD).ElectricField(:,:,iST)=Block(iD,iD).SolutionLocal(:,1+Sizes(iD).NumSpaceDim+...
        1+Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+(1:Sizes(iD).NumSpaceDim));
      Results(iD).MagneticFieldAux(:,:,iST)=Block(iD,iD).SolutionLocal(:,1+Sizes(iD).NumSpaceDim+...
        1+Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim+...
        (1:Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2));
      Results(iD).ElectricFieldAux(:,:,iST)=Block(iD,iD).SolutionLocal(:,1+Sizes(iD).NumSpaceDim+...
        1+Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim+...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+(1:Sizes(iD).NumSpaceDim));
      if strcmp(Parameters(iD).PostProcessingHDG,'yes') && ...
         matchField(Block(iD,iD),'SolutionPost')
        Results(iD).ElectricFieldPost=Block(iD,iD).SolutionPost;
      end
    end
    
    %% Compute radiation
    function [Results]=computeRadiation(~,iD,Results,Elements,Mesh,Faces,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      ResultsFFTElem=mat2cell(...
        [Results(iD).MagneticFieldFFT(Mesh(iD).Elements(:),:),...
         Results(iD).ElectricFieldFFT(Mesh(iD).Elements(:),:)],...
        ones(Sizes(iD).NumElements,1)*Sizes(iD).NumElementNodes,...
        Sizes(iD).NumSpaceDim*(Sizes(iD).NumSpaceDim-1)/2+Sizes(iD).NumSpaceDim);
      Radiation=0;
      for iFace=1:size(Faces(iD,iD).Radiation,1)
        iElem=Faces(iD,iD).Radiation(iFace,1);
        iElemFace=Faces(iD,iD).Radiation(iFace,2);
        [RadiationElem]=...
          computeRadiationElement(iElemFace,NodesElem{iElem},ResultsFFTElem{iElem},...
          RefElement,Sizes);
        Radiation=Radiation+RadiationElem;
      end
      Results(iD).Radiation=Radiation;
    end
    
    %% Compute divergence
    function [Results]=computeDivergence(~,iD,Results,Elements,Parameters,RefElement,Sizes)
      NodesElem=Elements(iD).Nodes;
      SolutionLocalElem=Elements(iD).SolutionLocal;
      MagneticFieldDiv=zeros(Sizes(iD).NumElementNodes,Sizes(iD).NumElements);
      ElectricFieldDiv=zeros(Sizes(iD).NumElementNodes,Sizes(iD).NumElements);
      parfor iElem=1:Sizes(iD).NumElements
        [MagneticFieldDivElem,ElectricFielDivdElem]=...
          computeDivergenceElement(NodesElem{iElem},...
          SolutionLocalElem{iElem},...
          Parameters,RefElement,Sizes);
        MagneticFieldDiv(:,iElem)=reshape(MagneticFieldDivElem',[],1);
        ElectricFieldDiv(:,iElem)=reshape(ElectricFielDivdElem',[],1);
      end
      Results(iD).MagneticFieldDiv=MagneticFieldDiv(:);
      Results(iD).ElectricFieldDiv=ElectricFieldDiv(:);
    end
    
  end
  
end

%% Build block element
function [LhsGlobal,RhsGlobal,MatLocal,VecLocal]=buildBlockElement(...
  Nodes,Faces,SolutionGlobal,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
isAxisymmetric=strcmp(Parameters.Axisymmetric,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementFaces=Sizes.NumElementFaces;
NumFaceNodes=Sizes.NumFaceNodes;
gamma=Parameters.SpecificHeatRatio;
mu=Parameters.MagneticPermeability;
epsilon=Parameters.ElectricPermittivity;
sigma=Parameters.ElectricConductivity;
q1=Parameters.ParticleCharge; if strcmp(Parameters.Electromagnetic2FluidCoupling,'no'); q1=0; end
q2=Parameters.ParticleCharge; if strcmp(Parameters.Fluid2ElectromagneticCoupling,'no'); q2=0; end
m=Parameters.ParticleMass;
Sigma=Parameters.DampingFunction;
rD=Parameters.Density;
wD=Parameters.Momentum;
gD=Parameters.Energy;
eD=Parameters.ElectricField;
Fr=Parameters.ForceDensity;
Fw=Parameters.ForceMomentum;
Fg=Parameters.ForceEnergy;
j=Parameters.CurrentDensity;
i=Parameters.IncidentField;
Xe=Nodes';
tauE=Parameters.StabElectricField;
t=Time.Time;
dt=Time.TimeStepSize;
BDFo=Time.BDFOrderEff;
alpha=Time.BDF1stDerEff;

% Get solution
re=reshape(SolutionLocal(:,1),[],1);
we=reshape(SolutionLocal(:,1+(1:nsd)),[],1);
ge=reshape(SolutionLocal(:,1+nsd+1),[],1);
He=reshape(SolutionLocal(:,1+nsd+1+(1:qsd)),[],1);
ee=reshape(SolutionLocal(:,1+nsd+1+qsd+(1:nsd)),[],1);
Me=reshape(SolutionLocal(:,1+nsd+1+qsd+nsd+(1:qsd)),[],1);
Je=reshape(SolutionLocal(:,1+nsd+1+qsd+nsd+qsd+(1:nsd)),[],1);
Ue=SolutionGlobal;
rolde=reshape(SolutionOld(:,1,:),[],BDFo);
wolde=reshape(SolutionOld(:,1+(1:nsd),:),[],BDFo);
golde=reshape(SolutionOld(:,1+nsd+1,:),[],BDFo);
Holde=reshape(SolutionOld(:,1+nsd+1+(1:qsd),:),[],BDFo);
eolde=reshape(SolutionOld(:,1+nsd+1+qsd+(1:nsd),:),[],BDFo);
Molde=reshape(SolutionOld(:,1+nsd+1+qsd+nsd+(1:qsd),:),[],BDFo);
Jolde=reshape(SolutionOld(:,1+nsd+1+qsd+nsd+qsd+(1:nsd),:),[],BDFo);

% Initialize lhs
Krr=zeros(NumElementNodes,NumElementNodes);
Krw=zeros(NumElementNodes,nsd*NumElementNodes);
KrR=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KrW=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
Kwr=zeros(nsd*NumElementNodes,NumElementNodes);
Kww=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kwg=zeros(nsd*NumElementNodes,NumElementNodes);
KwH=zeros(nsd*NumElementNodes,qsd*NumElementNodes);
Kwe=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KwR=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
KwW=zeros(nsd*NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KwG=zeros(nsd*NumElementNodes,NumElementFaces*NumFaceNodes);
Kgr=zeros(NumElementNodes,NumElementNodes);
Kgw=zeros(NumElementNodes,nsd*NumElementNodes);
Kgg=zeros(NumElementNodes,NumElementNodes);
Kge=zeros(NumElementNodes,nsd*NumElementNodes);
KgR=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KgW=zeros(NumElementNodes,nsd*NumElementFaces*NumFaceNodes);
KgG=zeros(NumElementNodes,NumElementFaces*NumFaceNodes);
KHH=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KHe=zeros(qsd*NumElementNodes,nsd*NumElementNodes);
KHM=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KHE=zeros(qsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
Kew=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
Kee=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KeJ=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KeE=zeros(nsd*NumElementNodes,(nsd-1)*NumElementFaces*NumFaceNodes);
KMH=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KMM=zeros(qsd*NumElementNodes,qsd*NumElementNodes);
KJe=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KJJ=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
KRr=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KRw=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KRg=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KRR=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KRW=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KRG=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KWr=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementNodes);
KWw=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KWg=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementNodes);
KWR=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KWW=zeros(nsd*NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KWG=zeros(nsd*NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KGr=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KGw=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KGg=zeros(NumElementFaces*NumFaceNodes,NumElementNodes);
KGR=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KGW=zeros(NumElementFaces*NumFaceNodes,nsd*NumElementFaces*NumFaceNodes);
KGG=zeros(NumElementFaces*NumFaceNodes,NumElementFaces*NumFaceNodes);
KEH=zeros((nsd-1)*NumElementFaces*NumFaceNodes,qsd*NumElementNodes);
KEe=zeros((nsd-1)*NumElementFaces*NumFaceNodes,nsd*NumElementNodes);
KEE=zeros((nsd-1)*NumElementFaces*NumFaceNodes,(nsd-1)*NumElementFaces*NumFaceNodes);

% Initialize rhs
fr=zeros(NumElementNodes,1);
fw=zeros(nsd*NumElementNodes,1);
fg=zeros(NumElementNodes,1);
fH=zeros(qsd*NumElementNodes,1);
fe=zeros(nsd*NumElementNodes,1);
fM=zeros(qsd*NumElementNodes,1);
fJ=zeros(nsd*NumElementNodes,1);
fR=zeros(NumElementFaces*NumFaceNodes,1);
fW=zeros(nsd*NumElementFaces*NumFaceNodes,1);
fG=zeros(NumElementFaces*NumFaceNodes,1);
fE=zeros((nsd-1)*NumElementFaces*NumFaceNodes,1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,weg]=mapShapeFunctions(1,RefElement,RefElement,Xe(1:nsd+1,:),nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
wxe=we(ne1);
wye=we(ne2);
if nsd==3
  wze=we(ne3);
end
if nsd==2
  Hze=He(ne1);
elseif nsd==3
  Hxe=He(ne1);
  Hye=He(ne2);
  Hze=He(ne3);
end
exe=ee(ne1);
eYe=ee(ne2);
if nsd==3
  eze=ee(ne3);
end
if nsd==2
  Mze=Me(ne1);
elseif nsd==3
  Mxe=Me(ne1);
  Mye=Me(ne2);
  Mze=Me(ne3);
end
Jxe=Je(ne1);
Jye=Je(ne2);
if nsd==3
  Jze=Je(ne3);
end
woldxe=wolde(ne1,:);
woldye=wolde(ne2,:);
if nsd==3
  woldze=wolde(ne3,:);
end
if nsd==2
  Holdze=Holde(ne1,:);
elseif nsd==3
  Holdxe=Holde(ne1,:);
  Holdye=Holde(ne2,:);
  Holdze=Holde(ne3,:);
end
eoldxe=eolde(ne1,:);
eoldye=eolde(ne2,:);
if nsd==3
  eoldze=eolde(ne3,:);
end
if nsd==2
  Moldze=Molde(ne1,:);
elseif nsd==3
  Moldxe=Molde(ne1,:);
  Moldye=Molde(ne2,:);
  Moldze=Molde(ne3,:);
end
Joldxe=Jolde(ne1,:);
Joldye=Jolde(ne2,:);
if nsd==3
  Joldze=Jolde(ne3,:);
end
Sigmae=Sigma(Xe(:,1),Xe(:,2),Xe(:,3));
sigmaxe=Sigmae(:,1);
sigmaye=Sigmae(:,2);
if nsd==3
  sigmaze=Sigmae(:,3);
end
isPML=max(max(abs(Sigmae)))>1e-9;

% Remove coupling inside PML region
if isPML
  q1=0;
  q2=0;
end

% Compute variables at Gauss points
Xeg=Ne*Xe;
reg=Ne*re;
wxeg=Ne*wxe;
wyeg=Ne*wye;
if nsd==3
  wzeg=Ne*wze;
end
geg=Ne*ge;
if nsd==2
  Hzeg=Ne*Hze;
elseif nsd==3
  Hxeg=Ne*Hxe;
  Hyeg=Ne*Hye;
  Hzeg=Ne*Hze;
end
exeg=Ne*exe;
eyeg=Ne*eYe;
if nsd==3
  ezeg=Ne*eze;
end
if nsd==2
  Mzeg=Ne*Mze;
elseif nsd==3
  Mxeg=Ne*Mxe;
  Myeg=Ne*Mye;
  Mzeg=Ne*Mze;
end
Jxeg=Ne*Jxe;
Jyeg=Ne*Jye;
if nsd==3
  Jzeg=Ne*Jze;
end
roldeg=Ne*rolde;
woldxeg=Ne*woldxe;
woldyeg=Ne*woldye;
if nsd==3
  woldzeg=Ne*woldze;
end
goldeg=Ne*golde;
if nsd==2
  Holdzeg=Ne*Holdze;
elseif nsd==3
  Holdxeg=Ne*Holdxe;
  Holdyeg=Ne*Holdye;
  Holdzeg=Ne*Holdze;
end
eoldxeg=Ne*eoldxe;
eoldyeg=Ne*eoldye;
if nsd==3
  eoldzeg=Ne*eoldze;
end
if nsd==2
  Moldzeg=Ne*Moldze;
elseif nsd==3
  Moldxeg=Ne*Moldxe;
  Moldyeg=Ne*Moldye;
  Moldzeg=Ne*Moldze;
end
Joldxeg=Ne*Joldxe;
Joldyeg=Ne*Joldye;
if nsd==3
  Joldzeg=Ne*Joldze;
end
sigmaxeg=Ne*sigmaxe;
sigmayeg=Ne*sigmaye;
if nsd==3
  sigmazeg=Ne*sigmaze;
end
Freg=Fr(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
Fweg=Fw(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
Fwxeg=Fweg(:,1);
Fwyeg=Fweg(:,2);
if nsd==3
  Fwzeg=Fweg(:,3);
end
Fgeg=Fg(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
jeg=j(Xeg(:,1),Xeg(:,2),Xeg(:,3),t);
jxeg=jeg(:,1);
jyeg=jeg(:,2);
if nsd==3
  jzeg=jeg(:,3);
end
if isa(mu,'function_handle')
  mueg=mu(Xeg(:,1),Xeg(:,2),Xeg(:,3));
  epsiloneg=epsilon(Xeg(:,1),Xeg(:,2),Xeg(:,3));
  sigmaeg=sigma(Xeg(:,1),Xeg(:,2),Xeg(:,3));
else
  mueg=mu;
  epsiloneg=epsilon;
  sigmaeg=sigma;
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

% Compute terms for speedup
vxeg=wxeg./reg;
vyeg=wyeg./reg;
if nsd==3
  vzeg=wzeg./reg;
end
if nsd==2
  v2eg=vxeg.^2+vyeg.^2;
elseif nsd==3
  v2eg=vxeg.^2+vyeg.^2+vzeg.^2;
end
g_reg=geg./reg;
peg=(gamma-1)*(geg-1/2*reg.*v2eg);
if nsd==2
  wHxeg=+wyeg.*Hzeg;
  wHyeg=-wxeg.*Hzeg;
elseif nsd==3
  wHxeg=+wyeg.*Hzeg-wzeg.*Hyeg;
  wHyeg=+wzeg.*Hxeg-wxeg.*Hzeg;
  wHzeg=+wxeg.*Hyeg-wyeg.*Hxeg;
end
if nsd==2
  weeg=wxeg.*exeg+wyeg.*eyeg;
elseif nsd==3
  weeg=wxeg.*exeg+wyeg.*eyeg+wzeg.*ezeg;
end
if nsd==2
  Sigma1xeg=sigmayeg-sigmaxeg;
  Sigma1yeg=sigmaxeg-sigmayeg;
  Sigma1zeg=sigmaxeg+sigmayeg;
elseif nsd==3
  Sigma1xeg=sigmayeg+sigmazeg-sigmaxeg;
  Sigma1yeg=sigmazeg+sigmaxeg-sigmayeg;
  Sigma1zeg=sigmaxeg+sigmayeg-sigmazeg;
end
Sigma2xeg=sigmaxeg;
Sigma2yeg=sigmayeg;
if nsd==3
  Sigma2zeg=sigmazeg;
end
if nsd==2
  Sigma3xeg=(sigmaxeg-sigmayeg).*(sigmaxeg);
  Sigma3yeg=(sigmayeg-sigmaxeg).*(sigmayeg);
  Sigma3zeg=(-sigmaxeg).*(-sigmayeg);
elseif nsd==3
  Sigma3xeg=(sigmaxeg-sigmayeg).*(sigmaxeg-sigmazeg);
  Sigma3yeg=(sigmayeg-sigmaxeg).*(sigmayeg-sigmazeg);
  Sigma3zeg=(sigmazeg-sigmaxeg).*(sigmazeg-sigmayeg);
end
if nsd==2
  gradWHxeg=-Ney*Hze;
  gradWHyeg=+Nex*Hze;
elseif nsd==3
  gradWHxeg=+Nez*Hye-Ney*Hze;
  gradWHyeg=+Nex*Hze-Nez*Hxe;
  gradWHzeg=+Ney*Hxe-Nex*Hye;
end

% Compute lhs
Krr(ne1,ne1)=alpha(1)/dt*Me;

Krw(ne1,ne1)=-Cxe;
Krw(ne1,ne2)=-Cye;
if nsd==3
  Krw(ne1,ne3)=-Cze;
end

if nsd==2
  Kww(ne1,ne1)=+alpha(1)/dt*Me...
               -NwexT*(((2-(gamma-1))*vxeg).*Ne)...
               -NweyT*((              vyeg).*Ne);
  Kww(ne2,ne1)=-NwexT*((              vyeg).*Ne)...
               -NweyT*((   -(gamma-1)*vxeg).*Ne)...
               -NweT*((-q1/m*mueg.*Hzeg).*Ne);
  Kww(ne1,ne2)=-NwexT*((   -(gamma-1)*vyeg).*Ne)...
               -NweyT*((              vxeg).*Ne)...
               -NweT*((+q1/m*mueg.*Hzeg).*Ne);
  Kww(ne2,ne2)=+alpha(1)/dt*Me...
               -NwexT*((              vxeg).*Ne)...
               -NweyT*(((2-(gamma-1))*vyeg).*Ne);
elseif nsd==3
  Kww(ne1,ne1)=+alpha(1)/dt*Me...
               -NwexT*(( (2-(gamma-1))*vxeg).*Ne)...
               -NweyT*((               vyeg).*Ne)...
               -NwezT*((               vzeg).*Ne);
  Kww(ne2,ne1)=-NwexT*((               vyeg).*Ne)...
               -NweyT*((    -(gamma-1)*vxeg).*Ne)...
               -NweT*((-q1/m*mueg.*Hzeg).*Ne);
  Kww(ne3,ne1)=-NwezT*((vzeg-(gamma-1)*vxeg).*Ne)...
               -NweT*((+q1/m*mueg.*Hyeg).*Ne);
  Kww(ne1,ne2)=-NwexT*((    -(gamma-1)*vyeg).*Ne)...
               -NweyT*((               vxeg).*Ne)...
               -NweT*((+q1/m*mueg.*Hzeg).*Ne);
  Kww(ne2,ne2)=+alpha(1)/dt*Me...
               -NwexT*((               vxeg).*Ne)...
               -NweyT*(( (2-(gamma-1))*vyeg).*Ne)...
               -NwezT*((               vzeg).*Ne);
  Kww(ne3,ne2)=-NwezT*((vzeg-(gamma-1)*vyeg).*Ne)...
               -NweT*((-q1/m*mueg.*Hxeg).*Ne);
  Kww(ne1,ne3)=-NwexT*((    -(gamma-1)*vzeg).*Ne)...
               -NwezT*((               vxeg).*Ne)...
               -NweT*((-q1/m*mueg.*Hyeg).*Ne);
  Kww(ne2,ne3)=-NweyT*((    -(gamma-1)*vzeg).*Ne)...
               -NwezT*((               vyeg).*Ne)...
               -NweT*((+q1/m*mueg.*Hxeg).*Ne);
  Kww(ne3,ne3)=+alpha(1)/dt*Me...
               -NwexT*((               vxeg).*Ne)...
               -NweyT*((               vyeg).*Ne)...
               -NwezT*(( (2-(gamma-1))*vzeg).*Ne);
end

if nsd==2
  Kwr(ne1,ne1)=+NwexT*((vxeg.*vxeg-(gamma-1)*1/2*v2eg).*Ne)...
               +NweyT*((vxeg.*vyeg).*Ne)...
               -NweT*((q1/m*exeg).*Ne);
  Kwr(ne2,ne1)=+NwexT*((vyeg.*vxeg).*Ne)...
               +NweyT*((vyeg.*vyeg-(gamma-1)*1/2*v2eg).*Ne)...
               -NweT*((q1/m*eyeg).*Ne);
elseif nsd==3
  Kwr(ne1,ne1)=+NwexT*((vxeg.*vxeg-(gamma-1)*1/2*v2eg).*Ne)...
               +NweyT*((vxeg.*vyeg).*Ne)...
               +NwezT*((vxeg.*vzeg).*Ne)...
               -NweT*((q1/m*exeg).*Ne);
  Kwr(ne2,ne1)=+NwexT*((vyeg.*vxeg).*Ne)...
               +NweyT*((vyeg.*vyeg-(gamma-1)*1/2*v2eg).*Ne)...
               +NwezT*((vyeg.*vzeg).*Ne)...
               -NweT*((q1/m*eyeg).*Ne);
  Kwr(ne3,ne1)=+NwexT*((vzeg.*vxeg).*Ne)...
               +NweyT*((vzeg.*vyeg).*Ne)...
               +NwezT*((vzeg.*vzeg-(gamma-1)*1/2*v2eg).*Ne)...
               -NweT*((q1/m*ezeg).*Ne);
end

Kwg(ne1,ne1)=-(gamma-1)*Cxe;
Kwg(ne2,ne1)=-(gamma-1)*Cye;
if nsd==3
  Kwg(ne3,ne1)=-(gamma-1)*Cze;
end

Kwe(ne1,ne1)=-NweT*((q1/m*reg).*Ne);
Kwe(ne2,ne2)=-NweT*((q1/m*reg).*Ne);
if nsd==3
  Kwe(ne3,ne3)=-NweT*((q1/m*reg).*Ne);
end

if nsd==2
  KwH(ne1,ne1)=-NweT*((+q1/m*mueg.*wyeg).*Ne);
  KwH(ne2,ne1)=-NweT*((-q1/m*mueg.*wxeg).*Ne);  
elseif nsd==3
  KwH(ne2,ne1)=-NweT*((+q1/m*mueg.*wzeg).*Ne);
  KwH(ne3,ne1)=-NweT*((-q1/m*mueg.*wyeg).*Ne);
  KwH(ne1,ne2)=-NweT*((-q1/m*mueg.*wzeg).*Ne);
  KwH(ne3,ne2)=-NweT*((+q1/m*mueg.*wxeg).*Ne);
  KwH(ne1,ne3)=-NweT*((+q1/m*mueg.*wyeg).*Ne);
  KwH(ne2,ne3)=-NweT*((-q1/m*mueg.*wxeg).*Ne);
end

Kgg(ne1,ne1)=+alpha(1)/dt*Me...
             -NwexT*((gamma*vxeg).*Ne)...
             -NweyT*((gamma*vyeg).*Ne);
if nsd==3
  Kgg(ne1,ne1)=Kgg(ne1,ne1)-NwezT*((gamma*vzeg).*Ne);
end

Kgr(ne1,ne1)=+NwexT*(((g_reg-(gamma-1)*(-g_reg+v2eg)).*vxeg).*Ne)...
             +NweyT*(((g_reg-(gamma-1)*(-g_reg+v2eg)).*vyeg).*Ne);
if nsd==3
  Kgr(ne1,ne1)=Kgr(ne1,ne1)+NwezT*(((g_reg-(gamma-1)*(-g_reg+v2eg)).*vzeg).*Ne);
end

if nsd==2
  Kgw(ne1,ne1)=-NwexT*((g_reg+(gamma-1)*(g_reg-1/2*(3*vxeg.^2+vyeg.^2))).*Ne)...
               -NweyT*((      (gamma-1)*(-vxeg.*vyeg)).*Ne)...
               -NweT*((q1/m*exeg).*Ne);
  Kgw(ne1,ne2)=-NwexT*((      (gamma-1)*(-vyeg.*vxeg)).*Ne)...
               -NweyT*((g_reg+(gamma-1)*(g_reg-1/2*(vxeg.^2+3*vyeg.^2))).*Ne)...
               -NweT*((q1/m*eyeg).*Ne);
elseif nsd==3
  Kgw(ne1,ne1)=-NwexT*((g_reg+(gamma-1)*(g_reg-1/2*(3*vxeg.^2+vyeg.^2+vzeg.^2))).*Ne)...
               -NweyT*((      (gamma-1)*(-vxeg.*vyeg)).*Ne)...
               -NwezT*((      (gamma-1)*(-vxeg.*vzeg)).*Ne)...
               -NweT*((q1/m*exeg).*Ne);
  Kgw(ne1,ne2)=-NwexT*((      (gamma-1)*(-vyeg.*vxeg)).*Ne)...
               -NweyT*((g_reg+(gamma-1)*(g_reg-1/2*(vxeg.^2+3*vyeg.^2+vzeg.^2))).*Ne)...
               -NwezT*((      (gamma-1)*(-vyeg.*vzeg)).*Ne)...
               -NweT*((q1/m*eyeg).*Ne);
  Kgw(ne1,ne3)=-NwexT*((      (gamma-1)*(-vzeg.*vxeg)).*Ne)...
               -NweyT*((      (gamma-1)*(-vzeg.*vyeg)).*Ne)...
               -NwezT*((g_reg+(gamma-1)*(g_reg-1/2*(vxeg.^2+vyeg.^2+3*vzeg.^2))).*Ne)...
               -NweT*((q1/m*ezeg).*Ne);
end

Kge(ne1,ne1)=-NweT*((q1/m*wxeg).*Ne);
Kge(ne1,ne2)=-NweT*((q1/m*wyeg).*Ne);
if nsd==3
  Kge(ne1,ne3)=-NweT*((q1/m*wzeg).*Ne);
end

if nsd==2
  KHH(ne1,ne1)=-NweT*((mueg.*(alpha(1)/dt+Sigma1zeg)).*Ne);
elseif nsd==3
  KHH(ne1,ne1)=-NweT*((mueg.*(alpha(1)/dt+Sigma1xeg)).*Ne);
  KHH(ne2,ne2)=-NweT*((mueg.*(alpha(1)/dt+Sigma1yeg)).*Ne);
  KHH(ne3,ne3)=-NweT*((mueg.*(alpha(1)/dt+Sigma1zeg)).*Ne);
end

if nsd==2
  KHe(ne1,ne1)=-Cye;
  KHe(ne1,ne2)=+Cxe;
elseif nsd==3
  KHe(ne2,ne1)=+Cze;
  KHe(ne3,ne1)=-Cye;
  KHe(ne1,ne2)=-Cze;
  KHe(ne3,ne2)=+Cxe;
  KHe(ne1,ne3)=+Cye;
  KHe(ne2,ne3)=-Cxe;
end

if isPML
  KHM(ne1,ne1)=-Me;
  if nsd==3
    KHM(ne2,ne2)=-Me;
    KHM(ne3,ne3)=-Me;
  end
end

Kee(ne1,ne1)=NweT*((epsiloneg.*(alpha(1)/dt+Sigma1xeg)+sigmaeg).*Ne);
Kee(ne2,ne2)=NweT*((epsiloneg.*(alpha(1)/dt+Sigma1yeg)+sigmaeg).*Ne);
if nsd==3
  Kee(ne3,ne3)=NweT*((epsiloneg.*(alpha(1)/dt+Sigma1zeg)+sigmaeg).*Ne);
end

if isPML
  KeJ(ne1,ne1)=Me;
  KeJ(ne2,ne2)=Me;
  if nsd==3
    KeJ(ne3,ne3)=Me;
  end
end

Kew(ne1,ne1)=q2/m*Me;
Kew(ne2,ne2)=q2/m*Me;
if nsd==3
  Kew(ne3,ne3)=q2/m*Me;
end

if isPML
  if nsd==2
    KMM(ne1,ne1)=alpha(1)/dt*Me;
  elseif nsd==3
    KMM(ne1,ne1)=alpha(1)/dt*Me+NweT*((Sigma2xeg).*Ne);
    KMM(ne2,ne2)=alpha(1)/dt*Me+NweT*((Sigma2yeg).*Ne);
    KMM(ne3,ne3)=alpha(1)/dt*Me+NweT*((Sigma2zeg).*Ne);
  end
  
  if nsd==2
    KMH(ne1,ne1)=-NweT*((mueg.*Sigma3zeg).*Ne);
  elseif nsd==3
    KMH(ne1,ne1)=-NweT*((mueg.*Sigma3xeg).*Ne);
    KMH(ne2,ne2)=-NweT*((mueg.*Sigma3yeg).*Ne);
    KMH(ne3,ne3)=-NweT*((mueg.*Sigma3zeg).*Ne);
  end
  
  KJJ(ne1,ne1)=-alpha(1)/dt*Me-NweT*((Sigma2xeg).*Ne);
  KJJ(ne2,ne2)=-alpha(1)/dt*Me-NweT*((Sigma2yeg).*Ne);
  if nsd==3
    KJJ(ne3,ne3)=-alpha(1)/dt*Me-NweT*((Sigma2zeg).*Ne);
  end
  
  KJe(ne1,ne1)=+NweT*((epsiloneg.*Sigma3xeg).*Ne);
  KJe(ne2,ne2)=+NweT*((epsiloneg.*Sigma3yeg).*Ne);
  if nsd==3
    KJe(ne3,ne3)=+NweT*((epsiloneg.*Sigma3zeg).*Ne);
  end
end

% Compute rhs
fr(ne1,1)=+NweT*(-[reg,roldeg]*alpha/dt+Freg)...
          +NwexT*(wxeg)...
          +NweyT*(wyeg);
if nsd==3
  fr(ne1,1)=fr(ne1,1)+NwezT*(wzeg);
end

if nsd==2
  fw(ne1,1)=+NweT*(-[wxeg,woldxeg]*alpha/dt+q1/m*(reg.*exeg+mueg.*wHxeg)+Fwxeg)...
            +NwexT*(wxeg.*vxeg+peg)+NweyT*(wxeg.*vyeg    );
  fw(ne2,1)=+NweT*(-[wyeg,woldyeg]*alpha/dt+q1/m*(reg.*eyeg+mueg.*wHyeg)+Fwyeg)...
            +NwexT*(wyeg.*vxeg    )+NweyT*(wyeg.*vyeg+peg);
elseif nsd==3
  fw(ne1,1)=+NweT*(-[wxeg,woldxeg]*alpha/dt+q1/m*(reg.*exeg+mueg.*wHxeg)+Fwxeg)...
            +NwexT*(wxeg.*vxeg+peg)+NweyT*(wxeg.*vyeg    )+NwezT*(wxeg.*vzeg    );
  fw(ne2,1)=+NweT*(-[wyeg,woldyeg]*alpha/dt+q1/m*(reg.*eyeg+mueg.*wHyeg)+Fwyeg)...
            +NwexT*(wyeg.*vxeg    )+NweyT*(wyeg.*vyeg+peg)+NwezT*(wyeg.*vzeg    );
  fw(ne3,1)=+NweT*(-[wzeg,woldzeg]*alpha/dt+q1/m*(reg.*ezeg+mueg.*wHzeg)+Fwzeg)...
            +NwexT*(wzeg.*vxeg    )+NweyT*(wzeg.*vyeg    )+NwezT*(wzeg.*vzeg+peg);
end

fg(ne1,1)=+NweT*(-[geg,goldeg]*alpha/dt+q1/m*weeg+Fgeg)...
          +NwexT*((geg+peg).*vxeg)...
          +NweyT*((geg+peg).*vyeg);
if nsd==3
  fg(ne1,1)=fg(ne1,1)+NwezT*((geg+peg).*vzeg);
end

if nsd==2
  fH(ne1,1)=+NweT*(+mueg.*[Hzeg,Holdzeg]*alpha/dt+mueg.*Sigma1zeg.*Hzeg+Mzeg)...
            -NwexT*(eyeg)...
            +NweyT*(exeg);
elseif nsd==3
  fH(ne1,1)=+NweT*(+mueg.*[Hxeg,Holdxeg]*alpha/dt+mueg.*Sigma1xeg.*Hxeg+Mxeg)...
            -NweyT*(ezeg)...
            +NwezT*(eyeg);
  fH(ne2,1)=+NweT*(+mueg.*[Hyeg,Holdyeg]*alpha/dt+mueg.*Sigma1yeg.*Hyeg+Myeg)...
            -NwezT*(exeg)...
            +NwexT*(ezeg);
  fH(ne3,1)=+NweT*(+mueg.*[Hzeg,Holdzeg]*alpha/dt+mueg.*Sigma1zeg.*Hzeg+Mzeg)...
            -NwexT*(eyeg)...
            +NweyT*(exeg);
end

fe(ne1,1)=NweT*(-epsiloneg.*[exeg,eoldxeg]*alpha/dt-sigmaeg.*exeg-gradWHxeg...
                -epsiloneg.*Sigma1xeg.*exeg-Jxeg-jxeg-q2/m*wxeg);
fe(ne2,1)=NweT*(-epsiloneg.*[eyeg,eoldyeg]*alpha/dt-sigmaeg.*eyeg-gradWHyeg...
                -epsiloneg.*Sigma1yeg.*eyeg-Jyeg-jyeg-q2/m*wyeg);
if nsd==3
  fe(ne3,1)=NweT*(-epsiloneg.*[ezeg,eoldzeg]*alpha/dt-sigmaeg.*ezeg-gradWHzeg...
                  -epsiloneg.*Sigma1zeg.*ezeg-Jzeg-jzeg-q2/m*wzeg);
end

if isPML
  if nsd==2
    fM(ne1,1)=NweT*(-[Mzeg,Moldzeg]*alpha/dt+mueg.*Sigma3zeg.*Hzeg);
  elseif nsd==3
    fM(ne1,1)=NweT*(-[Mxeg,Moldxeg]*alpha/dt-Sigma2xeg.*Mxeg+mueg.*Sigma3xeg.*Hxeg);
    fM(ne2,1)=NweT*(-[Myeg,Moldyeg]*alpha/dt-Sigma2yeg.*Myeg+mueg.*Sigma3yeg.*Hyeg);
    fM(ne3,1)=NweT*(-[Mzeg,Moldzeg]*alpha/dt-Sigma2zeg.*Mzeg+mueg.*Sigma3zeg.*Hzeg);
  end
  
  fJ(ne1,1)=NweT*(+[Jxeg,Joldxeg]*alpha/dt+Sigma2xeg.*Jxeg-epsiloneg.*Sigma3xeg.*exeg);
  fJ(ne2,1)=NweT*(+[Jyeg,Joldyeg]*alpha/dt+Sigma2yeg.*Jyeg-epsiloneg.*Sigma3yeg.*eyeg);
  if nsd==3
    fJ(ne3,1)=NweT*(+[Jzeg,Joldzeg]*alpha/dt+Sigma2zeg.*Jzeg-epsiloneg.*Sigma3zeg.*ezeg);
  end
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
    [Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf(1:nsd,:),nsd);
    
    % Check boundary
    isInterior=Faces.Interior(1,iFace);
    isDirichlet_r=Faces.Dirichlet_r(iFace);
    isDirichlet_w=Faces.Dirichlet_w(iFace);
    isDirichlet_g=Faces.Dirichlet_g(iFace);
    isDirichlet_e=Faces.Dirichlet_e(iFace);
    isFarField=Faces.FarField(iFace);
    isInviscidWall=Faces.InviscidWall(iFace);
    isAbsorbing=Faces.Absorbing(iFace);
    
    % Indices
    nf1=FaceNodes(iFace,:);
    nf2=nf1+NumElementNodes;
    nf3=nf2+NumElementNodes;
    nefU1=(iFace-1)*(1+nsd+1+nsd-1)*NumFaceNodes+(1:NumFaceNodes);
    nefU2=nefU1+NumFaceNodes;
    nefU3=nefU2+NumFaceNodes;
    nefU4=nefU3+NumFaceNodes;
    nefU5=nefU4+NumFaceNodes;
    nefU6=nefU5+NumFaceNodes;
    nefU7=nefU6+NumFaceNodes;
    nefR1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefW1=(iFace-1)*nsd*NumFaceNodes+(1:NumFaceNodes);
    nefW2=nefW1+NumFaceNodes;
    nefW3=nefW2+NumFaceNodes;
    nefG1=(iFace-1)*NumFaceNodes+(1:NumFaceNodes);
    nefE1=(iFace-1)*(nsd-1)*NumFaceNodes+(1:NumFaceNodes);
    nefE2=nefE1+NumFaceNodes;
    nefR=1:nsd;
    
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
      nefU6=nefU6(order);
      nefU7=nefU7(order);
      nefR1=nefR1(order);
      nefW1=nefW1(order);
      nefW2=nefW2(order);
      nefW3=nefW3(order);
      nefG1=nefG1(order);
      nefE1=nefE1(order);
      nefE2=nefE2(order);
      nefR=circshift(flip(nefR(1:nsd)),Node2Match1stNode1);
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
        nefU3=nefU3(order);
        nefU4=nefU4(order);
        nefU5=nefU5(order);
        nefU6=nefU6(order);
        nefU7=nefU7(order);
        nefR1=nefR1(order);
        nefW1=nefW1(order);
        nefW2=nefW2(order);
        nefW3=nefW3(order);
        nefG1=nefG1(order);
        nefE1=nefE1(order);
        nefE2=nefE2(order);
        nefR=circshift(flip(nefR(1:nsd)),Node2Match1stNode1);
      end
    end
    
    % Gram-Schmidt orthonormalization (straight faces)
    V=Xf(nefR,1:nsd)';
    W=(V(:,2:end)-V(:,1))./vecnorm(V(:,2:end)-V(:,1));
    R(:,1)=W(:,1);
    if nsd==3
      R(:,2)=(W(:,2)-(R(:,1)'*W(:,2))*R(:,1))/norm(W(:,2)-(R(:,1)'*W(:,2))*R(:,1));
    end
    
    % Get rotation operator
    if nsd==2
      Rxt1=R(1,1);
      Ryt1=R(2,1);
    elseif nsd==3
      Rxt1=R(1,1); Rxt2=R(1,2);
      Ryt1=R(2,1); Ryt2=R(2,2);
      Rzt1=R(3,1); Rzt2=R(3,2);
    end
    
    % Compute variables at nodes
    rf=re(nf1);
    wxf=wxe(nf1);
    wyf=wye(nf1);
    if nsd==3
      wzf=wze(nf1);
    end
    gf=ge(nf1);
    if nsd==2
      Hzf=Hze(nf1);
    elseif nsd==3
      Hxf=Hxe(nf1);
      Hyf=Hye(nf1);
      Hzf=Hze(nf1);
    end
    exf=exe(nf1);
    eyf=eYe(nf1);
    if nsd==3
      ezf=eze(nf1);
    end
    Rf=Ue(nefU1);
    Wxf=Ue(nefU2);
    Wyf=Ue(nefU3);
    if nsd==3
      Wzf=Ue(nefU4);
    end
    if nsd==2
      Gf=Ue(nefU4);
    elseif nsd==3
      Gf=Ue(nefU5);
    end
    if nsd==2
      Et1f=Ue(nefU5);
    elseif nsd==3
      Et1f=Ue(nefU6);
      Et2f=Ue(nefU7);
    end
    
    % Compute variables at Gauss points
    Xfg=Nf*Xf;
    rfg=Nf*rf;
    wxfg=Nf*wxf;
    wyfg=Nf*wyf;
    if nsd==3
      wzfg=Nf*wzf;
    end
    gfg=Nf*gf;
    if nsd==2
      Hzfg=Nf*Hzf;
    elseif nsd==3
      Hxfg=Nf*Hxf;
      Hyfg=Nf*Hyf;
      Hzfg=Nf*Hzf;
    end
    exfg=Nf*exf;
    eyfg=Nf*eyf;
    if nsd==3
      ezfg=Nf*ezf;
    end
    Rfg=Nf*Rf;
    Wxfg=Nf*Wxf;
    Wyfg=Nf*Wyf;
    if nsd==3
      Wzfg=Nf*Wzf;
    end
    Gfg=Nf*Gf;
    Et1fg=Nf*Et1f;
    if nsd==3
      Et2fg=Nf*Et2f;
    end
    if nsd==2
      Etxfg=Rxt1*Et1fg;
      Etyfg=Ryt1*Et1fg;
    elseif nsd==3
      Etxfg=Rxt1*Et1fg+Rxt2*Et2fg;
      Etyfg=Ryt1*Et1fg+Ryt2*Et2fg;
      Etzfg=Rzt1*Et1fg+Rzt2*Et2fg;
    end
    if isDirichlet_r || isFarField
      rDfg=rD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isDirichlet_w || isFarField
      wDfg=wD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      wDxfg=wDfg(:,1);
      wDyfg=wDfg(:,2);
      if nsd==3
        wDzfg=wDfg(:,3);
      end
    end
    if isDirichlet_g || isFarField
      gDfg=gD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
    end
    if isDirichlet_e
      eDfg=eD(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      if nsd==2
        eDt1fg=Rxt1*eDfg(:,1)+Ryt1*eDfg(:,2);
      elseif nsd==3
        eDt1fg=Rxt1*eDfg(:,1)+Ryt1*eDfg(:,2)+Rzt1*eDfg(:,3);
        eDt2fg=Rxt2*eDfg(:,1)+Ryt2*eDfg(:,2)+Rzt2*eDfg(:,3);
      end
      if nsd==2
        eDtxfg=Rxt1*eDt1fg;
        eDtyfg=Ryt1*eDt1fg;
      elseif nsd==3
        eDtxfg=Rxt1*eDt1fg+Rxt2*eDt2fg;
        eDtyfg=Ryt1*eDt1fg+Ryt2*eDt2fg;
        eDtzfg=Rzt1*eDt1fg+Rzt2*eDt2fg;
      end
    elseif isAbsorbing
      ifg=i(Xfg(:,1),Xfg(:,2),Xfg(:,3),t);
      ixfg=ifg(:,1);
      iyfg=ifg(:,2);
      if nsd==3
        izfg=ifg(:,3);
      end
    end
    if isa(mu,'function_handle')
      mufg=mu(Xfg(:,1),Xfg(:,2),Xfg(:,3));
      epsilonfg=epsilon(Xfg(:,1),Xfg(:,2),Xfg(:,3));
    else
      mufg=mu;
      epsilonfg=epsilon;
    end
    
    % Compute basic matrices
    NwfT=(wfg.*Nf)';
    Mf=NwfT*Nf;
    Mf=(Mf+Mf')/2;
    Mfnx=Mf*nx;
    Mfny=Mf*ny;
    if nsd==3
      Mfnz=Mf*nz;
    end
    if nsd==2
      Mfnxnx=Mf*nx*nx;
      Mfnxny=Mf*nx*ny;
      Mfnyny=Mf*ny*ny;
    else
      Mfnxnx=Mf*nx*nx;
      Mfnxny=Mf*nx*ny;
      Mfnxnz=Mf*nx*nz;
      Mfnyny=Mf*ny*ny;
      Mfnynz=Mf*ny*nz;
      Mfnznz=Mf*nz*nz;
    end
    
    % Compute terms for speedup
    Vxfg=Wxfg./Rfg;
    Vyfg=Wyfg./Rfg;
    if nsd==3
      Vzfg=Wzfg./Rfg;
    end
    if nsd==2
      V2fg=Vxfg.^2+Vyfg.^2;
      Vnfg=Vxfg*nx+Vyfg*ny;
      Wnfg=Wxfg*nx+Wyfg*ny;
    elseif nsd==3
      V2fg=Vxfg.^2+Vyfg.^2+Vzfg.^2;
      Vnfg=Vxfg*nx+Vyfg*ny+Vzfg*nz;
      Wnfg=Wxfg*nx+Wyfg*ny+Wzfg*nz;
    end
    G_Rfg=Gfg./Rfg;
    Pfg=(gamma-1)*(Gfg-1/2*Rfg.*V2fg);
    
    % Lax-Friedrichs flux
    Cfg=abs(sqrt(gamma*Pfg./Rfg));
    lamdbaMax=max(max(abs(Vnfg)+Cfg),1e-3);
    tauR=lamdbaMax;
    tauW=lamdbaMax;
    tauG=lamdbaMax;
    
    % Compute spectral decomposition
    if isFarField
      if nsd==2
        Wfg=[Wxfg,Wyfg];
        n=[nx,ny];
      elseif nsd==3
        Wfg=[Wxfg,Wyfg,Wzfg];
        n=[nx,ny,nz];
      end
      [Anp,Anm]=computeSpectralDecompositionEulerEquations(gamma,Rfg,Wfg,Gfg,n);
    end
    
    % Compute lhs
    KrW(nf1,nefW1)=KrW(nf1,nefW1)+Mfnx;
    KrW(nf1,nefW2)=KrW(nf1,nefW2)+Mfny;
    if nsd==3
      KrW(nf1,nefW3)=KrW(nf1,nefW3)+Mfnz;
    end
    
    Krr(nf1,nf1)=Krr(nf1,nf1)+tauR*Mf;
    
    KrR(nf1,nefR1)=KrR(nf1,nefR1)-tauR*Mf;
    
    KwR(nf1,nefR1)=KwR(nf1,nefR1)+NwfT*((-Vnfg.*Vxfg+(gamma-1)*1/2*V2fg*nx).*Nf);
    KwR(nf2,nefR1)=KwR(nf2,nefR1)+NwfT*((-Vnfg.*Vyfg+(gamma-1)*1/2*V2fg*ny).*Nf);
    if nsd==3
      KwR(nf3,nefR1)=KwR(nf3,nefR1)+NwfT*((-Vnfg.*Vzfg+(gamma-1)*1/2*V2fg*nz).*Nf);
    end
    
    if nsd==2
      KwW(nf1,nefW1)=KwW(nf1,nefW1)+NwfT*((2*Vxfg*nx+Vyfg*ny-(gamma-1)*Vxfg*nx-tauW).*Nf);
      KwW(nf2,nefW1)=KwW(nf2,nefW1)+NwfT*((Vyfg*nx          -(gamma-1)*Vxfg*ny).*Nf);
      KwW(nf1,nefW2)=KwW(nf1,nefW2)+NwfT*((Vxfg*ny          -(gamma-1)*Vyfg*nx).*Nf);
      KwW(nf2,nefW2)=KwW(nf2,nefW2)+NwfT*((Vxfg*nx+2*Vyfg*ny-(gamma-1)*Vyfg*ny-tauW).*Nf);
    elseif nsd==3
      KwW(nf1,nefW1)=KwW(nf1,nefW1)+NwfT*((2*Vxfg*nx+Vyfg*ny+Vzfg*nz-(gamma-1)*Vxfg*nx-tauW).*Nf);
      KwW(nf2,nefW1)=KwW(nf2,nefW1)+NwfT*((Vyfg*nx                  -(gamma-1)*Vxfg*ny).*Nf);
      KwW(nf3,nefW1)=KwW(nf3,nefW1)+NwfT*((Vzfg*nx                  -(gamma-1)*Vxfg*nz).*Nf);
      KwW(nf1,nefW2)=KwW(nf1,nefW2)+NwfT*((Vxfg*ny                  -(gamma-1)*Vyfg*nx).*Nf);
      KwW(nf2,nefW2)=KwW(nf2,nefW2)+NwfT*((Vxfg*nx+2*Vyfg*ny+Vzfg*nz-(gamma-1)*Vyfg*ny-tauW).*Nf);
      KwW(nf3,nefW2)=KwW(nf3,nefW2)+NwfT*((Vzfg*ny                  -(gamma-1)*Vyfg*nz).*Nf);
      KwW(nf1,nefW3)=KwW(nf1,nefW3)+NwfT*((Vxfg*nz                  -(gamma-1)*Vzfg*nx).*Nf);
      KwW(nf2,nefW3)=KwW(nf2,nefW3)+NwfT*((Vyfg*nz                  -(gamma-1)*Vzfg*ny).*Nf);
      KwW(nf3,nefW3)=KwW(nf3,nefW3)+NwfT*((Vxfg*nx+Vyfg*ny+2*Vzfg*nz-(gamma-1)*Vzfg*nz-tauW).*Nf);
    end
    
    KwG(nf1,nefG1)=KwG(nf1,nefG1)+(gamma-1)*Mfnx;
    KwG(nf2,nefG1)=KwG(nf2,nefG1)+(gamma-1)*Mfny;
    if nsd==3
      KwG(nf3,nefG1)=KwG(nf3,nefG1)+(gamma-1)*Mfnz;
    end
    
    Kww(nf1,nf1)=Kww(nf1,nf1)+tauW*Mf;
    Kww(nf2,nf2)=Kww(nf2,nf2)+tauW*Mf;
    if nsd==3
      Kww(nf3,nf3)=Kww(nf3,nf3)+tauW*Mf;
    end
    
    KgR(nf1,nefR1)=KgR(nf1,nefR1)+NwfT*(((-G_Rfg+(gamma-1)*(-G_Rfg+V2fg)).*Vnfg).*Nf);
    
    KgW(nf1,nefW1)=KgW(nf1,nefW1)+NwfT*((G_Rfg*nx+(gamma-1)*((G_Rfg-1/2*V2fg)*nx-Vnfg.*Vxfg)).*Nf);
    KgW(nf1,nefW2)=KgW(nf1,nefW2)+NwfT*((G_Rfg*ny+(gamma-1)*((G_Rfg-1/2*V2fg)*ny-Vnfg.*Vyfg)).*Nf);
    if nsd==3
     KgW(nf1,nefW3)=KgW(nf1,nefW3)+NwfT*((G_Rfg*nz+(gamma-1)*((G_Rfg-1/2*V2fg)*nz-Vnfg.*Vzfg)).*Nf);
    end
    
    KgG(nf1,nefG1)=KgG(nf1,nefG1)+NwfT*((gamma*Vnfg-tauG).*Nf);
    
    Kgg(nf1,nf1)=Kgg(nf1,nf1)+tauG*Mf;
    
    if nsd==2
      KHE(nf1,nefE1)=KHE(nf1,nefE1)+Rxt1*Mfny-Ryt1*Mfnx;
    else
      KHE(nf1,nefE1)=KHE(nf1,nefE1)-Rzt1*Mfny+Ryt1*Mfnz;
      KHE(nf2,nefE1)=KHE(nf2,nefE1)-Rxt1*Mfnz+Rzt1*Mfnx;
      KHE(nf3,nefE1)=KHE(nf3,nefE1)-Ryt1*Mfnx+Rxt1*Mfny;
      KHE(nf1,nefE2)=KHE(nf1,nefE2)-Rzt2*Mfny+Ryt2*Mfnz;
      KHE(nf2,nefE2)=KHE(nf2,nefE2)-Rxt2*Mfnz+Rzt2*Mfnx;
      KHE(nf3,nefE2)=KHE(nf3,nefE2)-Ryt2*Mfnx+Rxt2*Mfny;
    end
    
    if nsd==2
      Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*Mfnyny;
      Kee(nf2,nf1)=Kee(nf2,nf1)-tauE*Mfnxny;
      Kee(nf1,nf2)=Kee(nf1,nf2)-tauE*Mfnxny;
      Kee(nf2,nf2)=Kee(nf2,nf2)+tauE*Mfnxnx;
    else
      Kee(nf1,nf1)=Kee(nf1,nf1)+tauE*(Mfnyny+Mfnznz);
      Kee(nf2,nf1)=Kee(nf2,nf1)-tauE*Mfnxny;
      Kee(nf3,nf1)=Kee(nf3,nf1)-tauE*Mfnxnz;
      Kee(nf1,nf2)=Kee(nf1,nf2)-tauE*Mfnxny;
      Kee(nf2,nf2)=Kee(nf2,nf2)+tauE*(Mfnxnx+Mfnznz);
      Kee(nf3,nf2)=Kee(nf3,nf2)-tauE*Mfnynz;
      Kee(nf1,nf3)=Kee(nf1,nf3)-tauE*Mfnxnz;
      Kee(nf2,nf3)=Kee(nf2,nf3)-tauE*Mfnynz;
      Kee(nf3,nf3)=Kee(nf3,nf3)+tauE*(Mfnxnx+Mfnyny);
    end
    
    if nsd==2
      KeE(nf1,nefE1)=KeE(nf1,nefE1)-Rxt1*tauE*Mf;
      KeE(nf2,nefE1)=KeE(nf2,nefE1)-Ryt1*tauE*Mf;
    else
      KeE(nf1,nefE1)=KeE(nf1,nefE1)-Rxt1*tauE*Mf;
      KeE(nf2,nefE1)=KeE(nf2,nefE1)-Ryt1*tauE*Mf;
      KeE(nf3,nefE1)=KeE(nf3,nefE1)-Rzt1*tauE*Mf;
      KeE(nf1,nefE2)=KeE(nf1,nefE2)-Rxt2*tauE*Mf;
      KeE(nf2,nefE2)=KeE(nf2,nefE2)-Ryt2*tauE*Mf;
      KeE(nf3,nefE2)=KeE(nf3,nefE2)-Rzt2*tauE*Mf;
    end
    
    if isInterior
      KRW(nefR1,nefW1)=-Mfnx;
      KRW(nefR1,nefW2)=-Mfny;
      if nsd==3
        KRW(nefR1,nefW3)=-Mfnz;
      end
      
      KRr(nefR1,nf1)=KRr(nefR1,nf1)-tauR*Mf;
      
      KRR(nefR1,nefR1)=tauR*Mf;
    elseif isDirichlet_r
      KRR(nefR1,nefR1)=Mf;
    elseif isInviscidWall
      KRR(nefR1,nefR1)=Mf;
      
      KRr(nefR1,nf1)=KRr(nefR1,nf1)-Mf;
    elseif isFarField
      if nsd==2
        KRR(nefR1,nefR1)=NwfT*((Anp(:,1,1)+Anm(:,1,1)).*Nf);
        KRW(nefR1,nefW1)=NwfT*((Anp(:,1,2)+Anm(:,1,2)).*Nf);
        KRW(nefR1,nefW2)=NwfT*((Anp(:,1,3)+Anm(:,1,3)).*Nf);
        KRG(nefR1,nefG1)=NwfT*((Anp(:,1,4)+Anm(:,1,4)).*Nf);
        
        KRr(nefR1,nf1)=KRr(nefR1,nf1)-NwfT*((Anp(:,1,1)).*Nf);
        KRw(nefR1,nf1)=KRw(nefR1,nf1)-NwfT*((Anp(:,1,2)).*Nf);
        KRw(nefR1,nf2)=KRw(nefR1,nf2)-NwfT*((Anp(:,1,3)).*Nf);
        KRg(nefR1,nf1)=KRg(nefR1,nf1)-NwfT*((Anp(:,1,4)).*Nf);
      elseif nsd==3
        KRR(nefR1,nefR1)=NwfT*((Anp(:,1,1)+Anm(:,1,1)).*Nf);
        KRW(nefR1,nefW1)=NwfT*((Anp(:,1,2)+Anm(:,1,2)).*Nf);
        KRW(nefR1,nefW2)=NwfT*((Anp(:,1,3)+Anm(:,1,3)).*Nf);
        KRW(nefR1,nefW3)=NwfT*((Anp(:,1,4)+Anm(:,1,4)).*Nf);
        KRG(nefR1,nefG1)=NwfT*((Anp(:,1,5)+Anm(:,1,5)).*Nf);
        
        KRr(nefR1,nf1)=KRr(nefR1,nf1)-NwfT*((Anp(:,1,1)).*Nf);
        KRw(nefR1,nf1)=KRw(nefR1,nf1)-NwfT*((Anp(:,1,2)).*Nf);
        KRw(nefR1,nf2)=KRw(nefR1,nf2)-NwfT*((Anp(:,1,3)).*Nf);
        KRw(nefR1,nf3)=KRw(nefR1,nf3)-NwfT*((Anp(:,1,4)).*Nf);
        KRg(nefR1,nf1)=KRg(nefR1,nf1)-NwfT*((Anp(:,1,5)).*Nf);
      end
    end
    
    if isInterior
      KWR(nefW1,nefR1)=NwfT*((Vnfg.*Vxfg-(gamma-1)*1/2*V2fg*nx).*Nf);
      KWR(nefW2,nefR1)=NwfT*((Vnfg.*Vyfg-(gamma-1)*1/2*V2fg*ny).*Nf);
      if nsd==3
        KWR(nefW3,nefR1)=NwfT*((Vnfg.*Vzfg-(gamma-1)*1/2*V2fg*nz).*Nf);
      end
      
      if nsd==2
        KWW(nefW1,nefW1)=-NwfT*((Vyfg*ny+(3-gamma)*Vxfg*nx-tauW).*Nf);
        KWW(nefW2,nefW1)=-NwfT*((Vyfg*nx+(1-gamma)*Vxfg*ny).*Nf);
        KWW(nefW1,nefW2)=-NwfT*((Vxfg*ny+(1-gamma)*Vyfg*nx).*Nf);
        KWW(nefW2,nefW2)=-NwfT*((Vxfg*nx+(3-gamma)*Vyfg*ny-tauW).*Nf);
      elseif nsd==3
        KWW(nefW1,nefW1)=-NwfT*((Vyfg*ny+Vzfg*nz+(3-gamma)*Vxfg*nx-tauW).*Nf);
        KWW(nefW2,nefW1)=-NwfT*((Vyfg*nx        +(1-gamma)*Vxfg*ny).*Nf);
        KWW(nefW3,nefW1)=-NwfT*((Vzfg*nx        +(1-gamma)*Vxfg*nz).*Nf);
        KWW(nefW1,nefW2)=-NwfT*((Vxfg*ny        +(1-gamma)*Vyfg*nx).*Nf);
        KWW(nefW2,nefW2)=-NwfT*((Vxfg*nx+Vzfg*nz+(3-gamma)*Vyfg*ny-tauW).*Nf);
        KWW(nefW3,nefW2)=-NwfT*((Vzfg*ny        +(1-gamma)*Vyfg*nz).*Nf);
        KWW(nefW1,nefW3)=-NwfT*((Vxfg*nz        +(1-gamma)*Vzfg*nx).*Nf);
        KWW(nefW2,nefW3)=-NwfT*((Vyfg*nz        +(1-gamma)*Vzfg*ny).*Nf);
        KWW(nefW3,nefW3)=-NwfT*((Vxfg*nx+Vyfg*ny+(3-gamma)*Vzfg*nz-tauW).*Nf);
      end
      
      KWG(nefW1,nefG1)=-(gamma-1)*Mfnx;
      KWG(nefW2,nefG1)=-(gamma-1)*Mfny;
      if nsd==3
        KWG(nefW3,nefG1)=-(gamma-1)*Mfnz;
      end
      
      KWw(nefW1,nf1)=KWw(nefW1,nf1)-tauW*Mf;
      KWw(nefW2,nf2)=KWw(nefW2,nf2)-tauW*Mf;
      if nsd==3
        KWw(nefW3,nf3)=KWw(nefW3,nf3)-tauW*Mf;
      end
    elseif isDirichlet_w
      KWW(nefW1,nefW1)=Mf;
      KWW(nefW2,nefW2)=Mf;
      if nsd==3
        KWW(nefW3,nefW3)=Mf;
      end
    elseif isInviscidWall
      KWW(nefW1,nefW1)=Mf;
      KWW(nefW2,nefW2)=Mf;
      if nsd==3
        KWW(nefW3,nefW3)=Mf;
      end
      
      if nsd==2
        KWw(nefW1,nf1)=KWw(nefW1,nf1)+Mfnxnx-Mf;
        KWw(nefW2,nf1)=KWw(nefW2,nf1)+Mfnxny;
        KWw(nefW1,nf2)=KWw(nefW1,nf2)+Mfnxny;
        KWw(nefW2,nf2)=KWw(nefW2,nf2)+Mfnyny-Mf;
      elseif nsd==3
        KWw(nefW1,nf1)=KWw(nefW1,nf1)+Mfnxnx-Mf;
        KWw(nefW2,nf1)=KWw(nefW2,nf1)+Mfnxny;
        KWw(nefW3,nf1)=KWw(nefW3,nf1)+Mfnxnz;
        KWw(nefW1,nf2)=KWw(nefW1,nf2)+Mfnxny;
        KWw(nefW2,nf2)=KWw(nefW2,nf2)+Mfnyny-Mf;
        KWw(nefW3,nf2)=KWw(nefW3,nf2)+Mfnynz;
        KWw(nefW1,nf3)=KWw(nefW1,nf3)+Mfnxnz;
        KWw(nefW2,nf3)=KWw(nefW2,nf3)+Mfnynz;
        KWw(nefW3,nf3)=KWw(nefW3,nf3)+Mfnznz-Mf;
      end
    elseif isFarField
      if nsd==2
        KWR(nefW1,nefR1)=NwfT*((Anp(:,2,1)+Anm(:,2,1)).*Nf);
        KWW(nefW1,nefW1)=NwfT*((Anp(:,2,2)+Anm(:,2,2)).*Nf);
        KWW(nefW1,nefW2)=NwfT*((Anp(:,2,3)+Anm(:,2,3)).*Nf);
        KWG(nefW1,nefG1)=NwfT*((Anp(:,2,4)+Anm(:,2,4)).*Nf);
        
        KWr(nefW1,nf1)=KWr(nefW1,nf1)-NwfT*((Anp(:,2,1)).*Nf);
        KWw(nefW1,nf1)=KWw(nefW1,nf1)-NwfT*((Anp(:,2,2)).*Nf);
        KWw(nefW1,nf2)=KWw(nefW1,nf2)-NwfT*((Anp(:,2,3)).*Nf);
        KWg(nefW1,nf1)=KWg(nefW1,nf1)-NwfT*((Anp(:,2,4)).*Nf);
        
        KWR(nefW2,nefR1)=NwfT*((Anp(:,3,1)+Anm(:,3,1)).*Nf);
        KWW(nefW2,nefW1)=NwfT*((Anp(:,3,2)+Anm(:,3,2)).*Nf);
        KWW(nefW2,nefW2)=NwfT*((Anp(:,3,3)+Anm(:,3,3)).*Nf);
        KWG(nefW2,nefG1)=NwfT*((Anp(:,3,4)+Anm(:,3,4)).*Nf);
        
        KWr(nefW2,nf1)=KWr(nefW2,nf1)-NwfT*((Anp(:,3,1)).*Nf);
        KWw(nefW2,nf1)=KWw(nefW2,nf1)-NwfT*((Anp(:,3,2)).*Nf);
        KWw(nefW2,nf2)=KWw(nefW2,nf2)-NwfT*((Anp(:,3,3)).*Nf);
        KWg(nefW2,nf1)=KWg(nefW2,nf1)-NwfT*((Anp(:,3,4)).*Nf);
      elseif nsd==3
        KWR(nefW1,nefR1)=NwfT*((Anp(:,2,1)+Anm(:,2,1)).*Nf);
        KWW(nefW1,nefW1)=NwfT*((Anp(:,2,2)+Anm(:,2,2)).*Nf);
        KWW(nefW1,nefW2)=NwfT*((Anp(:,2,3)+Anm(:,2,3)).*Nf);
        KWW(nefW1,nefW3)=NwfT*((Anp(:,2,4)+Anm(:,2,4)).*Nf);
        KWG(nefW1,nefG1)=NwfT*((Anp(:,2,5)+Anm(:,2,5)).*Nf);
        
        KWr(nefW1,nf1)=KWr(nefW1,nf1)-NwfT*((Anp(:,2,1)).*Nf);
        KWw(nefW1,nf1)=KWw(nefW1,nf1)-NwfT*((Anp(:,2,2)).*Nf);
        KWw(nefW1,nf2)=KWw(nefW1,nf2)-NwfT*((Anp(:,2,3)).*Nf);
        KWw(nefW1,nf3)=KWw(nefW1,nf3)-NwfT*((Anp(:,2,4)).*Nf);
        KWg(nefW1,nf1)=KWg(nefW1,nf1)-NwfT*((Anp(:,2,5)).*Nf);
        
        KWR(nefW2,nefR1)=NwfT*((Anp(:,3,1)+Anm(:,3,1)).*Nf);
        KWW(nefW2,nefW1)=NwfT*((Anp(:,3,2)+Anm(:,3,2)).*Nf);
        KWW(nefW2,nefW2)=NwfT*((Anp(:,3,3)+Anm(:,3,3)).*Nf);
        KWW(nefW2,nefW3)=NwfT*((Anp(:,3,4)+Anm(:,3,4)).*Nf);
        KWG(nefW2,nefG1)=NwfT*((Anp(:,3,5)+Anm(:,3,5)).*Nf);
        
        KWr(nefW2,nf1)=KWr(nefW2,nf1)-NwfT*((Anp(:,3,1)).*Nf);
        KWw(nefW2,nf1)=KWw(nefW2,nf1)-NwfT*((Anp(:,3,2)).*Nf);
        KWw(nefW2,nf2)=KWw(nefW2,nf2)-NwfT*((Anp(:,3,3)).*Nf);
        KWw(nefW2,nf3)=KWw(nefW2,nf3)-NwfT*((Anp(:,3,4)).*Nf);
        KWg(nefW2,nf1)=KWg(nefW2,nf1)-NwfT*((Anp(:,3,5)).*Nf);
        
        KWR(nefW3,nefR1)=NwfT*((Anp(:,4,1)+Anm(:,4,1)).*Nf);
        KWW(nefW3,nefW1)=NwfT*((Anp(:,4,2)+Anm(:,4,2)).*Nf);
        KWW(nefW3,nefW2)=NwfT*((Anp(:,4,3)+Anm(:,4,3)).*Nf);
        KWW(nefW3,nefW3)=NwfT*((Anp(:,4,4)+Anm(:,4,4)).*Nf);
        KWG(nefW3,nefG1)=NwfT*((Anp(:,4,5)+Anm(:,4,5)).*Nf);
        
        KWr(nefW3,nf1)=KWr(nefW3,nf1)-NwfT*((Anp(:,4,1)).*Nf);
        KWw(nefW3,nf1)=KWw(nefW3,nf1)-NwfT*((Anp(:,4,2)).*Nf);
        KWw(nefW3,nf2)=KWw(nefW3,nf2)-NwfT*((Anp(:,4,3)).*Nf);
        KWw(nefW3,nf3)=KWw(nefW3,nf3)-NwfT*((Anp(:,4,4)).*Nf);
        KWg(nefW3,nf1)=KWg(nefW3,nf1)-NwfT*((Anp(:,4,5)).*Nf);
      end
    end
    
    if isInterior
      KGR(nefG1,nefR1)=-NwfT*(((-G_Rfg+(gamma-1)*(-G_Rfg+V2fg)).*Vnfg).*Nf);
      
      KGW(nefG1,nefW1)=-NwfT*((G_Rfg*nx+(gamma-1)*((G_Rfg-1/2*V2fg)*nx-Vnfg.*Vxfg)).*Nf);
      KGW(nefG1,nefW2)=-NwfT*((G_Rfg*ny+(gamma-1)*((G_Rfg-1/2*V2fg)*ny-Vnfg.*Vyfg)).*Nf);
      if nsd==3
        KGW(nefG1,nefW3)=-NwfT*((G_Rfg*nz+(gamma-1)*((G_Rfg-1/2*V2fg)*nz-Vnfg.*Vzfg)).*Nf);
      end
      
      KGG(nefG1,nefG1)=-NwfT*((gamma*Vnfg-tauG).*Nf);
      
      KGg(nefG1,nf1)=KGg(nefG1,nf1)-tauG*Mf;
    elseif isDirichlet_g
      KGG(nefG1,nefG1)=Mf;
    elseif isInviscidWall
      KGG(nefG1,nefG1)=Mf;
      
      KGg(nefG1,nf1)=KGg(nefG1,nf1)-Mf;
    elseif isFarField
      if nsd==2
        KGR(nefG1,nefR1)=NwfT*((Anp(:,4,1)+Anm(:,4,1)).*Nf);
        KGW(nefG1,nefW1)=NwfT*((Anp(:,4,2)+Anm(:,4,2)).*Nf);
        KGW(nefG1,nefW2)=NwfT*((Anp(:,4,3)+Anm(:,4,3)).*Nf);
        KGG(nefG1,nefG1)=NwfT*((Anp(:,4,4)+Anm(:,4,4)).*Nf);
        
        KGr(nefG1,nf1)=KGr(nefG1,nf1)-NwfT*((Anp(:,4,1)).*Nf);
        KGw(nefG1,nf1)=KGw(nefG1,nf1)-NwfT*((Anp(:,4,2)).*Nf);
        KGw(nefG1,nf2)=KGw(nefG1,nf2)-NwfT*((Anp(:,4,3)).*Nf);
        KGg(nefG1,nf1)=KGg(nefG1,nf1)-NwfT*((Anp(:,4,4)).*Nf);
      elseif nsd==3
        KGR(nefG1,nefR1)=NwfT*((Anp(:,5,1)+Anm(:,5,1)).*Nf);
        KGW(nefG1,nefW1)=NwfT*((Anp(:,5,2)+Anm(:,5,2)).*Nf);
        KGW(nefG1,nefW2)=NwfT*((Anp(:,5,3)+Anm(:,5,3)).*Nf);
        KGW(nefG1,nefW3)=NwfT*((Anp(:,5,4)+Anm(:,5,4)).*Nf);
        KGG(nefG1,nefG1)=NwfT*((Anp(:,5,5)+Anm(:,5,5)).*Nf);
        
        KGr(nefG1,nf1)=KGr(nefG1,nf1)-NwfT*((Anp(:,5,1)).*Nf);
        KGw(nefG1,nf1)=KGw(nefG1,nf1)-NwfT*((Anp(:,5,2)).*Nf);
        KGw(nefG1,nf2)=KGw(nefG1,nf2)-NwfT*((Anp(:,5,3)).*Nf);
        KGw(nefG1,nf3)=KGw(nefG1,nf3)-NwfT*((Anp(:,5,4)).*Nf);
        KGg(nefG1,nf1)=KGg(nefG1,nf1)-NwfT*((Anp(:,5,5)).*Nf);
      end
    end
    
    if isInterior
      if nsd==2
        KEH(nefE1,nf1)=KEH(nefE1,nf1)+Rxt1*Mfny-Ryt1*Mfnx;
      else
        KEH(nefE1,nf1)=KEH(nefE1,nf1)-Rzt1*Mfny+Ryt1*Mfnz;
        KEH(nefE2,nf1)=KEH(nefE2,nf1)-Rzt2*Mfny+Ryt2*Mfnz;
        KEH(nefE1,nf2)=KEH(nefE1,nf2)-Rxt1*Mfnz+Rzt1*Mfnx;
        KEH(nefE2,nf2)=KEH(nefE2,nf2)-Rxt2*Mfnz+Rzt2*Mfnx;
        KEH(nefE1,nf3)=KEH(nefE1,nf3)-Ryt1*Mfnx+Rxt1*Mfny;
        KEH(nefE2,nf3)=KEH(nefE2,nf3)-Ryt2*Mfnx+Rxt2*Mfny;
      end
      
      if nsd==2
        KEe(nefE1,nf1)=KEe(nefE1,nf1)-Rxt1*tauE*Mf;
        KEe(nefE1,nf2)=KEe(nefE1,nf2)-Ryt1*tauE*Mf;
      else
        KEe(nefE1,nf1)=KEe(nefE1,nf1)-Rxt1*tauE*Mf;
        KEe(nefE2,nf1)=KEe(nefE2,nf1)-Rxt2*tauE*Mf;
        KEe(nefE1,nf2)=KEe(nefE1,nf2)-Ryt1*tauE*Mf;
        KEe(nefE2,nf2)=KEe(nefE2,nf2)-Ryt2*tauE*Mf;
        KEe(nefE1,nf3)=KEe(nefE1,nf3)-Rzt1*tauE*Mf;
        KEe(nefE2,nf3)=KEe(nefE2,nf3)-Rzt2*tauE*Mf;
      end
      
      if nsd==2
        KEE(nefE1,nefE1)=(Rxt1*Rxt1+Ryt1*Ryt1)*tauE*Mf;
      else
        KEE(nefE1,nefE1)=(Rxt1*Rxt1+Ryt1*Ryt1+Rzt1*Rzt1)*tauE*Mf;
        KEE(nefE2,nefE1)=(Rxt1*Rxt2+Ryt1*Ryt2+Rzt1*Rzt2)*tauE*Mf;
        KEE(nefE1,nefE2)=(Rxt1*Rxt2+Ryt1*Ryt2+Rzt1*Rzt2)*tauE*Mf;
        KEE(nefE2,nefE2)=(Rxt2*Rxt2+Ryt2*Ryt2+Rzt2*Rzt2)*tauE*Mf;
      end
    elseif isDirichlet_e
      if nsd==2
        KEE(nefE1,nefE1)=(Rxt1*Rxt1+Ryt1*Ryt1)*Mf;
      else
        KEE(nefE1,nefE1)=(Rxt1*Rxt1+Ryt1*Ryt1+Rzt1*Rzt1)*Mf;
        KEE(nefE2,nefE1)=(Rxt1*Rxt2+Ryt1*Ryt2+Rzt1*Rzt2)*Mf;
        KEE(nefE1,nefE2)=(Rxt1*Rxt2+Ryt1*Ryt2+Rzt1*Rzt2)*Mf;
        KEE(nefE2,nefE2)=(Rxt2*Rxt2+Ryt2*Ryt2+Rzt2*Rzt2)*Mf;
      end
    elseif isAbsorbing
      if nsd==2
        KEH(nefE1,nf1)=KEH(nefE1,nf1)+Rxt1*Mfny-Ryt1*Mfnx;
      else
        KEH(nefE1,nf1)=KEH(nefE1,nf1)-Rzt1*Mfny+Ryt1*Mfnz;
        KEH(nefE2,nf1)=KEH(nefE2,nf1)-Rzt2*Mfny+Ryt2*Mfnz;
        KEH(nefE1,nf2)=KEH(nefE1,nf2)-Rxt1*Mfnz+Rzt1*Mfnx;
        KEH(nefE2,nf2)=KEH(nefE2,nf2)-Rxt2*Mfnz+Rzt2*Mfnx;
        KEH(nefE1,nf3)=KEH(nefE1,nf3)-Ryt1*Mfnx+Rxt1*Mfny;
        KEH(nefE2,nf3)=KEH(nefE2,nf3)-Ryt2*Mfnx+Rxt2*Mfny;
      end
      
      if nsd==2
        KEE(nefE1,nefE1)=NwfT*((sqrt(epsilonfg./mufg)*(Rxt1*Rxt1+Ryt1*Ryt1)).*Nf);
      else
        KEE(nefE1,nefE1)=NwfT*((sqrt(epsilonfg./mufg)*(Rxt1*Rxt1+Ryt1*Ryt1+Rzt1*Rzt1)).*Nf);
        KEE(nefE2,nefE1)=NwfT*((sqrt(epsilonfg./mufg)*(Rxt1*Rxt2+Ryt1*Ryt2+Rzt1*Rzt2)).*Nf);
        KEE(nefE1,nefE2)=NwfT*((sqrt(epsilonfg./mufg)*(Rxt1*Rxt2+Ryt1*Ryt2+Rzt1*Rzt2)).*Nf);
        KEE(nefE2,nefE2)=NwfT*((sqrt(epsilonfg./mufg)*(Rxt2*Rxt2+Ryt2*Ryt2+Rzt2*Rzt2)).*Nf);
      end
    end
    
    % Compute rhs
    fr(nf1,1)=fr(nf1,1)-NwfT*(Wnfg+tauR*(rfg-Rfg));
    
    fw(nf1,1)=fw(nf1,1)-NwfT*(Wnfg.*Vxfg+Pfg*nx+tauW*(wxfg-Wxfg));
    fw(nf2,1)=fw(nf2,1)-NwfT*(Wnfg.*Vyfg+Pfg*ny+tauW*(wyfg-Wyfg));
    if nsd==3
      fw(nf3,1)=fw(nf3,1)-NwfT*(Wnfg.*Vzfg+Pfg*nz+tauW*(wzfg-Wzfg));
    end
    
    fg(nf1,1)=fg(nf1,1)-NwfT*((Gfg+Pfg).*Vnfg+tauG*(gfg-Gfg));
    
    if nsd==2
      fH(nf1,1)=fH(nf1,1)+NwfT*(nx*Etyfg-ny*Etxfg);
    elseif nsd==3
      fH(nf1,1)=fH(nf1,1)+NwfT*(ny*Etzfg-nz*Etyfg);
      fH(nf2,1)=fH(nf2,1)+NwfT*(nz*Etxfg-nx*Etzfg);
      fH(nf3,1)=fH(nf3,1)+NwfT*(nx*Etyfg-ny*Etxfg);
    end
    
    if nsd==2
      fe(nf1,1)=fe(nf1,1)-NwfT*(tauE*((-ny*(nx*eyfg-ny*exfg))-Etxfg));
      fe(nf2,1)=fe(nf2,1)-NwfT*(tauE*((+nx*(nx*eyfg-ny*exfg))-Etyfg));
    elseif nsd==3
      fe(nf1,1)=fe(nf1,1)-NwfT*(tauE*((-ny*(nx*eyfg-ny*exfg)+nz*(nz*exfg-nx*ezfg))-Etxfg));
      fe(nf2,1)=fe(nf2,1)-NwfT*(tauE*((-nz*(ny*ezfg-nz*eyfg)+nx*(nx*eyfg-ny*exfg))-Etyfg));
      fe(nf3,1)=fe(nf3,1)-NwfT*(tauE*((-nx*(nz*exfg-nx*ezfg)+ny*(ny*ezfg-nz*eyfg))-Etzfg));
    end
    
    if isInterior
      fR(nefR1,1)=NwfT*(Wnfg+tauR*(rfg-Rfg));
    elseif isDirichlet_r
      fR(nefR1,1)=-NwfT*(Rfg-rDfg);
    elseif isInviscidWall
      fR(nefR1,1)=-NwfT*(Rfg-rfg);
    elseif isFarField
      if nsd==2
        fR(nefR1,1)=-NwfT*(Anp(:,1,1).*(Rfg -rfg )+Anm(:,1,1).*(Rfg -rDfg )+...
                           Anp(:,1,2).*(Wxfg-wxfg)+Anm(:,1,2).*(Wxfg-wDxfg)+...
                           Anp(:,1,3).*(Wyfg-wyfg)+Anm(:,1,3).*(Wyfg-wDyfg)+...
                           Anp(:,1,4).*(Gfg -gfg )+Anm(:,1,4).*(Gfg -gDfg ));
      elseif nsd==3
        fR(nefR1,1)=-NwfT*(Anp(:,1,1).*(Rfg -rfg )+Anm(:,1,1).*(Rfg -rDfg )+...
                           Anp(:,1,2).*(Wxfg-wxfg)+Anm(:,1,2).*(Wxfg-wDxfg)+...
                           Anp(:,1,3).*(Wyfg-wyfg)+Anm(:,1,3).*(Wyfg-wDyfg)+...
                           Anp(:,1,4).*(Wzfg-wzfg)+Anm(:,1,4).*(Wzfg-wDzfg)+...
                           Anp(:,1,5).*(Gfg -gfg )+Anm(:,1,5).*(Gfg -gDfg ));
      end
    end
    
    if isInterior
      fW(nefW1,1)=NwfT*(Wnfg.*Vxfg+Pfg*nx+tauW*(wxfg-Wxfg));
      fW(nefW2,1)=NwfT*(Wnfg.*Vyfg+Pfg*ny+tauW*(wyfg-Wyfg));
      if nsd==3
        fW(nefW3,1)=NwfT*(Wnfg.*Vzfg+Pfg*nz+tauW*(wzfg-Wzfg));
      end
    elseif isDirichlet_w
      fW(nefW1,1)=-NwfT*(Wxfg-wDxfg);
      fW(nefW2,1)=-NwfT*(Wyfg-wDyfg);
      if nsd==3
        fW(nefW3,1)=-NwfT*(Wzfg-wDzfg);
      end
    elseif isInviscidWall
      if nsd==2
        fW(nefW1,1)=-NwfT*(Wxfg-(+(1-nx*nx)*wxfg-nx*ny*wyfg));
        fW(nefW2,1)=-NwfT*(Wyfg-(-ny*nx*wxfg+(1-ny*ny)*wyfg));
      elseif nsd==3
        fW(nefW1,1)=-NwfT*(Wxfg-(+(1-nx*nx)*wxfg-nx*ny*wyfg-nx*nz*wzfg));
        fW(nefW2,1)=-NwfT*(Wyfg-(-ny*nx*wxfg+(1-ny*ny)*wyfg-ny*nz*wzfg));
        fW(nefW3,1)=-NwfT*(Wzfg-(-nz*nx*wxfg-nz*ny*wyfg+(1-nz*nz)*wzfg));
      end
    elseif isFarField
      if nsd==2
        fW(nefW1,1)=-NwfT*(Anp(:,2,1).*(Rfg -rfg )+Anm(:,2,1).*(Rfg -rDfg )+...
                           Anp(:,2,2).*(Wxfg-wxfg)+Anm(:,2,2).*(Wxfg-wDxfg)+...
                           Anp(:,2,3).*(Wyfg-wyfg)+Anm(:,2,3).*(Wyfg-wDyfg)+...
                           Anp(:,2,4).*(Gfg -gfg )+Anm(:,2,4).*(Gfg -gDfg ));
        fW(nefW2,1)=-NwfT*(Anp(:,3,1).*(Rfg -rfg )+Anm(:,3,1).*(Rfg -rDfg )+...
                           Anp(:,3,2).*(Wxfg-wxfg)+Anm(:,3,2).*(Wxfg-wDxfg)+...
                           Anp(:,3,3).*(Wyfg-wyfg)+Anm(:,3,3).*(Wyfg-wDyfg)+...
                           Anp(:,3,4).*(Gfg -gfg )+Anm(:,3,4).*(Gfg -gDfg ));
      elseif nsd==3
        fW(nefW1,1)=-NwfT*(Anp(:,2,1).*(Rfg -rfg )+Anm(:,2,1).*(Rfg -rDfg )+...
                           Anp(:,2,2).*(Wxfg-wxfg)+Anm(:,2,2).*(Wxfg-wDxfg)+...
                           Anp(:,2,3).*(Wyfg-wyfg)+Anm(:,2,3).*(Wyfg-wDyfg)+...
                           Anp(:,2,4).*(Wzfg-wzfg)+Anm(:,2,4).*(Wzfg-wDzfg)+...
                           Anp(:,2,5).*(Gfg -gfg )+Anm(:,2,5).*(Gfg -gDfg ));
        fW(nefW2,1)=-NwfT*(Anp(:,3,1).*(Rfg -rfg )+Anm(:,3,1).*(Rfg -rDfg )+...
                           Anp(:,3,2).*(Wxfg-wxfg)+Anm(:,3,2).*(Wxfg-wDxfg)+...
                           Anp(:,3,3).*(Wyfg-wyfg)+Anm(:,3,3).*(Wyfg-wDyfg)+...
                           Anp(:,3,4).*(Wzfg-wzfg)+Anm(:,3,4).*(Wzfg-wDzfg)+...
                           Anp(:,3,5).*(Gfg -gfg )+Anm(:,3,5).*(Gfg -gDfg ));
        fW(nefW3,1)=-NwfT*(Anp(:,4,1).*(Rfg -rfg )+Anm(:,4,1).*(Rfg -rDfg )+...
                           Anp(:,4,2).*(Wxfg-wxfg)+Anm(:,4,2).*(Wxfg-wDxfg)+...
                           Anp(:,4,3).*(Wyfg-wyfg)+Anm(:,4,3).*(Wyfg-wDyfg)+...
                           Anp(:,4,4).*(Wzfg-wzfg)+Anm(:,4,4).*(Wzfg-wDzfg)+...
                           Anp(:,4,5).*(Gfg -gfg )+Anm(:,4,5).*(Gfg -gDfg ));
      end
    end
    
    if isInterior
      fG(nefG1,1)=NwfT*((Gfg+Pfg).*Vnfg+tauG*(gfg-Gfg));
    elseif isDirichlet_g
      fG(nefG1,1)=-NwfT*(Gfg-gDfg);
    elseif isInviscidWall
      fG(nefG1,1)=-NwfT*(Gfg-gfg);
    elseif isFarField
      if nsd==2
        fG(nefG1,1)=-NwfT*(Anp(:,4,1).*(Rfg -rfg )+Anm(:,4,1).*(Rfg -rDfg )+...
                           Anp(:,4,2).*(Wxfg-wxfg)+Anm(:,4,2).*(Wxfg-wDxfg)+...
                           Anp(:,4,3).*(Wyfg-wyfg)+Anm(:,4,3).*(Wyfg-wDyfg)+...
                           Anp(:,4,4).*(Gfg -gfg )+Anm(:,4,4).*(Gfg -gDfg ));
      elseif nsd==3
        fG(nefG1,1)=-NwfT*(Anp(:,5,1).*(Rfg -rfg )+Anm(:,5,1).*(Rfg -rDfg )+...
                           Anp(:,5,2).*(Wxfg-wxfg)+Anm(:,5,2).*(Wxfg-wDxfg)+...
                           Anp(:,5,3).*(Wyfg-wyfg)+Anm(:,5,3).*(Wyfg-wDyfg)+...
                           Anp(:,5,4).*(Wzfg-wzfg)+Anm(:,5,4).*(Wzfg-wDzfg)+...
                           Anp(:,5,5).*(Gfg -gfg )+Anm(:,5,5).*(Gfg -gDfg ));
      end
    end
    
    if isInterior
      if nsd==2
        fE(nefE1,1)=NwfT*(Rxt1*(-ny*Hzfg+tauE*(exfg-Etxfg))+...
                          Ryt1*(+nx*Hzfg+tauE*(eyfg-Etyfg)));
      elseif nsd==3
        fE(nefE1,1)=NwfT*(Rxt1*(nz*Hyfg-ny*Hzfg+tauE*(exfg-Etxfg))+...
                          Ryt1*(nx*Hzfg-nz*Hxfg+tauE*(eyfg-Etyfg))+...
                          Rzt1*(ny*Hxfg-nx*Hyfg+tauE*(ezfg-Etzfg)));
        fE(nefE2,1)=NwfT*(Rxt2*(nz*Hyfg-ny*Hzfg+tauE*(exfg-Etxfg))+...
                          Ryt2*(nx*Hzfg-nz*Hxfg+tauE*(eyfg-Etyfg))+...
                          Rzt2*(ny*Hxfg-nx*Hyfg+tauE*(ezfg-Etzfg)));
      end
    elseif isDirichlet_e
      if nsd==2
        fE(nefE1,1)=-NwfT*(Rxt1*(Etxfg-eDtxfg)+...
                           Ryt1*(Etyfg-eDtyfg));
      elseif nsd==3
        fE(nefE1,1)=-NwfT*(Rxt1*(Etxfg-eDtxfg)+...
                           Ryt1*(Etyfg-eDtyfg)+...
                           Rzt1*(Etzfg-eDtzfg));
        fE(nefE2,1)=-NwfT*(Rxt2*(Etxfg-eDtxfg)+...
                           Ryt2*(Etyfg-eDtyfg)+...
                           Rzt2*(Etzfg-eDtzfg));
      end
    elseif isAbsorbing
      if nsd==2
        fE(nefE1,1)=NwfT*(Rxt1*(-ny*Hzfg-sqrt(epsilonfg./mufg).*Etxfg-ixfg)+...
                          Ryt1*(+nx*Hzfg-sqrt(epsilonfg./mufg).*Etyfg-iyfg));
      elseif nsd==3
        fE(nefE1,1)=NwfT*(Rxt1*(nz*Hyfg-ny*Hzfg-sqrt(epsilonfg./mufg).*Etxfg-ixfg)+...
                          Ryt1*(nx*Hzfg-nz*Hxfg-sqrt(epsilonfg./mufg).*Etyfg-iyfg)+...
                          Rzt1*(ny*Hxfg-nx*Hyfg-sqrt(epsilonfg./mufg).*Etzfg-izfg));
        fE(nefE2,1)=NwfT*(Rxt2*(nz*Hyfg-ny*Hzfg-sqrt(epsilonfg./mufg).*Etxfg-ixfg)+...
                          Ryt2*(nx*Hzfg-nz*Hxfg-sqrt(epsilonfg./mufg).*Etyfg-iyfg)+...
                          Rzt2*(ny*Hxfg-nx*Hyfg-sqrt(epsilonfg./mufg).*Etzfg-izfg));
      end
    end
  end
end

% Compute elemental contributions to lhs and rhs

% Indices
ir=1:NumElementNodes;
iw=ir(end)+(1:nsd*NumElementNodes);
ig=iw(end)+(1:NumElementNodes);
iH=ig(end)+(1:qsd*NumElementNodes);
ie=iH(end)+(1:nsd*NumElementNodes);
iM=ie(end)+(1:qsd*NumElementNodes);
iJ=iM(end)+(1:nsd*NumElementNodes);
iR=reshape((0:NumElementFaces-1)*(1+nsd+1+nsd-1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[]);
iW=reshape((0:NumElementFaces-1)*(1+nsd+1+nsd-1)*NumFaceNodes+repmat((1:nsd*NumFaceNodes)',...
  1,NumElementFaces),1,[])+NumFaceNodes;
iG=reshape((0:NumElementFaces-1)*(1+nsd+1+nsd-1)*NumFaceNodes+repmat((1:NumFaceNodes)',...
  1,NumElementFaces),1,[])+(1+nsd)*NumFaceNodes;
iE=reshape((0:NumElementFaces-1)*(1+nsd+1+nsd-1)*NumFaceNodes+repmat((1:(nsd-1)*NumFaceNodes)',...
  1,NumElementFaces),1,[])+(1+nsd+1)*NumFaceNodes;

% Initialization of lhs and rhs
LhsLL=zeros((1+nsd+1+qsd+nsd+(qsd+nsd)*isPML)*NumElementNodes,...
            (1+nsd+1+qsd+nsd+(qsd+nsd)*isPML)*NumElementNodes);
LhsLG=zeros((1+nsd+1+qsd+nsd+(qsd+nsd)*isPML)*NumElementNodes,...
            (1+nsd+1+nsd-1)*NumElementFaces*NumFaceNodes);
LhsGL=zeros((1+nsd+1+nsd-1)*NumElementFaces*NumFaceNodes,...
            (1+nsd+1+qsd+nsd+(qsd+nsd)*isPML)*NumElementNodes);
LhsGG=zeros((1+nsd+1+nsd-1)*NumElementFaces*NumFaceNodes,...
             (1+nsd+1+nsd-1)*NumElementFaces*NumFaceNodes);
RhsL=zeros((1+nsd+1+qsd+nsd+(qsd+nsd)*isPML)*NumElementNodes,1);
RhsG=zeros((1+nsd+1+nsd-1)*NumElementFaces*NumFaceNodes,1);

% Lhs local-local
LhsLL(ir,ir)=Krr;
LhsLL(ir,iw)=Krw;
LhsLL(iw,ir)=Kwr;
LhsLL(iw,iw)=Kww;
LhsLL(iw,ig)=Kwg;
LhsLL(iw,iH)=KwH;
LhsLL(iw,ie)=Kwe;
LhsLL(ig,ir)=Kgr;
LhsLL(ig,iw)=Kgw;
LhsLL(ig,ig)=Kgg;
LhsLL(ig,ie)=Kge;
LhsLL(iH,iH)=KHH;
LhsLL(iH,ie)=KHe;
LhsLL(ie,iw)=Kew;
LhsLL(ie,iH)=KHe';
LhsLL(ie,ie)=Kee;
if isPML
  LhsLL(iH,iM)=KHM;
  LhsLL(ie,iJ)=KeJ;
  LhsLL(iM,iH)=KMH;
  LhsLL(iM,iM)=KMM;
  LhsLL(iJ,ie)=KJe;
  LhsLL(iJ,iJ)=KJJ;
end

% Lhs local-global
LhsLG(ir,iR)=KrR;
LhsLG(ir,iW)=KrW;
LhsLG(iw,iR)=KwR;
LhsLG(iw,iW)=KwW;
LhsLG(iw,iG)=KwG;
LhsLG(ig,iR)=KgR;
LhsLG(ig,iW)=KgW;
LhsLG(ig,iG)=KgG;
LhsLG(iH,iE)=KHE;
LhsLG(ie,iE)=KeE;

% Rhs local
RhsL(ir,1)=fr;
RhsL(iw,1)=fw;
RhsL(ig,1)=fg;
RhsL(iH,1)=fH;
RhsL(ie,1)=fe;
if isPML
  RhsL(iM,1)=fM;
  RhsL(iJ,1)=fJ;
end

% Lhs global-local
LhsGL(iR,ir)=KRr;
LhsGL(iR,iw)=KRw;
LhsGL(iR,ig)=KRg;
LhsGL(iW,ir)=KWr;
LhsGL(iW,iw)=KWw;
LhsGL(iW,ig)=KWg;
LhsGL(iG,ir)=KGr;
LhsGL(iG,iw)=KGw;
LhsGL(iG,ig)=KGg;
LhsGL(iE,iH)=KEH;
LhsGL(iE,ie)=KEe;

% Lhs global-global
LhsGG(iR,iR)=KRR;
LhsGG(iR,iW)=KRW;
LhsGG(iR,iG)=KRG;
LhsGG(iW,iR)=KWR;
LhsGG(iW,iW)=KWW;
LhsGG(iW,iG)=KWG;
LhsGG(iG,iR)=KGR;
LhsGG(iG,iW)=KGW;
LhsGG(iG,iG)=KGG;
LhsGG(iE,iE)=KEE;

% Rhs global
RhsG(iR,1)=fR;
RhsG(iW,1)=fW;
RhsG(iG,1)=fG;
RhsG(iE,1)=fE;

% Add missing terms for axisymmetric formulation
if isAxisymmetric
  yeg=Ne*Xe(:,2);
  
  Krw_axisym=zeros(NumElementNodes,nsd*NumElementNodes);
  Krw_axisym(ne1,ne2)=+NweT*((1./yeg).*Ne);
  LhsLL(ir,iw)=LhsLL(ir,iw)+Krw_axisym;
  
  fr_axisym=zeros(NumElementNodes,1);
  fr_axisym(ne1,1)=-NweT*(1./yeg.*wyeg);
  RhsL(ir,1)=RhsL(ir,1)+fr_axisym;
  
  Kwr_axisym=zeros(nsd*NumElementNodes,NumElementNodes);
  Kwr_axisym(ne1,ne1)=-NweT*((1./yeg.*vxeg.*vyeg).*Ne);
  Kwr_axisym(ne2,ne1)=-NweT*((1./yeg.*vyeg.*vyeg).*Ne);
  LhsLL(iw,ir)=LhsLL(iw,ir)+Kwr_axisym;
  
  Kww_axisym=zeros(nsd*NumElementNodes,nsd*NumElementNodes);
  Kww_axisym(ne1,ne1)=+NweT*((1./yeg.*vyeg).*Ne);
  Kww_axisym(ne1,ne2)=+NweT*((1./yeg.*vxeg).*Ne);
  Kww_axisym(ne2,ne2)=+NweT*((2./yeg.*vyeg).*Ne);
  LhsLL(iw,iw)=LhsLL(iw,iw)+Kww_axisym;
  
  fw_axisym=zeros(nsd*NumElementNodes,1);
  fw_axisym(ne1,1)=-NweT*(1./yeg.*wxeg.*vyeg);
  fw_axisym(ne2,1)=-NweT*(1./yeg.*wyeg.*vyeg);
  RhsL(iw,1)=RhsL(iw,1)+fw_axisym;
  
  Kgr_axisym=zeros(NumElementNodes,NumElementNodes);
  Kgw_axisym=zeros(NumElementNodes,nsd*NumElementNodes);
  Kgg_axisym=zeros(NumElementNodes,NumElementNodes);
  Kgr_axisym(ne1,ne1)=+NweT*((1./yeg.*(-gamma*g_reg+(gamma-1)*(vxeg.^2+vyeg.^2)).*vyeg).*Ne);
  Kgw_axisym(ne1,ne1)=+NweT*((-1./yeg*(gamma-1).*wxeg.*wyeg./reg.^2).*Ne);
  Kgw_axisym(ne1,ne2)=+NweT*((1./yeg.*(gamma*g_reg-(gamma-1)*1/2*(vxeg.^2+3*vyeg.^2))).*Ne);
  Kgg_axisym(ne1,ne1)=+NweT*((1./yeg*gamma.*vyeg).*Ne);
  LhsLL(ig,ir)=LhsLL(ig,ir)+Kgr_axisym;
  LhsLL(ig,iw)=LhsLL(ig,iw)+Kgw_axisym;
  LhsLL(ig,ig)=LhsLL(ig,ig)+Kgg_axisym;
  
  fg_axisym=zeros(NumElementNodes,1);
  fg_axisym(ne1,1)=-NweT*(1./yeg.*(geg+(gamma-1)*(geg-1/2*(vxeg.^2+vyeg.^2).*reg)).*vyeg);
  RhsL(ig,1)=RhsL(ig,1)+fg_axisym;
  
  KeH_axisym=zeros(nsd*NumElementNodes,qsd*NumElementNodes);
  KeH_axisym(ne1,ne1)=-NweT*((1./yeg).*Ne);
  LhsLL(ie,iH)=LhsLL(ie,iH)+KeH_axisym;
  
  fe_axisym=zeros(nsd*NumElementNodes,1);
  fe_axisym(ne1,1)=+NweT*(1./yeg.*Hzeg);
  RhsL(ie,1)=RhsL(ie,1)+fe_axisym;
end

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

% Add zero entries if outside PML region
if not(isPML)
  MatLocal(end+(1:(qsd+nsd)*NumElementNodes),:)=0;
  VecLocal(end+(1:(qsd+nsd)*NumElementNodes),1)=0;
end

end

%% Do post-process element
function [LhsPost,RhsPost]=doPostProcessElement(...
  Nodes,SolutionLocal,SolutionOld,Parameters,Time,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
NumElementNodesPost=size(RefElement.Post.NodesCoordElem,1);
NumElementNodesPostHigh=size(RefElement.PostHighHigh.NodesCoordElem,1);
mu=Parameters.MagneticPermeability;
Sigma=Parameters.DampingFunction;
Xe=Nodes';
dt=Time.TimeStepSize;
BDFo=Time.BDFOrderEff;
alpha=Time.BDF1stDerEff;

% Get solution
He=reshape(SolutionLocal(:,1+nsd+1+(1:qsd)),[],1);
ee=reshape(SolutionLocal(:,1+nsd+1+qsd+(1:nsd)),[],1);
Me=reshape(SolutionLocal(:,1+nsd+1+qsd+nsd+(1:qsd)),[],1);
Holde=reshape(SolutionOld(:,1+nsd+1+(1:qsd),:),[],BDFo);

% Initialize lhs
Kpp=zeros(nsd*NumElementNodesPost,nsd*NumElementNodesPost);
Khp=zeros(NumElementNodesPostHigh,nsd*NumElementNodesPost);

% Initialize rhs
fp=zeros(nsd*NumElementNodesPost,1);
fh=zeros(NumElementNodesPostHigh,1);

% Get reference data
Nhe=RefElement.PostHigh.ShapeFunctionsElem;

% Compute weights at Gauss points
[~,Nex,Ney,Nez,weg,Nle]=mapShapeFunctions(1,RefElement.PostLow,RefElement.Post,Xe,nsd);
[~,Nhhex,Nhhey,Nhhez,wheg,Nlhe]=mapShapeFunctions(1,RefElement.PostLowHigh,...
                                                    RefElement.PostHighHigh,Xe,nsd);

% Indices
ne1=1:NumElementNodesPost;
ne2=ne1+NumElementNodesPost;
ne3=ne2+NumElementNodesPost;
nle1=1:NumElementNodes;
nle2=nle1+NumElementNodes;
nle3=nle2+NumElementNodes;
nhe1=1:NumElementNodesPostHigh;

% Compute variables at nodes
if nsd==2
  Hze=He(nle1);
elseif nsd==3
  Hxe=He(nle1);
  Hye=He(nle2);
  Hze=He(nle3);
end
exe=ee(nle1);
eYe=ee(nle2);
if nsd==3
  eze=ee(nle3);
end
if nsd==2
  Mze=Me(nle1);
elseif nsd==3
  Mxe=Me(nle1);
  Mye=Me(nle2);
  Mze=Me(nle3);
end
if nsd==2
  Holdze=Holde(nle1,:);
elseif nsd==3
  Holdxe=Holde(nle1,:);
  Holdye=Holde(nle2,:);
  Holdze=Holde(nle3,:);
end
Sigmae=Sigma(Xe(nle1,1),Xe(nle1,2),Xe(nle1,3));
sigmaxe=Sigmae(:,1);
sigmaye=Sigmae(:,2);
if nsd==3
  sigmaze=Sigmae(:,3);
end

% Compute variables at Gauss points
Xeg=Nle*Xe;
if nsd==2
  Hzeg=Nle*Hze;
elseif nsd==3
  Hxeg=Nle*Hxe;
  Hyeg=Nle*Hye;
  Hzeg=Nle*Hze;
end
exeg=Nlhe*exe;
eyeg=Nlhe*eYe;
if nsd==3
  ezeg=Nlhe*eze;
end
if nsd==2
  Mzeg=Nle*Mze;
elseif nsd==3
  Mxeg=Nle*Mxe;
  Myeg=Nle*Mye;
  Mzeg=Nle*Mze;
end
if nsd==2
  Holdzeg=Nle*Holdze;
elseif nsd==3
  Holdxeg=Nle*Holdxe;
  Holdyeg=Nle*Holdye;
  Holdzeg=Nle*Holdze;
end
sigmaxeg=Nle*sigmaxe;
sigmayeg=Nle*sigmaye;
if nsd==3
  sigmazeg=Nle*sigmaze;
end
if isa(mu,'function_handle')
  mueg=mu(Xeg(:,1),Xeg(:,2),Xeg(:,3));
else
  mueg=mu;
end

% Compute basic matrices
NwexT=(weg.*Nex)';
NweyT=(weg.*Ney)';
if nsd==3
  NwezT=(weg.*Nez)';
end
NwhhexT=(wheg.*Nhhex)';
NwhheyT=(wheg.*Nhhey)';
if nsd==3
  NwhhezT=(wheg.*Nhhez)';
end
if nsd==2
  Kxxe=NwexT*Nex; Kxye=NwexT*Ney;
  Kyxe=NweyT*Nex; Kyye=NweyT*Ney;
elseif nsd==3
  Kxxe=NwexT*Nex; Kxye=NwexT*Ney; Kxze=NwexT*Nez;
  Kyxe=NweyT*Nex; Kyye=NweyT*Ney; Kyze=NweyT*Nez;
  Kzxe=NwezT*Nex; Kzye=NwezT*Ney; Kzze=NwezT*Nez;
end

% Compute lhs
if nsd==2
  Kpp(ne1,ne1)=-Kyye;
  Kpp(ne1,ne2)=+Kyxe;
  Kpp(ne2,ne1)=+Kxye;
  Kpp(ne2,ne2)=-Kxxe;
elseif nsd==3
  Kpp(ne1,ne1)=-Kyye-Kzze;
  Kpp(ne1,ne2)=+Kyxe;
  Kpp(ne1,ne3)=+Kzxe;
  Kpp(ne2,ne1)=+Kxye;
  Kpp(ne2,ne2)=-Kxxe-Kzze;
  Kpp(ne2,ne3)=+Kzye;
  Kpp(ne3,ne1)=+Kxze;
  Kpp(ne3,ne2)=+Kyze;
  Kpp(ne3,ne3)=-Kxxe-Kyye;
end

Khp(nhe1,ne1)=NwhhexT*(Nhe);
Khp(nhe1,ne2)=NwhheyT*(Nhe);
if nsd==3
  Khp(nhe1,ne3)=NwhhezT*(Nhe);
end

% Compute rhs
if nsd==2
  fp(ne1,1)=-NweyT*(mueg.*(Hzeg*alpha(1)+Holdzeg*alpha(2:BDFo+1,1))/dt);
  fp(ne2,1)=+NwexT*(mueg.*(Hzeg*alpha(1)+Holdzeg*alpha(2:BDFo+1,1))/dt);
elseif nsd==3
  fp(ne1,1)=+NwezT*(mueg.*(Hyeg*alpha(1)+Holdyeg*alpha(2:BDFo+1,1))/dt)...
            -NweyT*(mueg.*(Hzeg*alpha(1)+Holdzeg*alpha(2:BDFo+1,1))/dt);
  fp(ne2,1)=+NwexT*(mueg.*(Hzeg*alpha(1)+Holdzeg*alpha(2:BDFo+1,1))/dt)...
            -NwezT*(mueg.*(Hxeg*alpha(1)+Holdxeg*alpha(2:BDFo+1,1))/dt);
  fp(ne3,1)=+NweyT*(mueg.*(Hxeg*alpha(1)+Holdxeg*alpha(2:BDFo+1,1))/dt)...
            -NwexT*(mueg.*(Hyeg*alpha(1)+Holdyeg*alpha(2:BDFo+1,1))/dt);
end

if nsd==2
  fp(ne1,1)=fp(ne1,1)-NweyT*(mueg.*(sigmaxeg+sigmayeg).*Hzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(mueg.*(sigmaxeg+sigmayeg).*Hzeg);
elseif nsd==3
  fp(ne1,1)=fp(ne1,1)+NwezT*(mueg.*(sigmazeg+sigmaxeg-sigmayeg).*Hyeg)...
                     -NweyT*(mueg.*(sigmaxeg+sigmayeg-sigmazeg).*Hzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(mueg.*(sigmaxeg+sigmayeg-sigmazeg).*Hzeg)...
                     -NwezT*(mueg.*(sigmayeg+sigmazeg-sigmaxeg).*Hxeg);
  fp(ne3,1)=fp(ne3,1)+NweyT*(mueg.*(sigmayeg+sigmazeg-sigmaxeg).*Hxeg)...
                     -NwexT*(mueg.*(sigmazeg+sigmaxeg-sigmayeg).*Hyeg);
end

if nsd==2
  fp(ne1,1)=fp(ne1,1)-NweyT*(Mzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(Mzeg);
elseif nsd==3
  fp(ne1,1)=fp(ne1,1)+NwezT*(Myeg)...
                     -NweyT*(Mzeg);
  fp(ne2,1)=fp(ne2,1)+NwexT*(Mzeg)...
                     -NwezT*(Mxeg);
  fp(ne3,1)=fp(ne3,1)+NweyT*(Mxeg)...
                     -NwexT*(Myeg);
end

fh(nhe1,1)=NwhhexT*(exeg)...
          +NwhheyT*(eyeg);
if nsd==3
  fh(nhe1,1)=fh(nhe1,1)+NwhhezT*(ezeg);
end

% Indices
ip=1:nsd*NumElementNodesPost;
ih=ip(end)+(1:NumElementNodesPostHigh);

% Initialization of lhs and rhs
LhsPost=zeros(nsd*NumElementNodesPost+NumElementNodesPostHigh,nsd*NumElementNodesPost);
RhsPost=zeros(nsd*NumElementNodesPost+NumElementNodesPostHigh,1);

% Lhs for post-processing
LhsPost(ip,ip)=Kpp;
LhsPost(ih,ip)=Khp;

% Rhs for post-processing
RhsPost(ip,1)=fp;
RhsPost(ih,1)=fh;

end

%% Compute radiation element
function [Radiation]=computeRadiationElement(...
  iFace,Nodes,ResultsFFT,RefElement,Sizes)

% Get general parameters
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
Xe=Nodes';

% Get solution
Hwe=reshape(ResultsFFT(:,1:qsd),[],1);
ewe=reshape(ResultsFFT(:,qsd+(1:nsd)),[],1);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Hwze=Hwe(ne1);
elseif nsd==3
  Hwxe=Hwe(ne1);
  Hwye=Hwe(ne2);
  Hwze=Hwe(ne3);
end
ewxe=ewe(ne1);
ewye=ewe(ne2);
if nsd==3
  ewze=ewe(ne3);
end

% Compute weights at Gauss points
FaceNodes=RefElement.FaceNodesElem;
Xf=Xe(FaceNodes(iFace,:),:);
[Nf,nx,ny,nz,wfg]=mapShapeFunctions(0,RefElement,RefElement,Xf(1:nsd,:),nsd);

% Indices
nf1=FaceNodes(iFace,:);

% Compute variables at nodes
if nsd==2
  Hwzf=Hwze(nf1);
elseif nsd==3
  Hwxf=Hwxe(nf1);
  Hwyf=Hwye(nf1);
  Hwzf=Hwze(nf1);
end
ewxf=ewxe(nf1);
ewyf=ewye(nf1);
if nsd==3
  ewzf=ewze(nf1);
end

% Compute variables at Gauss points
if nsd==2
  Hwzfg=Nf*Hwzf;
elseif nsd==3
  Hwxfg=Nf*Hwxf;
  Hwyfg=Nf*Hwyf;
  Hwzfg=Nf*Hwzf;
end
ewxfg=Nf*ewxf;
ewyfg=Nf*ewyf;
if nsd==3
  ewzfg=Nf*ewzf;
end

% Poynting vector at Gauss points
if nsd==2
  Sxfg=real(+ewyfg.*conj(Hwzfg));
  Syfg=real(-ewxfg.*conj(Hwzfg));
elseif nsd==3
  Sxfg=real(ewyfg.*conj(Hwzfg)-ewzfg.*conj(Hwyfg));
  Syfg=real(ewzfg.*conj(Hwxfg)-ewxfg.*conj(Hwzfg));
  Szfg=real(ewxfg.*conj(Hwyfg)-ewyfg.*conj(Hwxfg));
end

% Compute radiation
if nsd==2
  Radiation=sum((Sxfg.*nx+Syfg.*ny).*wfg);
elseif nsd==3
  Radiation=sum((Sxfg.*nx+Syfg.*ny+Szfg.*nz).*wfg);
end

end

%% Compute divergence element
function [MagneticFieldDiv,ElectricFieldDiv]=computeDivergenceElement(...
  Nodes,SolutionLocal,Parameters,RefElement,Sizes)

% Get general parameters
isAxisymmetric=strcmp(Parameters.Axisymmetric,'yes');
nsd=Sizes.NumSpaceDim;
msd=nsd*(nsd+1)/2;
qsd=msd-nsd;
NumElementNodes=Sizes.NumElementNodes;
Xe=Nodes';

% Get solution
He=reshape(SolutionLocal(:,1+nsd+1+(1:qsd)),[],1);
ee=reshape(SolutionLocal(:,1+nsd+1+qsd+(1:nsd)),[],1);

% Compute weights at Gauss points
[Ne,Nex,Ney,Nez,~,~,pinvNe]=mapShapeFunctions(1,RefElement,RefElement,Xe,nsd);

% Indices
ne1=1:NumElementNodes;
ne2=ne1+NumElementNodes;
ne3=ne2+NumElementNodes;

% Compute variables at nodes
if nsd==2
  Hze=He(ne1);
elseif nsd==3
  Hxe=He(ne1);
  Hye=He(ne2);
  Hze=He(ne3);
end
exe=ee(ne1);
eYe=ee(ne2);
if nsd==3
  eze=ee(ne3);
end

% Compute divergence of magnetic field
if nsd==2
  MagneticFieldDiv=zeros(NumElementNodes,1);
elseif nsd==3
  MagneticFieldDiv=pinvNe*(Nex*Hxe+Ney*Hye+Nez*Hze);
end

% Compute divergence of electric field
if nsd==2
  ElectricFieldDiv=pinvNe*(Nex*exe+Ney*eYe);
elseif nsd==3
  ElectricFieldDiv=pinvNe*(Nex*exe+Ney*eYe+Nez*eze);
end

% Add missing term for axisymmetric formulation
if isAxisymmetric
  ElectricFieldDiv=ElectricFieldDiv+pinvNe*(1./(Ne*Xe(:,2)).*(Ne*eYe));
end

end