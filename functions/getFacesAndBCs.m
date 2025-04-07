function [Faces,Elements,BCs]=getFacesAndBCs(...
         Simulation,Geometry,Mesh,Boundaries,Regions,Parameters,RefElement,Sizes)
         % Get faces and BCs

% Initialize structures
Faces=struct();
Elements=struct();
BCs=struct();
for iD=1:Simulation.NumDiscretizations
  % Get info
  FaceNodesElem=RefElement(iD,iD).Degree1.FaceNodesElem;
  NumElements=size(Mesh(iD).Degree1.Elements,2);
  NumElemNodes=size(Mesh(iD).Degree1.Elements,1);
  NumElemFaces=size(FaceNodesElem,1);
  NumFaceNodes=size(FaceNodesElem,2);
  if Sizes(iD).NumSpaceDim==2
    NumFaces=Geometry(iD).NumEdges;
  elseif Sizes(iD).NumSpaceDim==3
    NumFaces=Geometry(iD).NumFaces;
  end
  NumGlobalComp=Sizes(iD).NumGlobalComp;
  NumFaceNodesHigh=size(RefElement(iD,iD).FaceNodesElem,2);
  
  % Get interior faces
  if strcmp(Parameters(iD).DiscretizationType,'HDG')
    Faces(iD,iD).Interior=zeros(NumElements*NumElemFaces,5);
    Elements(iD).Faces.Interior=mat2cell(repmat(zeros(2,NumElemFaces),1,NumElements),2,...
                                                               repmat(NumElemFaces,1,NumElements))';
    iFaceInt=0;
    MeshElementsFast=Mesh(iD).Degree1.Elements';
    for iElem1=1:NumElements
      MeshElementsFast(1,:)=[];
      for iElemFace1=1:NumElemFaces
        FaceNodes1=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFace1,:),iElem1)';
        [Filter1,~]=find(MeshElementsFast==FaceNodes1(1));
        [Filter2,~]=find(MeshElementsFast(Filter1,:)==FaceNodes1(2));
        if Sizes(iD).NumSpaceDim==3
          [Filter3,~]=find(MeshElementsFast(Filter1(Filter2),:)==FaceNodes1(3));
        else
          Filter3=1:length(Filter2);
        end
        iElem2=iElem1+Filter1(Filter2(Filter3));
        if not(isempty(iElem2))
          iFaceInt=iFaceInt+1;
          [iFaceNodes2,~]=find(Mesh(iD).Degree1.Elements(:,iElem2)==FaceNodes1);
          iElemFace2=find(all(sort(FaceNodesElem,2)==sort(iFaceNodes2)',2));
          FaceNodes2=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFace2,:),iElem2)';
          Node2Match1stNode1=find(FaceNodes2==FaceNodes1(1));
          Faces(iD,iD).Interior(iFaceInt,:)=[iElem1,iElemFace1,...
                                             iElem2,iElemFace2,Node2Match1stNode1];
          Elements(iD).Faces.Interior{iElem1,1}(:,iElemFace1)=[1;0];
          Elements(iD).Faces.Interior{iElem2,1}(:,iElemFace2)=[1;Node2Match1stNode1];
        end
      end
    end
    Faces(iD,iD).Interior=Faces(iD,iD).Interior(1:iFaceInt,:);
  end
  
  % Get exterior faces and BCs
  Faces(iD,iD).Exterior=double.empty(0,2);
  Elements(iD).Faces.Exterior=...
        mat2cell(repmat(zeros(1,NumElemFaces),1,NumElements),1,repmat(NumElemFaces,1,NumElements))';
  Elements(iD).Faces.Boundary=...
        mat2cell(repmat(zeros(1,NumElemFaces),1,NumElements),1,repmat(NumElemFaces,1,NumElements))';
  BoundaryNames=fieldnames(Boundaries(iD));
  for iBoundaryName=1:length(BoundaryNames)
    Faces(iD,iD).(BoundaryNames{iBoundaryName})=[];
    Elements(iD).Faces.(BoundaryNames{iBoundaryName})=...
        mat2cell(repmat(zeros(1,NumElemFaces),1,NumElements),1,repmat(NumElemFaces,1,NumElements))';  
    BCs(iD).(BoundaryNames{iBoundaryName})=[];
  end
  for iBoundaryFace=1:NumFaces
    if Sizes(iD).NumSpaceDim==2
      Nodes=(findNodes(Mesh(iD).Degree1,'Region','Edge',iBoundaryFace))';
    elseif Sizes(iD).NumSpaceDim==3
      Nodes=(findNodes(Mesh(iD).Degree1,'Region','Face',iBoundaryFace))';
    end
    ElementNodes=ismember(Mesh(iD).Degree1.Elements',Nodes);
    Element=find(sum(ElementNodes,2)==NumFaceNodes);
    A=(ElementNodes(Element,:)*diag(1:NumElemNodes))'; A(A==0)=[];
    A=reshape(A,size(FaceNodesElem,2),length(Element))'; B=sort(FaceNodesElem,2);
    [~,ElementFaces]=ismember(A,B,'rows');
    if strcmp(Parameters(iD).DiscretizationType,'HDG')
      ElementMore=repelem(find(sum(ElementNodes,2)>NumFaceNodes),NumElemFaces,1);
      ElementFacesMore=repmat((1:NumElemFaces)',numel(ElementMore)/NumElemFaces,1);
      Element=[Element;ElementMore];                %#ok
      ElementFaces=[ElementFaces;ElementFacesMore]; %#ok
      Interior=ismember([Element,ElementFaces],Faces(iD,iD).Interior(:,1:2),'rows') | ...
               ismember([Element,ElementFaces],Faces(iD,iD).Interior(:,3:4),'rows');
      Element=Element(not(Interior));
      ElementFaces=ElementFaces(not(Interior));
    end
    Faces(iD,iD).Exterior=[Faces(iD,iD).Exterior;Element,ElementFaces];
    for iElem=1:length(Element)
      Elements(iD).Faces.Exterior{Element(iElem),1}(ElementFaces(iElem))=1;
      Elements(iD).Faces.Boundary{Element(iElem),1}(ElementFaces(iElem))=iBoundaryFace;
    end
    for iBoundaryName=1:length(BoundaryNames)
      if ismember(iBoundaryFace,Boundaries(iD).(BoundaryNames{iBoundaryName}))
        Faces(iD,iD).(BoundaryNames{iBoundaryName})=...
                                 [Faces(iD,iD).(BoundaryNames{iBoundaryName});Element,ElementFaces];
        for iElem=1:length(Element)
          Elements(iD).Faces.(BoundaryNames{iBoundaryName}){Element(iElem),1}(ElementFaces(iElem))=1;
        end
        BCs(iD).(BoundaryNames{iBoundaryName})=[BCs(iD).(BoundaryNames{iBoundaryName});Nodes];
      end
    end
  end
  Faces(iD,iD).Exterior=sortrows(Faces(iD,iD).Exterior);
  for iBoundaryName=1:length(BoundaryNames)
    if not(isempty(Faces(iD,iD).(BoundaryNames{iBoundaryName})))
      [~,OrderFaces]=sort(Faces(iD,iD).(BoundaryNames{iBoundaryName})(:,1));
      Faces(iD,iD).(BoundaryNames{iBoundaryName})=...
                                          Faces(iD,iD).(BoundaryNames{iBoundaryName})(OrderFaces,:);
      BCs(iD).(BoundaryNames{iBoundaryName})=unique(BCs(iD).(BoundaryNames{iBoundaryName}));
    else
      Faces(iD,iD).(BoundaryNames{iBoundaryName})=double.empty(0,2);
    end
  end
  
  % Get faces of specific regions
  if not(isempty(Regions))
    RegionsNames=fieldnames(Regions(iD));
    for iRegionName=1:length(RegionsNames)
      Faces(iD,iD).(RegionsNames{iRegionName})=[];
      Elements(iD).Faces.(RegionsNames{iRegionName})=...
        mat2cell(repmat(zeros(1,NumElemFaces),1,NumElements),1,repmat(NumElemFaces,1,NumElements))';
      BCs(iD).(RegionsNames{iRegionName})=[];
      if Sizes(iD).NumSpaceDim==2
        Nodes=find(Regions(iD).(RegionsNames{iRegionName}).Definition(...
          Mesh(iD).Degree1.Nodes(1,:),Mesh(iD).Degree1.Nodes(2,:)))';
      elseif Sizes(iD).NumSpaceDim==3
        Nodes=find(Regions(iD).(RegionsNames{iRegionName}).Definition(...
          Mesh(iD).Degree1.Nodes(1,:),Mesh(iD).Degree1.Nodes(2,:),Mesh(iD).Degree1.Nodes(3,:)))';                                             
      end                  
      ElementNodes=ismember(Mesh(iD).Degree1.Elements',Nodes);
      ElementSharing=find(sum(ElementNodes,2)==NumFaceNodes);
      Element=ElementSharing;
      for iElem=1:length(ElementSharing)
        ElemNodes=Mesh(iD).Degree1.Elements(:,ElementSharing(iElem));
        ElemNodesCoord=Mesh(iD).Degree1.Nodes(:,ElemNodes);
        ElemCenterCoord=mean(ElemNodesCoord,2);
        if Sizes(iD).NumSpaceDim==2
          if not(Regions(iD).(RegionsNames{iRegionName}).Condition(...
              ElemCenterCoord(1),ElemCenterCoord(2)))
            Element(Element==ElementSharing(iElem))=[];
          end
        elseif Sizes(iD).NumSpaceDim==3
          if not(Regions(iD).(RegionsNames{iRegionName}).Condition(...
              ElemCenterCoord(1),ElemCenterCoord(2),ElemCenterCoord(3)))
            Element(Element==ElementSharing(iElem))=[];
          end
        end
      end
      A=(ElementNodes(Element,:)*diag(1:NumElemNodes))'; A(A==0)=[];
      A=reshape(A,size(FaceNodesElem,2),length(Element))'; B=sort(FaceNodesElem,2);
      [~,ElementFaces]=ismember(A,B,'rows');
      Faces(iD,iD).(RegionsNames{iRegionName})=...
        [Faces(iD,iD).(RegionsNames{iRegionName});Element,ElementFaces];
      for iElem=1:length(Element)
        Elements(iD).Faces.(RegionsNames{iRegionName}){Element(iElem),1}(ElementFaces(iElem))=1;
      end
      BCs(iD).(RegionsNames{iRegionName})=[BCs(iD).(RegionsNames{iRegionName});Nodes];
    end
  end
  
  % Couple faces at the interface
  if matchField(Boundaries(iD),'Interface')
    for iFace=1:size(Faces(iD,iD).Interface,1)
      iElem=Faces(iD,iD).Interface(iFace,1);
      iElemFace=Faces(iD,iD).Interface(iFace,2);
      Nodes=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFace,:),iElem);
      NodesCoord=Mesh(iD).Degree1.Nodes(:,Nodes)';
      iDOthers=setdiff(1:Simulation.NumDiscretizations,iD);
      for iDAux=1:length(iDOthers)
        iDC=iDOthers(iDAux);
        FaceNodesElemC=RefElement(iDC,iDC).Degree1.FaceNodesElem;
        NodesCoordC=NodesCoord;
        Tolerance=Mesh(iDC).Degree1.MinElementSize/1000;
        [~,NodesC]=ismembertol(NodesCoordC,Mesh(iDC).Degree1.Nodes',Tolerance,'ByRows',true);
        ElementNodesC=ismember(Mesh(iDC).Degree1.Elements',NodesC);
        iElemC=find(sum(ElementNodesC,2)==Sizes(iDC).NumSpaceDim);
        A=(ElementNodesC(iElemC,:)*diag(1:size(ElementNodesC,2)));
        A(A==0)=[]; B=sort(FaceNodesElemC,2);
        [~,iElemFaceC]=ismember(A,B,'rows');
        FaceNodesC=Mesh(iDC).Degree1.Elements(FaceNodesElemC(iElemFaceC,:),iElemC)';
        Node2Match1stNode1=find(FaceNodesC==NodesC(1));
        Faces(iD,iDC).Interface(iFace,:)=[iElem,iElemFace,iElemC,iElemFaceC,Node2Match1stNode1];
        Elements(iD).Faces.Interface{iElem,1}(2:3,iElemFace)=[iElemFaceC;Node2Match1stNode1];
      end
    end
  end
  
  % Couple faces for periodic boundary conditions
  if matchField(Boundaries(iD),'PeriodicMaster')
    Faces(iD,iD).Periodic=zeros(size(Faces(iD,iD).PeriodicMaster,1),5);
    Elements(iD).Faces.Periodic=mat2cell(repmat(zeros(2,NumElemFaces),1,NumElements),2,...
                                                               repmat(NumElemFaces,1,NumElements))';
    Tolerance=Mesh(iD).Degree1.MinElementSize/1000;
    for iFaceM=1:size(Faces(iD,iD).PeriodicMaster,1)
      iElemM=Faces(iD,iD).PeriodicMaster(iFaceM,1);
      iElemFaceM=Faces(iD,iD).PeriodicMaster(iFaceM,2);
      NodesM=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFaceM,:),iElemM);
      NodesCoordMxyz=Mesh(iD).Degree1.Nodes(:,NodesM)';
      [~,Dim]=min(max(NodesCoordMxyz)-min(NodesCoordMxyz));
      NodesCoordM=NodesCoordMxyz(:,setdiff(1:Sizes(iD).NumSpaceDim,Dim));
      for iFaceS=1:size(Faces(iD,iD).PeriodicSlave,1)
        iElemS=Faces(iD,iD).PeriodicSlave(iFaceS,1);
        iElemFaceS=Faces(iD,iD).PeriodicSlave(iFaceS,2);
        NodesS=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFaceS,:),iElemS);
        NodesCoordS=Mesh(iD).Degree1.Nodes(setdiff(1:Sizes(iD).NumSpaceDim,Dim),NodesS)';
        [NodesMatch,NodesOrder]=ismembertol(NodesCoordS,NodesCoordM,Tolerance,'ByRows',true);
        if all(NodesMatch) && length(unique(NodesOrder))==length(NodesMatch)
          break
        end
      end
      Node2Match1stNode1=NodesOrder(1);
      Faces(iD,iD).Periodic(iFaceM,:)=[iElemM,iElemFaceM,iElemS,iElemFaceS,Node2Match1stNode1];
      Elements(iD).Faces.Periodic{iElemM,1}(:,iElemFaceM)=[1;0];
      Elements(iD).Faces.Periodic{iElemS,1}(:,iElemFaceS)=[1;Node2Match1stNode1];
    end
    if any(Faces(iD,iD).Periodic(:,5)==0)
      error('Not all periodic faces are properly coupled!');
    end
  end
  
  % Build faces connectivity
  if strcmp(Parameters(iD).DiscretizationType,'HDG')
    Faces(iD,iD).Connectivity=zeros(NumElements,NumElemFaces);
    for iFace=1:size(Faces(iD,iD).Interior,1)
      Int=Faces(iD,iD).Interior(iFace,:);
      Faces(iD,iD).Connectivity(Int(1),Int(2))=iFace;
      Faces(iD,iD).Connectivity(Int(3),Int(4))=iFace;
    end
    for iFace=1:size(Faces(iD,iD).Exterior,1)
      Ext=Faces(iD,iD).Exterior(iFace,:);
      Faces(iD,iD).Connectivity(Ext(1),Ext(2))=iFace+size(Faces(iD,iD).Interior,1);
    end
  end
  
  % Modify faces connectivity for periodic boundary conditions
  if strcmp(Parameters(iD).DiscretizationType,'HDG') && matchField(Faces(iD),'Periodic')
   for iFace=1:size(Faces(iD,iD).Periodic,1)
     Per=Faces(iD,iD).Periodic(iFace,:);
     Faces(iD,iD).Connectivity(Per(3),Per(4))=Faces(iD,iD).Connectivity(Per(1),Per(2));
   end
   Faces(iD,iD).Connectivity=changem(Faces(iD,iD).Connectivity,...
     1:length(unique(Faces(iD,iD).Connectivity)),unique(Faces(iD,iD).Connectivity));
  end
  
  % Build global-local map
  Faces(iD,iD).GlobalLocal=cell(NumElements,1);
  if strcmp(Parameters(iD).DiscretizationType,'HDG')
    Faces(iD,iD).GlobalLocal=mat2cell(reshape(((repmat(reshape(Faces(iD,iD).Connectivity',[],1),...
        1,NumGlobalComp*NumFaceNodesHigh)-1)*NumGlobalComp*NumFaceNodesHigh...
        +(1:NumGlobalComp*NumFaceNodesHigh))',[],1),...
        repmat(NumGlobalComp*NumElemFaces*NumFaceNodesHigh,NumElements,1));
  end
end

% Convert element faces into a better structure
for iD=1:Simulation.NumDiscretizations
  Elements(iD).Faces=table2struct(struct2table(Elements(iD).Faces));
end

end