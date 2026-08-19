function [Faces,Elements,BCs]=getFacesAndBCs(Simulation,Geometry,Mesh,Boundaries,Regions,...
  Parameters,RefElement,Sizes)
% Get the interior, exterior, boundary, region, interface and periodic faces of each discretization
%
% Faces collects, for each set of faces, the pairs [element,local face] it is made of, Elements
% flags the same information element by element for an effective use of parfor, and BCs collects
% the mesh nodes lying on each boundary and region.

Faces=struct();
Elements=struct();
BCs=struct();
for iD=1:Simulation.NumDiscretizations

  % Mesh and reference element data
  FaceNodesElem=RefElement(iD,iD).Degree1.FaceNodesElem;
  NumElements=size(Mesh(iD).Degree1.Elements,2);
  NumElementNodes=size(Mesh(iD).Degree1.Elements,1);
  NumElementFaces=size(FaceNodesElem,1);
  NumFaceNodes=size(FaceNodesElem,2);
  NumFaceNodesHigh=size(RefElement(iD,iD).NodesCoordFace,1);
  NumGlobalComp=Sizes(iD).NumGlobalComp;
  if Sizes(iD).NumSpaceDim==2
    NumBoundaries=Geometry(iD).NumEdges;
  elseif Sizes(iD).NumSpaceDim==3
    NumBoundaries=Geometry(iD).NumFaces;
  end
  HDG=strcmp(Parameters(iD).DiscretizationType,'HDG');
  Flags=repmat({zeros(1,NumElementFaces)},NumElements,1);
  FlagsPaired=repmat({zeros(2,NumElementFaces)},NumElements,1);

  % Interior faces, shared by two elements
  if HDG
    Faces(iD,iD).Interior=zeros(NumElements*NumElementFaces,5);
    Elements(iD).Faces.Interior=FlagsPaired;
    iFaceInt=0;
    ElementsAhead=Mesh(iD).Degree1.Elements';
    for iElem1=1:NumElements
      ElementsAhead(1,:)=[];
      for iElemFace1=1:NumElementFaces
        FaceNodes1=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFace1,:),iElem1)';
        [Match1,~]=find(ElementsAhead==FaceNodes1(1));
        [Match2,~]=find(ElementsAhead(Match1,:)==FaceNodes1(2));
        if Sizes(iD).NumSpaceDim==3
          [Match3,~]=find(ElementsAhead(Match1(Match2),:)==FaceNodes1(3));
        else
          Match3=1:length(Match2);
        end
        iElem2=iElem1+Match1(Match2(Match3));
        if not(isempty(iElem2))
          iFaceInt=iFaceInt+1;
          [iFaceNodes2,~]=find(Mesh(iD).Degree1.Elements(:,iElem2)==FaceNodes1);
          iElemFace2=find(all(sort(FaceNodesElem,2)==sort(iFaceNodes2)',2));
          FaceNodes2=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFace2,:),iElem2)';
          Node2Match1stNode1=find(FaceNodes2==FaceNodes1(1));
          Faces(iD,iD).Interior(iFaceInt,:)=...
            [iElem1,iElemFace1,iElem2,iElemFace2,Node2Match1stNode1];
          Elements(iD).Faces.Interior{iElem1,1}(:,iElemFace1)=[1;0];
          Elements(iD).Faces.Interior{iElem2,1}(:,iElemFace2)=[1;Node2Match1stNode1];
        end
      end
    end
    Faces(iD,iD).Interior=Faces(iD,iD).Interior(1:iFaceInt,:);
  end

  % Exterior faces, split into the boundaries named in the input deck
  Faces(iD,iD).Exterior=zeros(0,2);
  Elements(iD).Faces.Exterior=Flags;
  Elements(iD).Faces.Boundary=Flags;
  BoundaryNames=fieldnames(Boundaries(iD));
  for iBN=1:numel(BoundaryNames)
    Faces(iD,iD).(BoundaryNames{iBN})=[];
    Elements(iD).Faces.(BoundaryNames{iBN})=Flags;
    BCs(iD).(BoundaryNames{iBN})=[];
  end
  for iBoundary=1:NumBoundaries
    if Sizes(iD).NumSpaceDim==2
      Nodes=(findNodes(Mesh(iD).Degree1,'Region','Edge',iBoundary))';
    elseif Sizes(iD).NumSpaceDim==3
      Nodes=(findNodes(Mesh(iD).Degree1,'Region','Face',iBoundary))';
    end
    ElementNodes=ismember(Mesh(iD).Degree1.Elements',Nodes);
    Element=find(sum(ElementNodes,2)==NumFaceNodes);
    A=(ElementNodes(Element,:)*diag(1:NumElementNodes))'; A(A==0)=[];
    A=reshape(A,NumFaceNodes,length(Element))';
    [~,ElementFaces]=ismember(A,sort(FaceNodesElem,2),'rows');
    if HDG
      % Elements with more nodes on the boundary own several exterior faces
      ElementMore=repelem(find(sum(ElementNodes,2)>NumFaceNodes),NumElementFaces,1);
      Element=[Element;ElementMore];                                                  %#ok
      ElementFaces=[ElementFaces;repmat((1:NumElementFaces)',...
                                        numel(ElementMore)/NumElementFaces,1)];       %#ok
      Interior=ismember([Element,ElementFaces],Faces(iD,iD).Interior(:,1:2),'rows') | ...
               ismember([Element,ElementFaces],Faces(iD,iD).Interior(:,3:4),'rows');
      Element=Element(not(Interior));
      ElementFaces=ElementFaces(not(Interior));
    end
    Faces(iD,iD).Exterior=[Faces(iD,iD).Exterior;Element,ElementFaces];
    for iElem=1:length(Element)
      Elements(iD).Faces.Exterior{Element(iElem),1}(ElementFaces(iElem))=1;
      Elements(iD).Faces.Boundary{Element(iElem),1}(ElementFaces(iElem))=iBoundary;
    end
    for iBN=1:numel(BoundaryNames)
      if ismember(iBoundary,Boundaries(iD).(BoundaryNames{iBN}))
        Faces(iD,iD).(BoundaryNames{iBN})=...
          [Faces(iD,iD).(BoundaryNames{iBN});Element,ElementFaces];
        for iElem=1:length(Element)
          Elements(iD).Faces.(BoundaryNames{iBN}){Element(iElem),1}(ElementFaces(iElem))=1;
        end
        BCs(iD).(BoundaryNames{iBN})=[BCs(iD).(BoundaryNames{iBN});Nodes];
      end
    end
  end
  Faces(iD,iD).Exterior=sortrows(Faces(iD,iD).Exterior);
  for iBN=1:numel(BoundaryNames)
    if isempty(Faces(iD,iD).(BoundaryNames{iBN}))
      Faces(iD,iD).(BoundaryNames{iBN})=zeros(0,2);
    else
      [~,OrderFaces]=sort(Faces(iD,iD).(BoundaryNames{iBN})(:,1));
      Faces(iD,iD).(BoundaryNames{iBN})=Faces(iD,iD).(BoundaryNames{iBN})(OrderFaces,:);
      BCs(iD).(BoundaryNames{iBN})=unique(BCs(iD).(BoundaryNames{iBN}));
    end
  end

  % Faces of the regions defined in the input deck
  if not(isempty(Regions))
    RegionNames=fieldnames(Regions(iD));
    for iRN=1:numel(RegionNames)
      Faces(iD,iD).(RegionNames{iRN})=[];
      Elements(iD).Faces.(RegionNames{iRN})=Flags;
      Coord=num2cell(Mesh(iD).Degree1.Nodes(1:Sizes(iD).NumSpaceDim,:),2);
      Nodes=find(Regions(iD).(RegionNames{iRN}).Definition(Coord{:}))';
      ElementNodes=ismember(Mesh(iD).Degree1.Elements',Nodes);
      Element=find(sum(ElementNodes,2)==NumFaceNodes);
      Keep=false(numel(Element),1);
      for iElem=1:numel(Element)
        Center=num2cell(mean(Mesh(iD).Degree1.Nodes(1:Sizes(iD).NumSpaceDim,...
          Mesh(iD).Degree1.Elements(:,Element(iElem))),2));
        Keep(iElem)=Regions(iD).(RegionNames{iRN}).Condition(Center{:});
      end
      Element=Element(Keep);
      A=(ElementNodes(Element,:)*diag(1:NumElementNodes))'; A(A==0)=[];
      A=reshape(A,NumFaceNodes,length(Element))';
      [~,ElementFaces]=ismember(A,sort(FaceNodesElem,2),'rows');
      Faces(iD,iD).(RegionNames{iRN})=[Element,ElementFaces];
      for iElem=1:length(Element)
        Elements(iD).Faces.(RegionNames{iRN}){Element(iElem),1}(ElementFaces(iElem))=1;
      end
      BCs(iD).(RegionNames{iRN})=Nodes;
    end
  end

  % Faces matching across the interface with the other discretizations
  if matchField(Boundaries(iD),'Interface')
    for iFace=1:size(Faces(iD,iD).Interface,1)
      iElem=Faces(iD,iD).Interface(iFace,1);
      iElemFace=Faces(iD,iD).Interface(iFace,2);
      NodesCoord=Mesh(iD).Degree1.Nodes(:,...
        Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFace,:),iElem))';
      for iDC=setdiff(1:Simulation.NumDiscretizations,iD)
        FaceNodesElemC=RefElement(iDC,iDC).Degree1.FaceNodesElem;
        Tolerance=Mesh(iDC).Degree1.MinElementSize/1000;
        [~,NodesC]=ismembertol(NodesCoord,Mesh(iDC).Degree1.Nodes',Tolerance,'ByRows',true);
        ElementNodesC=ismember(Mesh(iDC).Degree1.Elements',NodesC);
        iElemC=find(sum(ElementNodesC,2)==Sizes(iDC).NumSpaceDim);
        A=(ElementNodesC(iElemC,:)*diag(1:size(ElementNodesC,2)));
        A(A==0)=[];
        [~,iElemFaceC]=ismember(A,sort(FaceNodesElemC,2),'rows');
        FaceNodesC=Mesh(iDC).Degree1.Elements(FaceNodesElemC(iElemFaceC,:),iElemC)';
        Node2Match1stNode1=find(FaceNodesC==NodesC(1));
        Faces(iD,iDC).Interface(iFace,:)=...
          [iElem,iElemFace,iElemC,iElemFaceC,Node2Match1stNode1];
        Elements(iD).Faces.Interface{iElem,1}(2:3,iElemFace)=[iElemFaceC;Node2Match1stNode1];
      end
    end
  end

  % Faces matching across the periodic boundaries
  if matchField(Boundaries(iD),'PeriodicMaster')
    Faces(iD,iD).Periodic=zeros(size(Faces(iD,iD).PeriodicMaster,1),5);
    Elements(iD).Faces.Periodic=FlagsPaired;
    Tolerance=Mesh(iD).Degree1.MinElementSize/1000;
    for iFaceM=1:size(Faces(iD,iD).PeriodicMaster,1)
      iElemM=Faces(iD,iD).PeriodicMaster(iFaceM,1);
      iElemFaceM=Faces(iD,iD).PeriodicMaster(iFaceM,2);
      NodesM=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFaceM,:),iElemM);
      NodesCoordM=Mesh(iD).Degree1.Nodes(:,NodesM)';
      [~,Dim]=min(max(NodesCoordM)-min(NodesCoordM));
      InPlane=setdiff(1:Sizes(iD).NumSpaceDim,Dim);
      for iFaceS=1:size(Faces(iD,iD).PeriodicSlave,1)
        iElemS=Faces(iD,iD).PeriodicSlave(iFaceS,1);
        iElemFaceS=Faces(iD,iD).PeriodicSlave(iFaceS,2);
        NodesS=Mesh(iD).Degree1.Elements(FaceNodesElem(iElemFaceS,:),iElemS);
        [NodesMatch,NodesOrder]=ismembertol(Mesh(iD).Degree1.Nodes(InPlane,NodesS)',...
          NodesCoordM(:,InPlane),Tolerance,'ByRows',true);
        if all(NodesMatch) && length(unique(NodesOrder))==length(NodesMatch)
          break
        end
      end
      Faces(iD,iD).Periodic(iFaceM,:)=...
        [iElemM,iElemFaceM,iElemS,iElemFaceS,NodesOrder(1)];
      Elements(iD).Faces.Periodic{iElemM,1}(:,iElemFaceM)=[1;0];
      Elements(iD).Faces.Periodic{iElemS,1}(:,iElemFaceS)=[1;NodesOrder(1)];
    end
    if any(Faces(iD,iD).Periodic(:,5)==0)
      error('Not all periodic faces are properly coupled!');
    end
  end

  % Faces connectivity, with the periodic faces sharing the same face number
  Faces(iD,iD).GlobalLocal=cell(NumElements,1);
  if HDG
    Faces(iD,iD).Connectivity=zeros(NumElements,NumElementFaces);
    for iFace=1:size(Faces(iD,iD).Interior,1)
      Int=Faces(iD,iD).Interior(iFace,:);
      Faces(iD,iD).Connectivity(Int(1),Int(2))=iFace;
      Faces(iD,iD).Connectivity(Int(3),Int(4))=iFace;
    end
    for iFace=1:size(Faces(iD,iD).Exterior,1)
      Ext=Faces(iD,iD).Exterior(iFace,:);
      Faces(iD,iD).Connectivity(Ext(1),Ext(2))=iFace+size(Faces(iD,iD).Interior,1);
    end
    if matchField(Faces(iD,iD),'Periodic')
      for iFace=1:size(Faces(iD,iD).Periodic,1)
        Per=Faces(iD,iD).Periodic(iFace,:);
        Faces(iD,iD).Connectivity(Per(3),Per(4))=Faces(iD,iD).Connectivity(Per(1),Per(2));
      end
      Faces(iD,iD).Connectivity=changem(Faces(iD,iD).Connectivity,...
        1:length(unique(Faces(iD,iD).Connectivity)),unique(Faces(iD,iD).Connectivity));
    end

    % Map from the global face degrees of freedom to the local ones
    Faces(iD,iD).GlobalLocal=mat2cell(reshape(...
      ((repmat(reshape(Faces(iD,iD).Connectivity',[],1),1,NumGlobalComp*NumFaceNodesHigh)-1)...
      *NumGlobalComp*NumFaceNodesHigh+(1:NumGlobalComp*NumFaceNodesHigh))',[],1),...
      repmat(NumGlobalComp*NumElementFaces*NumFaceNodesHigh,NumElements,1));
  end
end

% Turn the flags of the element faces into a structure array
for iD=1:Simulation.NumDiscretizations
  Elements(iD).Faces=table2struct(struct2table(Elements(iD).Faces));
end

end
