function [Faces]=getFacesRegions(...
         Simulation,Mesh,Faces,Regions,RefElement,Sizes)
         % Get faces and BCs

% Initialize structures
for iD=1:Simulation.NumDiscretizations
  % Get info
  FaceNodesElem=RefElement(iD,iD).Degree1.FaceNodesElem;
  NumElemNodes=size(Mesh(iD).Degree1.Elements,1);
  NumFaceNodes=size(FaceNodesElem,2);
    
  % Get faces of specific regions
  if not(isempty(Regions))
    RegionsNames=fieldnames(Regions(iD));
    for iRegionName=1:length(RegionsNames)
      Faces(iD,iD).(RegionsNames{iRegionName})=[];
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
    end
  end
end

end