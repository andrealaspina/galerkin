function [Nodes]=getBoundaryNodes(Mesh,RefElement,Sizes,BoundaryFaces)
% Get the (high-order) mesh nodes lying on the given geometry edges (2D) or faces (3D)

% Get info
FaceNodesElem=RefElement.FaceNodesElem;
FaceNodesElemDeg1=RefElement.Degree1.FaceNodesElem;
NumElemNodesDeg1=size(Mesh.Degree1.Elements,1);
NumFaceNodesDeg1=size(FaceNodesElemDeg1,2);
if Sizes.NumSpaceDim==2
  Region='Edge';
elseif Sizes.NumSpaceDim==3
  Region='Face';
end

% Collect the nodes of the element faces lying on each boundary face
Nodes=[];
for iBoundaryFace=BoundaryFaces
  NodesDeg1=findNodes(Mesh.Degree1,'Region',Region,iBoundaryFace)';
  ElementNodesDeg1=ismember(Mesh.Degree1.Elements',NodesDeg1);
  Elements=find(sum(ElementNodesDeg1,2)==NumFaceNodesDeg1);
  A=(ElementNodesDeg1(Elements,:)*diag(1:NumElemNodesDeg1))'; A(A==0)=[];
  A=reshape(A,NumFaceNodesDeg1,length(Elements))'; B=sort(FaceNodesElemDeg1,2);
  [~,ElementsFaces]=ismember(A,B,'rows');
  for iElem=1:length(Elements)
    Nodes=[Nodes,Mesh.Elements(FaceNodesElem(ElementsFaces(iElem),:),Elements(iElem))']; %#ok
  end
end
Nodes=unique(Nodes);

end
