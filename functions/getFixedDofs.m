function [Block]=getFixedDofs(...
         Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes)
         % Get fixed DOFs

% Initialize
for iD=1:Simulation.NumDiscretizations
  Block(iD,iD).DOFsFixed=[];
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
     matchField(Boundaries(iD),'Fixed') && not(isempty(Boundaries(iD).Fixed))
    % Get info
    FaceNodesElem=RefElement(iD,iD).FaceNodesElem;
    FaceNodesElemDeg1=RefElement(iD,iD).Degree1.FaceNodesElem;
    NumElemNodesDeg1=size(Mesh(iD).Degree1.Elements,1);
    NumFaceNodesDeg1=size(FaceNodesElemDeg1,2);
    NumGlobalComp=Sizes(iD).NumGlobalComp;
    NumNodes=Sizes(iD).NumNodes;
    
    % Get fixed nodes
    NodesFixed=[];
    for iFaceFixed=Boundaries(iD).Fixed
      if Sizes(iD).NumSpaceDim==2
        NodesDeg1=(findNodes(Mesh(iD).Degree1,'Region','Edge',iFaceFixed))';
      elseif Sizes(iD).NumSpaceDim==3
        NodesDeg1=(findNodes(Mesh(iD).Degree1,'Region','Face',iFaceFixed))';
      end
      ElementNodesDeg1=ismember(Mesh(iD).Degree1.Elements',NodesDeg1);
      Elements=find(sum(ElementNodesDeg1,2)==NumFaceNodesDeg1);
      A=(ElementNodesDeg1(Elements,:)*diag(1:NumElemNodesDeg1))'; A(A==0)=[];
      A=reshape(A,size(FaceNodesElemDeg1,2),length(Elements))'; B=sort(FaceNodesElemDeg1,2);
      [~,ElementsFaces]=ismember(A,B,'rows');
      for iElem=1:length(Elements)
        Element=Elements(iElem);
        ElementNodes=FaceNodesElem(ElementsFaces(iElem),:);
        NodesFixed=[NodesFixed,Mesh(iD).Elements(ElementNodes,Element)']; %#ok
      end
    end
    NodesFixed=unique(NodesFixed);
    
    % Get fixed DOFs
    Block(iD,iD).DOFsFixed=reshape(repmat(NodesFixed,NumGlobalComp,1)'...
                                                              +((0:NumGlobalComp-1)*NumNodes),[],1);
  end
end