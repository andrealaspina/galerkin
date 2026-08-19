function Block=getFixedDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes)
% Get the degrees of freedom fixed on the boundaries of the CG discretizations

for iD=1:Simulation.NumDiscretizations
  Block(iD,iD).DOFsFixed=[];
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
     matchField(Boundaries(iD),'Fixed') && not(isempty(Boundaries(iD).Fixed))

    % Reference element data
    FaceNodesElem=RefElement(iD,iD).FaceNodesElem;
    FaceNodesElem1=RefElement(iD,iD).Degree1.FaceNodesElem;
    NumElementNodes1=size(Mesh(iD).Degree1.Elements,1);
    NumFaceNodes1=size(FaceNodesElem1,2);

    % Nodes lying on the fixed boundaries
    NodesFixed=[];
    for iFace=Boundaries(iD).Fixed
      if Sizes(iD).NumSpaceDim==2
        Nodes1=(findNodes(Mesh(iD).Degree1,'Region','Edge',iFace))';
      elseif Sizes(iD).NumSpaceDim==3
        Nodes1=(findNodes(Mesh(iD).Degree1,'Region','Face',iFace))';
      end
      ElementNodes1=ismember(Mesh(iD).Degree1.Elements',Nodes1);
      Elements=find(sum(ElementNodes1,2)==NumFaceNodes1);
      A=(ElementNodes1(Elements,:)*diag(1:NumElementNodes1))'; A(A==0)=[];
      A=reshape(A,NumFaceNodes1,length(Elements))';
      [~,ElementFaces]=ismember(A,sort(FaceNodesElem1,2),'rows');
      for iElem=1:length(Elements)
        NodesFixed=[NodesFixed,...
          Mesh(iD).Elements(FaceNodesElem(ElementFaces(iElem),:),Elements(iElem))']; %#ok
      end
    end
    NodesFixed=unique(NodesFixed);

    % Degrees of freedom of the fixed nodes, for every component
    Block(iD,iD).DOFsFixed=reshape(repmat(NodesFixed,Sizes(iD).NumGlobalComp,1)'...
      +(0:Sizes(iD).NumGlobalComp-1)*Sizes(iD).NumNodes,[],1);
  end
end

end
