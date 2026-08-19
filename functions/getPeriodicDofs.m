function Block=getPeriodicDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes)
% Get the master and slave degrees of freedom of the periodic boundaries of the CG discretizations

Sides={'Master','Slave'};
for iD=1:Simulation.NumDiscretizations
  Block(iD,iD).DOFsMaster=[];
  Block(iD,iD).DOFsSlave=[];

  % Periodicity requires both a master and a slave boundary
  Periodic=strcmp(Parameters(iD).DiscretizationType,'CG');
  for iSide=1:2
    Periodic=Periodic && matchField(Boundaries(iD),['Periodic',Sides{iSide}]) && ...
             not(isempty(Boundaries(iD).(['Periodic',Sides{iSide}])));
  end
  if Periodic

    % Reference element data
    FaceNodesElem=RefElement(iD,iD).FaceNodesElem;
    FaceNodesElem1=RefElement(iD,iD).Degree1.FaceNodesElem;
    NumElementNodes1=size(Mesh(iD).Degree1.Elements,1);
    NumFaceNodes1=size(FaceNodesElem1,2);

    % Nodes lying on the master and on the slave boundaries
    Nodes=cell(1,2);
    for iSide=1:2
      for iFace=Boundaries(iD).(['Periodic',Sides{iSide}])
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
          Nodes{iSide}=[Nodes{iSide},...
            Mesh(iD).Elements(FaceNodesElem(ElementFaces(iElem),:),Elements(iElem))'];
        end
      end
      Nodes{iSide}=unique(Nodes{iSide});
    end

    % Match the slave nodes to the master ones along the direction of periodicity
    Tolerance=Mesh(iD).Degree1.MinElementSize/1000;
    NodesCoordMaster=Mesh(iD).Nodes(:,Nodes{1})';
    [Thickness,Dim]=min(max(NodesCoordMaster(:,1:Sizes(iD).NumSpaceDim))...
                       -min(NodesCoordMaster(:,1:Sizes(iD).NumSpaceDim)));
    if Thickness>Tolerance
      error('Periodic nodes aligned only in one direction are implemented for now!');
    end
    InPlane=setdiff(1:Sizes(iD).NumSpaceDim,Dim);
    [NodesMatch,NodesOrder]=ismembertol(NodesCoordMaster(:,InPlane),...
      Mesh(iD).Nodes(InPlane,Nodes{2})',Tolerance,'ByRows',true);
    if not(all(NodesMatch))
      error('Not all periodic nodes are properly coupled!');
    end
    Nodes{2}=Nodes{2}(NodesOrder);

    % Degrees of freedom of the master and of the slave nodes, for every component
    for iSide=1:2
      Block(iD,iD).(['DOFs',Sides{iSide}])=...
        reshape(repmat(Nodes{iSide},Sizes(iD).NumGlobalComp,1)'...
        +(0:Sizes(iD).NumGlobalComp-1)*Sizes(iD).NumNodes,[],1);
    end
  end
end

end
