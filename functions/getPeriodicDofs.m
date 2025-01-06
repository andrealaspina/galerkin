function [Block]=getPeriodicDofs(...
         Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes)
         % Get periodic DOFs

% Initialize
for iD=1:Simulation.NumDiscretizations
  Block(iD,iD).DOFsMaster=[];
  Block(iD,iD).DOFsSlave=[];
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
   (matchField(Boundaries(iD),'PeriodicSlave') && not(isempty(Boundaries(iD).PeriodicSlave))) || ...
   (matchField(Boundaries(iD),'PeriodicMaster') && not(isempty(Boundaries(iD).PeriodicMaster)))
    % Get info
    FaceNodesElem=RefElement(iD,iD).FaceNodesElem;
    FaceNodesElemDeg1=RefElement(iD,iD).Degree1.FaceNodesElem;
    NumElemNodesDeg1=size(Mesh(iD).Degree1.Elements,1);
    NumFaceNodesDeg1=size(FaceNodesElemDeg1,2);
    NumGlobalComp=Sizes(iD).NumGlobalComp;
    NumNodes=Sizes(iD).NumNodes;
  end
  
  % Master nodes
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
     matchField(Boundaries(iD),'PeriodicMaster') && not(isempty(Boundaries(iD).PeriodicMaster))
    NodesM=[];
    for iFacePeriodicMaster=Boundaries(iD).PeriodicMaster
      if Sizes(iD).NumSpaceDim==2
        NodesDeg1=(findNodes(Mesh(iD).Degree1,'Region','Edge',iFacePeriodicMaster))';
      elseif Sizes(iD).NumSpaceDim==3
        NodesDeg1=(findNodes(Mesh(iD).Degree1,'Region','Face',iFacePeriodicMaster))';
      end
      ElementNodesDeg1=ismember(Mesh(iD).Degree1.Elements',NodesDeg1);
      Elements=find(sum(ElementNodesDeg1,2)==NumFaceNodesDeg1);
      A=(ElementNodesDeg1(Elements,:)*diag(1:NumElemNodesDeg1))'; A(A==0)=[];
      A=reshape(A,size(FaceNodesElemDeg1,2),length(Elements))'; B=sort(FaceNodesElemDeg1,2);
      [~,ElementsFaces]=ismember(A,B,'rows');
      for iElem=1:length(Elements)
        Element=Elements(iElem);
        ElementNodes=FaceNodesElem(ElementsFaces(iElem),:);
        NodesM=[NodesM,Mesh(iD).Elements(ElementNodes,Element)']; %#ok
      end
    end
    NodesM=unique(NodesM);
  end

  % Slave nodes
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
     matchField(Boundaries(iD),'PeriodicSlave') && not(isempty(Boundaries(iD).PeriodicSlave))
    NodesS=[];
    for iFacePeriodicSlave=Boundaries(iD).PeriodicSlave
      if Sizes(iD).NumSpaceDim==2
        NodesDeg1=(findNodes(Mesh(iD).Degree1,'Region','Edge',iFacePeriodicSlave))';
      elseif Sizes(iD).NumSpaceDim==3
        NodesDeg1=(findNodes(Mesh(iD).Degree1,'Region','Face',iFacePeriodicSlave))';
      end
      ElementNodesDeg1=ismember(Mesh(iD).Degree1.Elements',NodesDeg1);
      Elements=find(sum(ElementNodesDeg1,2)==NumFaceNodesDeg1);
      A=(ElementNodesDeg1(Elements,:)*diag(1:NumElemNodesDeg1))'; A(A==0)=[];
      A=reshape(A,size(FaceNodesElemDeg1,2),length(Elements))'; B=sort(FaceNodesElemDeg1,2);
      [~,ElementsFaces]=ismember(A,B,'rows');
      for iElem=1:length(Elements)
        Element=Elements(iElem);
        ElementNodes=FaceNodesElem(ElementsFaces(iElem),:);
        NodesS=[NodesS,Mesh(iD).Elements(ElementNodes,Element)']; %#ok
      end
    end
    NodesS=unique(NodesS);
  end

  % Get coupled master and slave nodes
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
     matchField(Boundaries(iD),'PeriodicSlave') && not(isempty(Boundaries(iD).PeriodicSlave)) && ...
     matchField(Boundaries(iD),'PeriodicMaster') && not(isempty(Boundaries(iD).PeriodicMaster))

    % Match slave nodes to master nodes
    Tolerance=Mesh(iD).Degree1.MinElementSize/1000;
    NodesCoordMxyz=Mesh(iD).Nodes(:,NodesM)';
    [Min,Dim]=min(max(NodesCoordMxyz(:,1:Sizes(iD).NumSpaceDim))...
                 -min(NodesCoordMxyz(:,1:Sizes(iD).NumSpaceDim)));
    if Min>Tolerance
      error('Periodic nodes aligned only in one direction are implemented for now!');
    end
    NodesCoordM=NodesCoordMxyz(:,setdiff(1:Sizes(iD).NumSpaceDim,Dim));
    NodesCoordS=Mesh(iD).Nodes(setdiff(1:Sizes(iD).NumSpaceDim,Dim),NodesS)';
    [NodesMatch,NodesOrder]=ismembertol(NodesCoordM,NodesCoordS,Tolerance,'ByRows',true);
    if not(all(NodesMatch))
      error('Not all periodic nodes are properly coupled!');
    end
    NodesS=NodesS(NodesOrder);
    
    % Get master DOFs
    Block(iD,iD).DOFsMaster=reshape(repmat(NodesM,NumGlobalComp,1)'...
                                                              +((0:NumGlobalComp-1)*NumNodes),[],1);

    % Get slave DOFs
    Block(iD,iD).DOFsSlave=reshape(repmat(NodesS,NumGlobalComp,1)'...
                                                              +((0:NumGlobalComp-1)*NumNodes),[],1);
  end
  
end