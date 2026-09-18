function [Block]=getPeriodicDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes)
% Get the master and slave DOFs of the CG discretizations with periodic boundary conditions

for iD=1:Simulation.NumDiscretizations
  Block(iD,iD).DOFsMaster=[];
  Block(iD,iD).DOFsSlave=[];
  if not(strcmp(Parameters(iD).DiscretizationType,'CG')) || ...
     not(matchField(Boundaries(iD),'PeriodicMaster')) || isempty(Boundaries(iD).PeriodicMaster) || ...
     not(matchField(Boundaries(iD),'PeriodicSlave')) || isempty(Boundaries(iD).PeriodicSlave)
    continue
  end

  % Get master and slave nodes
  NodesM=getBoundaryNodes(Mesh(iD),RefElement(iD,iD),Sizes(iD),Boundaries(iD).PeriodicMaster);
  NodesS=getBoundaryNodes(Mesh(iD),RefElement(iD,iD),Sizes(iD),Boundaries(iD).PeriodicSlave);

  % Match slave nodes to master nodes along the direction in which the master nodes are aligned
  nsd=Sizes(iD).NumSpaceDim;
  Tolerance=Mesh(iD).Degree1.MinElementSize/1000;
  NodesCoordMxyz=Mesh(iD).Nodes(:,NodesM)';
  [Min,Dim]=min(max(NodesCoordMxyz(:,1:nsd))-min(NodesCoordMxyz(:,1:nsd)));
  if Min>Tolerance
    error('Periodic nodes aligned only in one direction are implemented for now!');
  end
  NodesCoordM=NodesCoordMxyz(:,setdiff(1:nsd,Dim));
  NodesCoordS=Mesh(iD).Nodes(setdiff(1:nsd,Dim),NodesS)';
  [NodesMatch,NodesOrder]=ismembertol(NodesCoordM,NodesCoordS,Tolerance,'ByRows',true);
  if not(all(NodesMatch))
    error('Not all periodic nodes are properly coupled!');
  end
  NodesS=NodesS(NodesOrder);

  % Get master and slave DOFs
  Block(iD,iD).DOFsMaster=reshape(repmat(NodesM,Sizes(iD).NumGlobalComp,1)'+...
    (0:Sizes(iD).NumGlobalComp-1)*Sizes(iD).NumNodes,[],1);
  Block(iD,iD).DOFsSlave=reshape(repmat(NodesS,Sizes(iD).NumGlobalComp,1)'+...
    (0:Sizes(iD).NumGlobalComp-1)*Sizes(iD).NumNodes,[],1);
end

end
