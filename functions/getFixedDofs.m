function [Block]=getFixedDofs(Block,Simulation,Mesh,Boundaries,Parameters,RefElement,Sizes)
% Get the DOFs of the CG discretizations that are strongly fixed on the boundary

for iD=1:Simulation.NumDiscretizations
  Block(iD,iD).DOFsFixed=[];
  if strcmp(Parameters(iD).DiscretizationType,'CG') && ...
     matchField(Boundaries(iD),'Fixed') && not(isempty(Boundaries(iD).Fixed))
    NodesFixed=getBoundaryNodes(Mesh(iD),RefElement(iD,iD),Sizes(iD),Boundaries(iD).Fixed);
    Block(iD,iD).DOFsFixed=reshape(repmat(NodesFixed,Sizes(iD).NumGlobalComp,1)'+...
      (0:Sizes(iD).NumGlobalComp-1)*Sizes(iD).NumNodes,[],1);
  end
end

end
