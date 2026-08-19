function MeshHigh=generateMeshHighOrder(MeshLow,RefElementGeom,RefElementApprox)
% Generate a hybrid mesh with a degree for the geometry and a higher one for the approximation

% Vertex nodes, shared by the two degrees
Cg=MeshLow.Elements';
Xg=MeshLow.Nodes';
NumElements=size(MeshLow.Elements,2);
NumVertexNodes=size(MeshLow.Nodes,1)+1;
NumElementNodes=size(RefElementApprox.NodesCoordElem,1);
C1=Cg(:,1:NumVertexNodes);
X1=Xg(1:length(unique(C1)),:);

% Coordinates of the nodes added inside each element
Nga=RefElementGeom.ShapeFunctionsElem;
pinvNaa=RefElementApprox.PseudoinverseShapeFunctionsElem;
XAdd=cell(NumElements,1);
parfor iElem=1:NumElements
  Xge=Xg(Cg(iElem,:),:); %#ok
  Xae=pinvNaa*(Nga*Xge);
  XAdd{iElem}=Xae(NumVertexNodes+1:end,:);
end
XAdd=vertcat(XAdd{:});

% Condense the added nodes shared by neighboring elements
[~,iUnique,iShared]=uniquetol(double(XAdd),1e-9,'ByRows',true);
CAdd=reshape(size(X1,1)+iShared,[NumElementNodes-NumVertexNodes,NumElements])';

% Merge the vertex nodes with the added ones
MeshHigh.Nodes=[X1;XAdd(iUnique,:)]';
MeshHigh.Elements=[C1,CAdd]';

end
