function [MeshHigh]=generateMeshHighOrder(MeshLow,RefElementGeometry,RefElementApproxim)
% Generate a high-order mesh from a low-order one by adding the nodes of the approximation

% Get general parameters
Cg=MeshLow.Elements';
Xg=MeshLow.Nodes';
nsd=size(MeshLow.Nodes,1);
NumElements=size(MeshLow.Elements,2);
NumElementNodes1=nsd+1;
NumElementNodesApp=size(RefElementApproxim.NodesCoordElem,1);
NumElementNodesAdd=NumElementNodesApp-NumElementNodes1;

% Linear mesh
C1=Cg(:,1:NumElementNodes1);
N1=length(unique(C1));
X1=Xg(1:N1,:);

% Get shape functions
Nga=RefElementGeometry.ShapeFunctionsElem;
pinvNaa=RefElementApproxim.PseudoinverseShapeFunctionsElem;

% Compute coordinates of the additional nodes
Xa_add=cell(NumElements,1);
parfor iElem=1:NumElements
  Cge=Cg(iElem,:);
  Xge=Xg(Cge,:); %#ok
  Xae=pinvNaa*(Nga*Xge);
  Xa_add{iElem}=Xae(NumElementNodes1+1:end,:);
end
Xa_add=vertcat(Xa_add{:});

% Condense shared nodes
tol=1e-9;
[~,ind1,ind2]=uniquetol(double(Xa_add),tol,'ByRows',true);
Xa_add_unique=Xa_add(ind1,:);
Ca_add_unique=reshape(N1+ind2,[NumElementNodesAdd,NumElements])';

% Append the additional nodes to the linear mesh
MeshHigh.Nodes=[X1;Xa_add_unique]';
MeshHigh.Elements=[C1,Ca_add_unique]';

end
