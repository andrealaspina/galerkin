function [Mesh]=importGmshMesh(...
         FileName)
         % Import Gmsh mesh

% Run .geo script
system(['/Applications/Gmsh.app/Contents/MacOS/gmsh ',FileName,'.geo -3']);

% Run exported .m script
run([FileName,'.m']);
fprintf('\n');

% Remove exported .m script
delete([FileName,'.m']);

% Number of spatial dimensions
if max(abs(msh.POS(:,3)))<1e-12
  nsd=2;
else
  nsd=3;
end

% Import nodes coordinates
X=msh.POS(:,1:nsd);

% Polynomial degree
if nsd==2
  if     matchField(msh,'TRIANGLES')
    k=1;
  elseif matchField(msh,'TRIANGLES6')
    k=2;
  elseif matchField(msh,'TRIANGLES10')
    k=3;
  elseif matchField(msh,'TRIANGLES15')
    k=4;
  elseif matchField(msh,'TRIANGLES21')
    k=5;
  elseif matchField(msh,'TRIANGLES28')
    k=6;
  elseif matchField(msh,'TRIANGLES36')
    k=7;
  elseif matchField(msh,'TRIANGLES45')
    k=8;
  end
elseif nsd==3
  if     matchField(msh,'TETS')
    k=1;
  elseif matchField(msh,'TETS10')
    k=2;
  elseif matchField(msh,'TETS20')
    k=3;
  elseif matchField(msh,'TETS35')
    k=4;
  elseif matchField(msh,'TETS56')
    k=5;
  elseif matchField(msh,'TETS84')
    k=6;
  end
end

% Import elements nodes with re-ordering
if nsd==2
  switch k
    case 1
      order=[1,2,3];
      C=msh.TRIANGLES(:,order);
    case 2
      order=[1,2,3,4,5,6];
      C=msh.TRIANGLES6(:,order);
    case 3
      order=[1,2,3,4,5,6,7,8,9,10];
      C=msh.TRIANGLES10(:,order);
    case 4
      order=[1,2,3,4,5,6,7,8,9,10,11,12,14,15,13];
      C=msh.TRIANGLES15(:,order);
    case 5
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,20,18,21,16,19];
      C=msh.TRIANGLES21(:,order);
    case 6
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,24,25,21,26,27,19,22,23,28];
      C=msh.TRIANGLES28(:,order);
    case 7
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,28,29,30,24,31,32,33,22,25,26,27,35,36,34];
      C=msh.TRIANGLES36(:,order);
    case 8
      order=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,32,33,34,35,27,36,37,38,39,25,28,29,30,31,41,44,42,45,40,43];
      C=msh.TRIANGLES45(:,order);
    otherwise
      error('Element not yet implemented');
  end
elseif nsd==3
  switch k
    case 1
      order=[2,3,1,4];
      C=msh.TETS(:,order);
    case 2
      order=[2,3,1,4,6,7,5,10,9,8];
      C=msh.TETS10(:,order);
    case 3
      order=[2,3,1,4,7,8,9,10,6,5,16,15,14,13,12,11,17,19,20,18];
      C=msh.TETS20(:,order);
    case 4
      order=[2,3,1,4,8,9,10,11,12,13,7,6,5,22,21,20,19,18,17,16,15,14,23,24,25,30,31,29,33,34,32,27,28,26,35];
      C=msh.TETS35(:,order);
    case 5
      order=[2,3,1,4,9,10,11,12,13,14,15,16,8,7,6,5,28,27,26,25,24,23,22,21,20,19,18,17,29,32,30,33,31,34,42,45,43,46,41,44,48,51,49,52,47,50,36,39,37,40,35,38,54,55,53,56];
      C=msh.TETS56(:,order);
    case 6
      order=[2,3,1,4,10,11,12,13,14,15,16,17,18,19,9,8,7,6,5,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,35,38,39,36,40,41,37,42,43,44,56,60,61,57,62,63,55,58,59,64,66,70,71,67,72,73,65,68,69,74,46,50,51,47,52,53,45,48,49,54,76,77,75,78,80,81,79,84,83,82];
      C=msh.TETS84(:,order);
    case 7
      error('Element not supported in Gmsh');
    case 8
      error('Element not supported in Gmsh');
    otherwise
      error('Element not yet implemented');
  end
end

% Remove unconnected nodes
UnconnectedNodes=sort(setdiff(1:size(X,1),unique(C(:))),'descend');
for iN=1:numel(UnconnectedNodes)
  NodeRemove=UnconnectedNodes(iN);
  X=X([1:NodeRemove-1,NodeRemove+1:end],:);
  C(C>NodeRemove)=C(C>NodeRemove)-1;
end

% Re-number nodes
C1=C(:,1:(nsd+1));
Ck=C(:,(nsd+2):end);
Nmax1=numel(unique(C1));
nAbove=unique(C1(C1>Nmax1));
nBelow=unique(Ck(Ck<=Nmax1));
for n=1:numel(nAbove)
  C1(C1==nAbove(n))=nBelow(n);
  Ck(Ck==nBelow(n))=nAbove(n);
  X([nBelow(n),nAbove(n)],:)=X([nAbove(n),nBelow(n)],:);
end
C=[C1,Ck];

% Store mesh
Mesh.Nodes=X';
Mesh.Elements=C';

% Get axis limits
MinMaxAux=minmax(double([Mesh.Nodes]));
XMinMaxAux=MinMaxAux(1,:)+[-1,+1]*(MinMaxAux(1,1)==MinMaxAux(1,2));
YMinMaxAux=MinMaxAux(2,:)+[-1,+1]*(MinMaxAux(2,1)==MinMaxAux(2,2));
if nsd==3
  ZMinMaxAux=MinMaxAux(3,:)+[-1,+1]*(MinMaxAux(3,1)==MinMaxAux(3,2));
else
  ZMinMaxAux=[-1,+1];
end

% Plot
figure('Color','w');

% Plot mesh
subplot(1,2,1);
if nsd==2
  pdeplot(Mesh.Nodes(1:nsd,:),Mesh.Elements(1:nsd+1,:),'EdgeColor','b');
elseif nsd==3
  pdeplot3D(Mesh.Nodes(1:nsd,:),Mesh.Elements(1:nsd+1,:),...
    'FaceColor','w','EdgeColor','b','FaceAlpha',1);
  ChildrenAux=get(gca,'Children'); delete(ChildrenAux([2,3]));
end
box on; axis equal;
title('Mesh');
set(gca,'Visible','on')
xlim(XMinMaxAux); xlabel('x'); ylim(YMinMaxAux); ylabel('y'); zlim(ZMinMaxAux); zlabel('z');
pause(eps);

% Plot geometry
subplot(1,2,2);
Geometry=geometryFromMesh(createpde(),...
  Mesh.Nodes(:,1:max(max(Mesh.Elements(1:nsd+1,:)))),Mesh.Elements(1:nsd+1,:));
if nsd==2
  Plots=pdegplot(Geometry,'EdgeLabels','on');
  set(Plots(1),'Color','b')
elseif nsd==3
  Plots=pdegplot(Geometry,'FaceLabels','on','FaceAlpha',0.25);
  set(Plots(1),'FaceColor','b')
  ChildrenAux=get(gca,'Children'); delete(ChildrenAux([4,5]));
end
grid on; box on; axis equal;
title('Geometry');
xlim(XMinMaxAux); xlabel('x'); ylim(YMinMaxAux); ylabel('y'); zlim(ZMinMaxAux); zlabel('z');
pause(eps);

end