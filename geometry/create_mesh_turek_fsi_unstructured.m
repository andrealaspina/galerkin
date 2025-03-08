close all; clear; clc;

%% Fluid mesh

clearvars -except MeshF MeshM MeshS

% Import mesh
[Mesh]=importGmshMesh('/Users/andrealaspina/Desktop/galerkin/geometry/turek_fsi_unstructured_fluid');

% Remove node at (0.2,0.2)
[~,NodeRemove]=ismember([0.2,0.2],Mesh.Nodes','Rows');
Mesh.Nodes=Mesh.Nodes(:,[1:NodeRemove-1,NodeRemove+1:end]);
Mesh.Elements(Mesh.Elements>NodeRemove)=Mesh.Elements(Mesh.Elements>NodeRemove)-1;

% Re-number nodes ----------------------------------------------------------------------------------
nsd=2;

X=Mesh.Nodes';
C=Mesh.Elements';

C1=C(:,1:(nsd+1));
Ck=C(:,(nsd+2):end);

Nmax1=length(unique(C1));

nAbove=unique(C1(C1(1:end)>Nmax1));
nBelow=unique(Ck(Ck(1:end)<=Nmax1));

for n=1:length(nAbove)
  C1(C1==nAbove(n))=nBelow(n);
  Ck(Ck==nBelow(n))=nAbove(n);
  X([nBelow(n),nAbove(n)],:)=X([nAbove(n),nBelow(n)],:);
end
C=[C1,Ck];

Mesh.Nodes=X';
Mesh.Elements=C';
% --------------------------------------------------------------------------------------------------

% Store mesh
MeshF=Mesh;

%% ALE mesh

clearvars -except MeshF MeshM MeshS

% Import mesh
[Mesh]=importGmshMesh('/Users/andrealaspina/Desktop/galerkin/geometry/turek_fsi_unstructured_mesh');

% Remove node at (0.2,0.2)
[~,NodeRemove]=ismember([0.2,0.2],Mesh.Nodes','Rows');
Mesh.Nodes=Mesh.Nodes(:,[1:NodeRemove-1,NodeRemove+1:end]);
Mesh.Elements(Mesh.Elements>NodeRemove)=Mesh.Elements(Mesh.Elements>NodeRemove)-1;

% Re-number nodes ----------------------------------------------------------------------------------
nsd=2;

X=Mesh.Nodes';
C=Mesh.Elements';

C1=C(:,1:(nsd+1));
Ck=C(:,(nsd+2):end);

Nmax1=length(unique(C1));

nAbove=unique(C1(C1(1:end)>Nmax1));
nBelow=unique(Ck(Ck(1:end)<=Nmax1));

for n=1:length(nAbove)
  C1(C1==nAbove(n))=nBelow(n);
  Ck(Ck==nBelow(n))=nAbove(n);
  X([nBelow(n),nAbove(n)],:)=X([nAbove(n),nBelow(n)],:);
end
C=[C1,Ck];

Mesh.Nodes=X';
Mesh.Elements=C';
% --------------------------------------------------------------------------------------------------

% Store mesh
MeshM=Mesh;

%% Structure mesh

clearvars -except MeshF MeshM MeshS

% Import mesh
[Mesh]=importGmshMesh('/Users/andrealaspina/Desktop/galerkin/geometry/turek_fsi_unstructured_structure');

% Re-number nodes ----------------------------------------------------------------------------------
nsd=2;

X=Mesh.Nodes';
C=Mesh.Elements';

C1=C(:,1:(nsd+1));
Ck=C(:,(nsd+2):end);

Nmax1=length(unique(C1));

nAbove=unique(C1(C1(1:end)>Nmax1));
nBelow=unique(Ck(Ck(1:end)<=Nmax1));

for n=1:length(nAbove)
  C1(C1==nAbove(n))=nBelow(n);
  Ck(Ck==nBelow(n))=nAbove(n);
  X([nBelow(n),nAbove(n)],:)=X([nAbove(n),nBelow(n)],:);
end
C=[C1,Ck];

Mesh.Nodes=X';
Mesh.Elements=C';
% --------------------------------------------------------------------------------------------------

% Store mesh
MeshS=Mesh;

%% Save meshes

% Merge meshes
clearvars -except MeshF MeshM MeshS
Mesh(1)=MeshF;
Mesh(2)=MeshM;
Mesh(3)=MeshS;
clear MeshF MeshM MeshS

% Plot mesh
figure
hold on;
if size(Mesh(1).Elements,1)<=6
  pdeplot(Mesh(1).Nodes,Mesh(1).Elements,'EdgeColor','b');
end
plot(Mesh(1).Nodes(1,:),Mesh(1).Nodes(2,:),'ko');
if size(Mesh(2).Elements,1)<=6
  pdeplot(Mesh(2).Nodes,Mesh(2).Elements,'EdgeColor','c');
end
plot(Mesh(2).Nodes(1,:),Mesh(2).Nodes(2,:),'k^');
if size(Mesh(3).Elements,1)<=6
  pdeplot(Mesh(3).Nodes,Mesh(3).Elements,'EdgeColor','r');
end
plot(Mesh(3).Nodes(1,:),Mesh(3).Nodes(2,:),'kx');
axis equal
xlim([0,2.50]);
ylim([0,0.41]);

% Plot geometry
for iD=1:3
  Geometry(iD)=geometryFromMesh(createpde(),Mesh(iD).Nodes,Mesh(iD).Elements(1:3,:));
  figure
  pdegplot(Geometry(iD),'EdgeLabels','on');
end
axis equal
xlim([0,2.50]);
ylim([0,0.41]);

% Save mesh
save('/Users/andrealaspina/Desktop/galerkin/geometry/Mesh_turek_fsi_unstructured.mat','Mesh');
fprintf('\nMesh_turek_fsi_unstructured.mat completed');