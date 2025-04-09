close all; clear; clc;

% Import mesh (fluid)
[Mesh(1)]=importGmshMesh('fsi_3d_tube_fluid');

% Import mesh (ALE)
[Mesh(2)]=importGmshMesh('fsi_3d_tube_mesh');

% Import mesh (structure)
[Mesh(3)]=importGmshMesh('fsi_3d_tube_structure');

% Save mesh
save('Mesh_fsi_3d_tube.mat','Mesh');