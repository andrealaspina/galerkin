close all; clear; clc;

% Import mesh (fluid)
[Mesh(1)]=importGmshMesh('fsi_3d_channel_fluid');

% Import mesh (ALE)
[Mesh(2)]=importGmshMesh('fsi_3d_channel_mesh');

% Import mesh (structure)
[Mesh(3)]=importGmshMesh('fsi_3d_channel_structure');

% Save mesh
save('Mesh_fsi_3d_channel.mat','Mesh');