close all; clear; clc;

% Import mesh (fluid)
[Mesh(1)]=importGmshMesh('turek_fsi_unstructured_fluid');

% Import mesh (ALE)
[Mesh(2)]=importGmshMesh('turek_fsi_unstructured_mesh');

% Import mesh (structure)
[Mesh(3)]=importGmshMesh('turek_fsi_unstructured_structure');

% Save mesh
save('Mesh_turek_fsi_unstructured.mat','Mesh');