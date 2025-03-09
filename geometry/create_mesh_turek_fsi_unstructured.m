close all; clear; clc;

% Import mesh (fluid)
[Mesh(1)]=importGmshMesh('/Users/andrealaspina/Desktop/galerkin/geometry/turek_fsi_unstructured_fluid');

% Import mesh (ALE)
[Mesh(2)]=importGmshMesh('/Users/andrealaspina/Desktop/galerkin/geometry/turek_fsi_unstructured_mesh');

% Import mesh (structure)
[Mesh(3)]=importGmshMesh('/Users/andrealaspina/Desktop/galerkin/geometry/turek_fsi_unstructured_structure');

% Save mesh
save('/Users/andrealaspina/Desktop/galerkin/geometry/Mesh_turek_fsi_unstructured.mat','Mesh');
fprintf('\nMesh_turek_fsi_unstructured.mat completed');