//Inputs
r = 1;
kF = 1;
kM = 1;
kS = 1;
L = 3.00;
W = 1.00;
H = 0.50;
l = 0.05;
w = 0.60;
h = 0.40;
d = 0.50;
o = 3; // (High-order options: 1=optimisation, 2=elastic+optimisation, 3=elastic)
a = 5; // (Meshing algorithm: 1=MeshAdapt, 2=Automatic, 3=InitialMeshOnly, 4=Automatic, 5=Delaunay, 6=Frontal-Delaunay)

// set factory
SetFactory("OpenCASCADE");

// structure
Box(1) = {d, -w/2, -H/2, l, w, h};
Transfinite Curve {1:8} = 2^(r-1)+1;
Transfinite Curve {9:12} = 1+1;
Transfinite Surface {1:6} Left;

// mesh structure
Mesh.ElementOrder      = kS;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh 1;
Mesh 2;
Mesh 3;

// save mesh structure
Save "fsi_3d_channel_structure.mat";

// fluid
Box(2) = {0, -W/2, -H/2, L, W, H};
BooleanDifference{Volume{2}; Delete;}{Volume{1}; Delete;}
Transfinite Curve {13, 15, 18, 23} = 2^(r-1)+1;
Transfinite Curve {14, 16, 21, 24} = 2^r+1;
Transfinite Curve {17, 19, 20, 22} = 3*2^r+1;
Transfinite Surface {7, 8, 9, 10, 12} Left;

// mesh fluid
Mesh.ElementOrder      = kF;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh.MeshSizeMin       = 1/2^r;
Mesh.MeshSizeMax       = 1/2^r;
Mesh.MeshSizeFactor    = 1.5;
Mesh 1;
Mesh 2;
Mesh 3;

// save mesh fluid
Save "fsi_3d_channel_fluid.m";

// mesh ALE
Mesh.ElementOrder      = kM;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh.MeshSizeMin       = 1/2^r;
Mesh.MeshSizeMax       = 1/2^r;
Mesh.MeshSizeFactor    = 1.5;
Mesh 1;
Mesh 2;
Mesh 3;

// save mesh ALE
Save "fsi_3d_channel_mesh.m";

// exit
Exit;