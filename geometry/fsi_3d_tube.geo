// Inputs
r = 1;
kF = 4;
kM = 4;
kS = 4;
L = 5.0;
R = 0.5;
T = 0.1;
o = 3; // (High-order options: 1=optimisation, 2=elastic+optimisation, 3=elastic)
a = 5; // (Meshing algorithm: 1=MeshAdapt, 2=Automatic, 3=InitialMeshOnly, 4=Automatic, 5=Delaunay, 6=Frontal-Delaunay)

// mesh fluid settings
Mesh.ElementOrder      = kF;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;

// fluid
Point(1) = {0, 0, 0, 1};
Point(2) = {L, 0, 0, 1};
Point(3) = {L, R, 0, 1};
Point(4) = {0, R, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Curve {1,  3} = 2^(r+1)+1;
Transfinite Curve {2, -4} = 2^(r-1)+1;
Transfinite Surface{1};
Mesh 1;
Mesh 2;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{1}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{21}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{38}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{55}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;

// save mesh fluid
Save "fsi_3d_tube_fluid.m";

// mesh ALE settings
Mesh.ElementOrder      = kM;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh 1;
Mesh 2;
Mesh 3;

// save mesh ALE
Save "fsi_3d_tube_mesh.m";

// delete mesh fluid
Delete{Volume{:}; Surface{:};}

// mesh structure settings
Mesh.ElementOrder      = kS;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;

// structure
Point(35) = {L, R+T, 0, 1};
Point(36) = {0, R+T, 0, 1};
Line(68) = { 3, 35};
Line(69) = {35, 36};
Line(70) = {36,  4};
Line Loop(2) = {-69, -70, 3, -68};
Plane Surface(72) = {2};
Transfinite Curve {69    } = 2^(r+1)+1;
Transfinite Curve {68, 70} = 1+1;
Transfinite Surface{72};
Mesh 1;
Mesh 2;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{72}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{94}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{116}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {Surface{138}; Layers{2^(r-1)};}
Mesh 1;
Mesh 2;
Mesh 3;

// save mesh structure
Save "fsi_3d_tube_structure.m";