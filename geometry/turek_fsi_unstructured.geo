//Inputs
r = 1;
kF = 4;
kM = 2;
kS = 4;
L = 2.50;
H = 0.41;
R = 0.05;
l = 0.35;
h = 0.02;
X = 0.20;
Y = 0.20;
o = 3; // (High-order options: 1=optimisation, 2=elastic+optimisation, 3=elastic)
a = 5; // (Meshing algorithm: 1=MeshAdapt, 2=Automatic, 3=InitialMeshOnly, 4=Automatic, 5=Delaunay, 6=Frontal-Delaunay)

// structure
Point(1) = {X+R*Cos(Asin(h/(2*R))), Y-h/2, 0, 1};
Point(2) = {X+R+l                 , Y-h/2, 0, 1};
Point(3) = {X+R+l                 , Y+h/2, 0, 1};
Point(4) = {X+R*Cos(Asin(h/(2*R))), Y+h/2, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(1) = {1, 2, 3, 4};
Transfinite Curve {1, 3} = 2^(r+1)+1;
Transfinite Curve {2, 4} = r+1;
Plane Surface(1) = {1};
Transfinite Surface {1} Alternate;

// mesh structure
Mesh.ElementOrder      = kS;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh 1;
Mesh 2;

// save mesh structure
Save "turek_fsi_unstructured_structure.m";

// delete structure
Delete{Surface{1}; Curve{4};}

// cylinder
Point(5) = {X  , Y  , 0, 1};
Point(6) = {X-R, Y  , 0, 1};
Point(7) = {X  , Y-R, 0, 1};
Point(8) = {X  , Y+R, 0, 1};
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 1};
Circle(7) = {4, 5, 8};
Circle(8) = {8, 5, 6};
Line Loop(2) = {5, 6, 1, 2, 3, 7, 8};
Transfinite Curve {5, 6, 7, 8} = 2^r+1;

// channel
Point( 9) = {0, 0, 0, 1};
Point(10) = {L, 0, 0, 1};
Point(11) = {L, H, 0, 1};
Point(12) = {0, H, 0, 1};
Line( 9) = { 9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12,  9};
Line Loop(3) = {9, 10, 11, 12};
Transfinite Curve {9, -11} = 2^(r+2)+1 Using Progression 2^(2^(-(r+2)));
Transfinite Curve {10}     = Ceil(3*2^(r-2))+1;
Transfinite Curve {12}     = 2^r+1;
Plane Surface(2) = {3, 2};

// mesh fluid
Mesh.ElementOrder      = kF;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh 1;
Mesh 2;

// save mesh fluid
Save "turek_fsi_unstructured_fluid.m";

// mesh ALE
Mesh.ElementOrder      = kM;
Mesh.HighOrderOptimize = o;
Mesh.Algorithm         = a;
Mesh 1;
Mesh 2;

// save mesh ALE
Save "turek_fsi_unstructured_mesh.m";

// exit
Exit;
