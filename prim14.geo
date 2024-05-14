// Gmsh project created on Mon Feb 26 17:05:30 2024
p = 0.0625;
p1 = 0.00625;
//+
Point(1) = {1, 1, 0, p};
//+
Point(2) = {1, 0, 0, p};
//+
Point(3) = {0, 1, 0, p};
//+
Point(4) = {0.1, 0, 0, p1};
//+
Point(5) = {0, 0.05, 0, p1};
//+
Point(6) = {-0.05, -0.175, 0, p1};

//+
Line(1) = {5, 3};
//+
Line(2) = {3, 1};
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 4};

//+
Circle(5) = {5, 6, 4};
//+
Line Loop(6) = {1, 2, 3, 4, -5};
//+
Plane Surface(7) = {6};
//+
Physical Line(8) = {1};
//+
Physical Line(9) = {2};
//+
Physical Line(10) = {3};
//+
Physical Line(11) = {4};
//+
Physical Line(12) = {5};
//+
Physical Surface(13) = {7};