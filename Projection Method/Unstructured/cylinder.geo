//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {3, 0, 0, 1.0};
//+
Point(3) = {3, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0.5, 0.5, 0, 1.0};
//+
Point(6) = {0.6, 0.5, 0, 1.0};
//+
Point(7) = {0.5, 0.6, 0, 1.0};
//+
Point(8) = {0.4, 0.5, 0, 1.0};
//+
Point(9) = {0.5, 0.4, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 8};
//+
Circle(7) = {8, 5, 9};
//+
Circle(8) = {9, 5, 6};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("INLET", 9) = {4};
//+
Physical Curve("OUTLET", 10) = {2};
//+
Physical Curve("WALL", 11) = {1, 3};
//+
Physical Curve("BODY", 12) = {5, 6, 7, 8};
