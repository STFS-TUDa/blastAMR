Point(1) = {0,    -0.01, -0.0001, 0.001};
Point(2) = {0.02, -0.01, -0.0001, 0.001};
Point(3) = {0.02, 0.01,  -0.0001, 0.001};
Point(4) = {0,    0.01,  -0.0001, 0.001};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(5) = {5};
Physical Volume("internal") = {1};
Extrude {0, 0, 0.0002} {
 Surface{5};
 Layers{1};
 Recombine;
}
Physical Surface("front", 28) = {5};
Physical Surface("back", 29) = {27};
Physical Surface("inlet", 30) = {14};
Physical Surface("outlet", 31) = {22};
Physical Surface("sides", 32) = {26, 18};
