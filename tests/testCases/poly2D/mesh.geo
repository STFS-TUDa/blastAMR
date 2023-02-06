  Point(1) = {0, 0, 0, 0.005};
  Point(2) = {0.1, 0, 0, 0.005};
  Point(3) = {0.1, 0.1, 0, 0.005};
  Point(4) = {0, 0.1, 0, 0.005};
  Line(1) = {1, 2};
  Line(2) = {2, 3};
  Line(3) = {3, 4};
  Line(4) = {4, 1};
  Line Loop(6) = {4, 1, 2, 3};
  Plane Surface(6) = {6};
  Physical Volume("internal") = {1};
  Extrude {0, 0, 0.01} {
   Surface{6};
   Layers{1};
   Recombine;
  }
  Physical Surface("front") = {28};
  Physical Surface("back") = {6};
  Physical Surface("bottom") = {27};
  Physical Surface("left") = {15};
  Physical Surface("top") = {19};
  Physical Surface("right") = {23};
