Point(1) = {0, 0, 0, 1.0};
Point(2) = {100, 0, 0, 1.0};
Point(3) = {100, 40, 0, 1.0};
Point(4) = {0, 40, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};

Transfinite Surface {1};
Transfinite Curve {1, 3} = 301 Using Progression 1;
Transfinite Curve {2, 4} = 121 Using Progression 1;
Recombine Surface {1};
