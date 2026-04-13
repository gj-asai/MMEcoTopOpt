L = 50;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {L, 0, 0, 1.0};
Point(3) = {2*L, 0, 0, 1.0};
Point(4) = {2*L, L, 0, 1.0};
Point(5) = {L, L, 0, 1.0};
Point(6) = {L, 2*L, 0, 1.0};
Point(7) = {0, 2*L, 0, 1.0};
Point(8) = {0, L, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 8};
Line(4) = {8, 1};
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 5};
Line(8) = {5, 6};
Line(9) = {6, 7};
Line(10) = {7, 8};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, -2};
Curve Loop(3) = {-3, 8, 9, 10};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Transfinite Surface {1, 2, 5, 8};
Transfinite Surface {2, 3, 4, 5};
Transfinite Surface {8, 5, 6, 7};
Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} = 51;
Recombine Surface {1, 2, 3};
