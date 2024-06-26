cl__1 = 1000000.0;
Point(0) = {475000.0, 0.0, 0.0, cl__1};
Point(1) = {475000.0, 0.0, 600000.0, cl__1};
Point(2) = {525000.0, 0.0, 0.0, cl__1};
Point(3) = {525000.0, 0.0, 600000.0, cl__1};
Point(4) = {475000.0, -10000.0, 0.0, cl__1};
Point(5) = {475000.0, -10000.0, 600000.0, cl__1};
Point(6) = {475000.0, -80000.0, 0.0, cl__1};
Point(7) = {475000.0, -80000.0, 600000.0, cl__1};
Point(8) = {543686.9375, -25000.0, 0.0, cl__1};
Point(9) = {543686.9375, -25000.0, 600000.0, cl__1};
Point(10) = {529949.5625, -100000.0, 0.0, cl__1};
Point(11) = {529949.5625, -100000.0, 600000.0, cl__1};
Point(12) = {749747.75, -100000.0, 0.0, cl__1};
Point(13) = {749747.75, -100000.0, 600000.0, cl__1};
Point(14) = {593686.9375, -25000.0, 0.0, cl__1};
Point(15) = {593686.9375, -25000.0, 600000.0, cl__1};
Point(16) = {621161.6875, -35000.0, 0.0, cl__1};
Point(17) = {621161.6875, -35000.0, 600000.0, cl__1};
Point(18) = {854697.3125, -120000.0, 0.0, cl__1};
Point(19) = {854697.3125, -120000.0, 600000.0, cl__1};
Point(20) = {804697.3125, -120000.0, 0.0, cl__1};
Point(21) = {804697.3125, -120000.0, 600000.0, cl__1};
Point(22) = {0.0, 0.0, 0.0, cl__1};
Point(23) = {0.0, 0.0, 600000.0, cl__1};
Point(24) = {0.0, -10000.0, 0.0, cl__1};
Point(25) = {0.0, -10000.0, 600000.0, cl__1};
Point(26) = {0.0, -80000.0, 0.0, cl__1};
Point(27) = {0.0, -80000.0, 600000.0, cl__1};
Point(28) = {0.0, -300000.0, 0.0, cl__1};
Point(29) = {0.0, -300000.0, 600000.0, cl__1};
Point(30) = {1000000.0, 0.0, 0.0, cl__1};
Point(31) = {1000000.0, 0.0, 600000.0, cl__1};
Point(32) = {1000000.0, -25000.0, 0.0, cl__1};
Point(33) = {1000000.0, -25000.0, 600000.0, cl__1};
Point(34) = {1000000.0, -35000.0, 0.0, cl__1};
Point(35) = {1000000.0, -35000.0, 600000.0, cl__1};
Point(36) = {1000000.0, -120000.0, 0.0, cl__1};
Point(37) = {1000000.0, -120000.0, 600000.0, cl__1};
Point(38) = {1000000.0, -300000.0, 0.0, cl__1};
Point(39) = {1000000.0, -300000.0, 600000.0, cl__1};
//+
Line(1) = {23, 23};
//+
Line(2) = {23, 25};
//+
Line(3) = {25, 27};
//+
Line(4) = {27, 29};
//+
Line(5) = {23, 1};
//+
Line(6) = {1, 3};
//+
Line(7) = {3, 31};
//+
Line(8) = {31, 33};
//+
Line(9) = {33, 35};
//+
Line(10) = {35, 39};
//+
Line(11) = {39, 29};
//+
Line(12) = {25, 5};
//+
Line(13) = {27, 7};
//+
Line(14) = {1, 9};
//+
Line(15) = {3, 15};
//+
Line(16) = {15, 17};
//+
Line(17) = {15, 33};
//+
Line(18) = {35, 17};
//+
Line(19) = {37, 19};
//+
Line(20) = {19, 21};
//+
Line(21) = {19, 17};
//+
Line(22) = {9, 13};
//+
Line(23) = {13, 21};
//+
Line(24) = {5, 9};
//+
Line(25) = {7, 11};
//+
Line(26) = {11, 13};
//+
Line(27) = {39, 38};
//+
Line(28) = {31, 30};
//+
Line(29) = {30, 32};
//+
Line(30) = {32, 33};
//+
Line(31) = {35, 34};
//+
Line(32) = {34, 32};
//+
Line(33) = {34, 38};
//+
Delete {
  Curve{10}; 
}
//+
Line(34) = {35, 37};
//+
Line(35) = {37, 39};
//+
Curve Loop(1) = {4, -11, -35, 19, 20, -23, -26, -25, -13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 25, 26, -22, -24, -12, 3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, 24, -14, -5, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {6, 15, 16, -21, 20, -23, -22, -14};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {19, 21, -18, 34};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {18, -16, 17, 9};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {17, -8, -7, 15};
//+
Plane Surface(7) = {7};
//+
Line(36) = {37, 36};
//+
Delete {
  Curve{33}; 
}
//+
Line(37) = {34, 36};
//+
Line(38) = {36, 38};
//+
Curve Loop(8) = {38, -27, -35, 36};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {37, -36, -34, 31};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {31, 32, 30, 9};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {29, 30, -8, 28};
//+
Plane Surface(11) = {11};
//+
Line(39) = {38, 28};
//+
Line(40) = {28, 26};
//+
Line(41) = {26, 24};
//+
Line(42) = {24, 22};
//+
Line(43) = {22, 0};
//+
Line(44) = {4, 24};
//+
Line(45) = {26, 6};
//+
Line(46) = {6, 10};
//+
Line(47) = {4, 8};
//+
Line(48) = {0, 8};
//+
Line(49) = {0, 2};
//+
Line(50) = {2, 14};
//+
Line(51) = {14, 16};
//+
Line(52) = {16, 18};
//+
Line(53) = {18, 36};
//+
Line(54) = {18, 20};
//+
Line(55) = {20, 12};
//+
Line(56) = {12, 8};
//+
Line(57) = {10, 12};
//+
Line(58) = {2, 30};
//+
Line(59) = {32, 14};
//+
Line(60) = {34, 16};
//+
Curve Loop(12) = {39, 40, 45, 46, 57, -55, -54, 53, 38};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {41, -44, 47, -56, -57, -46, -45};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {52, 53, -37, 60};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {60, -51, -59, -32};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {50, -59, -29, -58};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {56, -48, 49, 50, 51, 52, 54, 55};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {44, 42, 43, 48, -47};
//+
Plane Surface(18) = {18};
//+
Line(61) = {28, 29};
//+
Line(62) = {26, 27};
//+
Line(63) = {24, 25};
//+
Line(64) = {22, 23};
//+
Curve Loop(19) = {61, -4, -62, -40};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {62, -3, -63, -41};
//+
Plane Surface(20) = {20};

//+
Curve Loop(21) = {64, 2, -63, 42};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {11, -61, -39, -27};
//+
Plane Surface(22) = {22};
//+
Line(65) = {6, 7};
//+
Line(66) = {11, 10};
//+
Curve Loop(23) = {62, 13, -65, -45};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {65, 25, 66, -46};
//+
Plane Surface(24) = {24};
//+
Line(67) = {12, 13};
//+
Line(68) = {21, 20};
//+
Curve Loop(25) = {26, -67, -57, -66};
//+
Surface(25) = {25};
//+
Curve Loop(26) = {23, 68, 55, 67};
//+
Surface(26) = {26};
//+
Line(69) = {19, 18};
//+
Curve Loop(27) = {68, -54, -69, 20};
//+
Plane Surface(27) = {27};
//+
Curve Loop(28) = {69, 53, -36, 19};
//+
Plane Surface(28) = {28};
//+
Surface Loop(1) = {12, 22, 1, 19, 8, 23, 24, 25, 26, 27, 28};
//+
Volume(1) = {1};
//+
Line(70) = {5, 4};
//+
Line(71) = {9, 8};
//+
Curve Loop(29) = {12, 70, 44, 63};
//+
Plane Surface(29) = {29};
//+
Curve Loop(30) = {70, 47, -71, -24};
//+
Plane Surface(30) = {30};
//+
Curve Loop(31) = {22, -67, 56, -71};
//+
Plane Surface(31) = {31};
//+
Surface Loop(2) = {29, 30, 31, 13, 20, 2, 23, 24, 25};
//+
Volume(2) = {2};
//+
Line(72) = {1, 0};
//+
Curve Loop(32) = {14, 71, -48, -72};
//+
Plane Surface(32) = {32};
//+
Curve Loop(33) = {72, -43, 64, 5};
//+
Plane Surface(33) = {33};
//+
Surface Loop(3) = {33, 32, 18, 21, 3, 29, 30};
//+
Volume(3) = {3};
//+
Line(73) = {17, 16};
//+
Curve Loop(34) = {21, 73, 52, -69};
//+
Plane Surface(34) = {34};
//+
Curve Loop(35) = {18, 73, -60, -31};
//+
Plane Surface(35) = {35};
//+
Surface Loop(4) = {35, 34, 5, 9, 14, 28};
//+
Volume(4) = {4};
//+
Line(74) = {15, 14};
//+
Curve Loop(36) = {16, 73, -51, -74};
//+
Plane Surface(36) = {36};
//+
Curve Loop(37) = {74, -59, 30, -17};
//+
Plane Surface(37) = {37};
//+
Surface Loop(5) = {6, 10, 15, 36, 37, 35};
//+
Volume(5) = {5};
//+
Line(75) = {3, 2};
//+
Curve Loop(38) = {15, 74, -50, -75};
//+
Plane Surface(38) = {38};
//+
Curve Loop(39) = {7, 28, -58, -75};
//+
Plane Surface(39) = {39};
//+
Curve Loop(40) = {6, 75, -49, -72};
//+
Plane Surface(40) = {40};
//+
Surface Loop(6) = {7, 11, 16, 39, 38, 37};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {40, 4, 17, 38, 36, 34, 31, 32, 27, 26};
//+
Volume(7) = {7};
//+
Physical Volume("Asthenosphere", 76) = {1};
//+
Physical Volume("Oceanic lithosphere", 77) = {2};
//+
Physical Volume("Oceanic crust", 78) = {3};
//+
Physical Volume("Continental lithosphere", 79) = {4};
//+
Physical Volume("Weak zone", 80) = {7};
//+
Physical Volume("Continental lower crust", 81) = {5};
//+
Physical Volume("Continental upper crust", 82) = {6};
//+
Physical Surface("Bottom", 83) = {22};
//+
Physical Surface("Asthenosphere", 84) = {19, 1, 8, 12};
//+
Physical Surface("Ocean dirichlet", 85) = {20, 21};
//+
Physical Surface("zmax", 86) = {2, 3, 4, 5, 7, 6};
//+
Physical Surface("Continent dirichlet", 87) = {9, 10, 11};
//+
Physical Surface("zmax", 86) += {13, 18, 17, 14, 16, 15};
