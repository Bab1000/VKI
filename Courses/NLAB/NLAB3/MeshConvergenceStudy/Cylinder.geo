// Mesh for a cylinder in a pipe using OpenCASCADE
SetFactory("OpenCASCADE");

// =========================================================
// GEOMETRY GENERATION
// =========================================================

// Define characteristic lengths for meshing
lc = 0.05;  // Fine mesh size
lcc = 0.07; // Coarser mesh size

// Define geometric dimensions
L = 2.2;  // Length of the domain
H = 0.41; // Height of the domain
D = 0.1;  // Diameter of the cylinder

// Define circle center coordinates
Xc = 0.2; // X-coordinate of the circle center
Yc = 0.2; // Y-coordinate of the circle center
R = D/2;  // Radius of the cylinder

// =========== Define Points ===========
Point(1) = {0, 0, 0, lc};   // Bottom-left corner
Point(2) = {L, 0, 0, lcc};  // Bottom-right corner
Point(3) = {0, H, 0, lc};   // Top-left corner
Point(4) = {L, H, 0, lcc};  // Top-right corner

// =========== Define Lines ===========
Line(1) = {1, 2}; // Bottom edge
Line(2) = {2, 4}; // Right edge
Line(3) = {4, 3}; // Top edge
Line(4) = {3, 1}; // Left edge

// =========== Define Circle ===========
Circle(5) = {Xc, Yc, 0, R, 0, 2*Pi}; // Full circle centered at (Xc, Yc)

// Force structured mesh on the cylinder boundary
Transfinite Curve{5} = 100; // Enforces 100 divisions along the circle

// =========== Define Surfaces ===========
// Outer rectangular boundary
Curve Loop(3) = {3, 4, 1, 2};       

// Inner circular boundary
Curve Loop(4) = {5};       

// Subtract the cylinder from the rectangle
Plane Surface(2) = {3, 4};

// =========================================================
// MESH GENERATION
// =========================================================

// Define physical groups for boundary conditions
Physical Line("Inlet") = {4};    // Left boundary
Physical Line("Outlet") = {2};   // Right boundary
Physical Line("Bottom") = {1};   // Bottom boundary
Physical Line("Top") = {3};      // Top boundary
Physical Line("Cylinder") = {5}; // Cylinder boundary (entire hole)

// =========================================================
// Boundary Layer Definition
// =========================================================

// Wall boundary layer for bottom and top boundaries
Field[1] = BoundaryLayer;
Field[1].EdgesList = {1, 3};  // Bottom and Top boundaries
Field[1].Size = 0.001;
Field[1].thickness = 0.01;
Field[1].ratio = 1.05;
Field[1].IntersectMetrics = 1;
Field[1].Quads = 1;
Field[1].PointsList = {1, 2, 3, 4};  // Points for the surface (optional based on your geometry)
BoundaryLayer Field = 1;

// Cylinder boundary layer
Field[2] = BoundaryLayer;
Field[2].EdgesList = {5};  // Cylinder boundary
Field[2].Size = 0.0001;
Field[2].thickness = 0.01;
Field[2].ratio = 1.05;
Field[2].IntersectMetrics = 1;
Field[2].Quads = 1;
BoundaryLayer Field = 2;

// =========================================================
// MESH SETTINGS
// =========================================================

// Use Delaunay algorithm (works well with OpenCASCADE)
Mesh.Algorithm = 6;
