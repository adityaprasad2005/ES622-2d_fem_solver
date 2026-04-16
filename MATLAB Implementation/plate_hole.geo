// --- Parameters (Matching feaforall.com tutorial) ---
L = 50.0;   // Half-Length of the plate (100mm / 2)
W = 25.0;   // Half-Width/Height of the plate (50mm / 2)
a = 20.0;   // Radius of the hole
lc = 0.5;   // Characteristic length (controls mesh density)

// --- Points ---
// Format: Point(id) = {x, y, z, mesh_size};
Point(1) = {a, 0, 0, lc/3};   // Bottom edge of hole (highly refined for stress concentration)
Point(2) = {L, 0, 0, lc};     // Bottom right corner
Point(3) = {L, W, 0, lc};     // Top right corner
Point(4) = {0, W, 0, lc};     // Top left corner
Point(5) = {0, a, 0, lc/3};   // Top edge of hole (highly refined for stress concentration)
Point(6) = {0, 0, 0, lc};     // Center of hole (used as the arc origin)

// --- Lines ---
Line(1) = {1, 2};             // Bottom symmetrical edge
Line(2) = {2, 3};             // Right edge (where tension is applied)
Line(3) = {3, 4};             // Top edge
Line(4) = {4, 5};             // Left symmetrical edge
Circle(5) = {5, 6, 1};        // Quarter circle arc for the hole

// --- Surface ---
Curve Loop(1) = {1, 2, 3, 4, 5}; 
Plane Surface(1) = {1};

// --- CRITICAL: Force 4-Node Quadrilateral Elements ---
// Gmsh defaults to triangles (TRI3/TRI6). The article mentions Quad elements
// converge better for this geometry, so we force Q4 mapping.
Recombine Surface {1}; 

// Make the mesh structurally smooth
Mesh.Algorithm = 8; // Frontal-Delaunay for Quads


Physical Surface("Plate") = {1};