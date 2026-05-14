// Foundation types for ImpliSolid Slang realisation.
//
// Sign convention (F-REP / R-function, matching the C++ reference):
//   f(p) > 0  →  p is INSIDE the solid
//   f(p) = 0  →  p is ON the surface
//   f(p) < 0  →  p is OUTSIDE the solid
//
// This is NOT signed-distance-function (SDF) convention.
// Gradients are not unit-length; normalize explicitly when needed for rendering normals.
//
// Reference: js_iteration_2/implicit_function/implicit_function.hpp

module types;

// Axis-aligned bounding box.
// Absolute world-space coordinates (not relative to any center).
struct BoundingBox {
    float3 minCorner;
    float3 maxCorner;
};
