// Foundation types for ImpliSolid — Slang realisation.
//
// Sign convention (F-REP / R-function, matching the C++ and Python references):
//   f(p) > 0  →  p is INSIDE the solid
//   f(p) = 0  →  p is ON the surface
//   f(p) < 0  →  p is OUTSIDE the solid
//
// This is NOT signed-distance-function (SDF) convention — gradient magnitude != 1.
// Matches neural implicit surface literature (NeRF, occupancy networks, DeepSDF).
//
// Note on Slang modules: struct fields default to `internal` access.
// All fields used outside this module must be explicitly `public`.
//
// References:
//   js_iteration_2/implicit_function/implicit_function.hpp
//   python_implicit/primitives.py

module types;

// Axis-aligned bounding box in absolute world-space coordinates.
public struct BoundingBox {
    public float3 minCorner;
    public float3 maxCorner;
};
