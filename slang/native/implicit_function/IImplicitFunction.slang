// IImplicitFunction: the core interface every implicit solid must satisfy.
//
// Slang realisation of mp5_implicit::implicit_function
// Reference: js_iteration_2/implicit_function/implicit_function.hpp
//
// Each implementing type represents one implicit solid (or CSG tree node).
// In the C++ reference, eval_implicit and eval_gradient operate on batches
// (vectorized_vect) for cache efficiency. Here each method works on a single
// point; GPU thread-level parallelism replaces the explicit batch loop.
//
// The two required methods correspond directly to:
//   C++ eval_implicit  →  eval(float3 p)
//   C++ eval_gradient  →  gradient(float3 p)

module IImplicitFunction;

import types;

interface IImplicitFunction {
    // Implicit function value at p.
    // f > 0: inside, f = 0: surface, f < 0: outside.
    float eval(float3 p);

    // Gradient ∇f at p. Not normalised — magnitude is not 1.
    // Direction points toward the solid interior (inward normal).
    // Normalise with normalize() when you need a surface normal for rendering.
    float3 gradient(float3 p);

    // Axis-aligned bounding box in world space.
    // The polygoniser uses this to set the grid extent.
    BoundingBox getBoundingBox();
};
