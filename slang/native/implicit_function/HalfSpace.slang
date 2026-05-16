// HalfSpace — half-space implicit function, Slang realisation.
//
// A half-space is the set of points on one side of a plane:
//   f(p) = dot(p − plane_point, normal)
//   f(p) > 0  on the side the normal points toward  (inside the half-space)
//   f(p) = 0  on the bounding plane
//   f(p) < 0  on the opposite side
//
// Math:
//   f(p)  = dot(p − plane_point, normal)   (linear — gradient is constant)
//   ∇f(p) = normal                         (everywhere, not just on the plane)
//
// C++ reference: js_iteration_2/implicit_function/half_plane.hpp
//   The C++ version includes a transformation matrix (rotation/translation) and
//   a sign flip in eval_gradient based on inside/outside. This Slang port is
//   the clean, transformation-free version: the normal direction IS the gradient.
//   Add a Transformed<T> wrapper for rotation/translation when needed.
//
// Usage examples:
//   makeZPositiveHalfSpace()  →  f(p) = p.z  (positive above z=0)
//   makeZNegativeHalfSpace()  →  f(p) = -p.z (positive below z=0)
//   makeHalfSpace(normal, pt) →  f(p) = dot(p-pt, normalize(normal))

module HalfSpace;

__exported import IImplicitFunction;  // re-exports gradient<T>, BoundingBox

public struct HalfSpace : IImplicitFunction {
    public float3 normal;       // outward unit normal (caller must normalize)
    public float3 plane_point;  // any point on the bounding plane

    // f(p) = dot(p - plane_point, normal)
    [Differentiable]
    [NoDiffThis]
    public float eval(float3 p) {
        return dot(p - plane_point, normal);
    }

    // ∇f = normal  (constant — gradient of a linear function is its coefficient vector).
    [BackwardDerivativeOf(eval)]
    [NoDiffThis]
    public void eval_bwd(inout DifferentialPair<float3> dp, float dResult) {
        dp = DifferentialPair<float3>(dp.p, dp.d + dResult * normal);
    }

    public BoundingBox getBoundingBox() {
        // Half-space is unbounded; callers that need a bounded domain
        // should intersect with a finite bounding box externally.
        const float INF = 1e9f;
        BoundingBox bbox;
        bbox.minCorner = float3(-INF, -INF, -INF);
        bbox.maxCorner = float3( INF,  INF,  INF);
        return bbox;
    }
}

// General constructor — normalises the provided normal direction.
public HalfSpace makeHalfSpace(float3 normal, float3 plane_point) {
    HalfSpace h;
    h.normal      = normalize(normal);
    h.plane_point = plane_point;
    return h;
}

// f(p) = p.z   — the half-space above z=0.
// Used to cut the lower hemisphere: hemisphere = sphere ∩ {p.z ≥ 0}.
public HalfSpace makeZPositiveHalfSpace() {
    HalfSpace h;
    h.normal      = float3(0.0f, 0.0f, 1.0f);
    h.plane_point = float3(0.0f, 0.0f, 0.0f);
    return h;
}

// f(p) = -p.z  — the half-space below z=0.
public HalfSpace makeZNegativeHalfSpace() {
    HalfSpace h;
    h.normal      = float3(0.0f, 0.0f, -1.0f);
    h.plane_point = float3(0.0f, 0.0f, 0.0f);
    return h;
}
