// UnitSphere: sphere implicit function.
//
// Named after the C++ reference class (js_iteration_2/implicit_function/unit_sphere.hpp),
// though it accepts arbitrary radius and center — not restricted to r=1, origin.
//
// Math:
//   f(p)  = r² − ‖p − center‖²
//   ∇f(p) = −2·(p − center)
//
// The gradient points toward center (inward), consistent with f > 0 = inside.
// It is NOT unit-length: ‖∇f‖ = 2‖p − center‖. Normalize when needed.
//
// This struct has no transformation matrix. The C++ reference applies
// inv_transf_matrix before evaluating; add a TransformableSphere wrapper
// (or a generic Transformed<T: IImplicitFunction>) when that is needed.
//
// Cross-testing reference:
//   js_iteration_2/implicit_function/unit_sphere.hpp  eval_implicit / eval_gradient
//   python_implicit/primitives.py  UnitSphere.implicitFunction / implicitGradient

module UnitSphere;

import types;
import IImplicitFunction;

struct UnitSphere : IImplicitFunction {
    float radius;
    float3 center;

    // f(p) = r² − dot(p − center, p − center)
    float eval(float3 p) {
        float3 d = p - center;
        return radius * radius - dot(d, d);
    }

    // ∇f(p) = −2·(p − center)
    float3 gradient(float3 p) {
        return -2.0f * (p - center);
    }

    BoundingBox getBoundingBox() {
        float margin = radius * 1.1f;
        BoundingBox bbox;
        bbox.minCorner = center - float3(margin, margin, margin);
        bbox.maxCorner = center + float3(margin, margin, margin);
        return bbox;
    }
};

// Convenience constructor: unit sphere at origin.
UnitSphere makeUnitSphere() {
    UnitSphere s;
    s.radius = 1.0f;
    s.center = float3(0.0f, 0.0f, 0.0f);
    return s;
}

// Convenience constructor: sphere at origin with given radius.
UnitSphere makeSphere(float radius) {
    UnitSphere s;
    s.radius = radius;
    s.center = float3(0.0f, 0.0f, 0.0f);
    return s;
}

// Convenience constructor: sphere at given center and radius.
UnitSphere makeSphere(float radius, float3 center) {
    UnitSphere s;
    s.radius = radius;
    s.center = center;
    return s;
}
