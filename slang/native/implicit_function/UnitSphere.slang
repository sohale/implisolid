// UnitSphere — sphere implicit function, Slang realisation.
//
// Named after the C++ reference class; accepts arbitrary radius and center.
//
// Math (identical to C++ and Python references):
//   f(p)  =  r² − ‖p − center‖²      (positive inside, zero on surface)
//   ∇f(p) =  −2·(p − center)         (inward, not unit-length)
//
// C++ reference:  js_iteration_2/implicit_function/unit_sphere.hpp
//                   eval_implicit  ↔  eval()
//                   eval_gradient  ↔  eval_bwd() via [BackwardDerivativeOf]
// Python reference: python_implicit/primitives.py  UnitSphere.implicitFunction
//
// KEY DIFFERENCES FROM THE C++ VERSION:
//
//   - No transformation matrix (transf_matrix / inv_transf_matrix).
//     Add a generic Transformed<T : IImplicitFunction> wrapper when needed.
//     The C++ version bakes the transform into the class; Slang's generics
//     make a separate wrapper cleaner and more composable.
//
//   - radius and center are `no_diff`. They are shape parameters, not AD
//     variables. The polygoniser differentiates f w.r.t. p (the query point),
//     not w.r.t. the shape parameters. Mark them `IDifferentiable` if you
//     later need gradient-based shape fitting (e.g. fitting a sphere to a
//     point cloud).
//
//   - gradient() is NOT a method here. Call gradient(sphere, p) from
//     IImplicitFunction.slang. The backward pass below is what Slang's AD
//     system uses when bwd_diff(sphere.eval) is called.
//
// AUTODIFF DESIGN:
//   eval() is [Differentiable]. eval_bwd() provides the hand-coded backward
//   pass via [BackwardDerivativeOf(eval)]. This replaces the auto-generated
//   backward pass with the analytically correct one — same result, no overhead.
//
//   Cross-check: the auto-generated backward pass of `r² - dot(d,d)` where
//   `d = p - center` should produce `d/dp = -2*(p - center)`, which matches
//   eval_bwd. Verify by temporarily removing [BackwardDerivativeOf(eval)] and
//   comparing outputs.

module UnitSphere;

import types;
import IImplicitFunction;

struct UnitSphere : IImplicitFunction {
    no_diff float  radius;   // shape parameter — not in the AD graph
    no_diff float3 center;   // shape parameter — not in the AD graph

    // f(p) = r² − dot(p − center, p − center)
    [Differentiable]
    float eval(float3 p) {
        float3 d = p - center;
        return radius * radius - dot(d, d);
    }

    // Hand-coded backward pass for eval().
    // ∂f/∂p = −2·(p − center)
    // Slang uses this instead of the auto-generated pass when bwd_diff(sphere.eval)
    // is called (via the free function gradient<T> in IImplicitFunction.slang).
    //
    // VERIFY WITH SLANGC: the exact signature for [BackwardDerivativeOf] on a
    // struct method depends on whether UnitSphere implements IDifferentiable.
    // Because UnitSphere does NOT implement IDifferentiable (radius and center are
    // no_diff), `this` does not appear as a DifferentialPair — it is accessed
    // as a plain UnitSphere inside the body. If slangc requires an explicit
    // `this` parameter or a different form, adjust accordingly.
    [BackwardDerivativeOf(eval)]
    void eval_bwd(inout DifferentialPair<float3> dp, float dResult) {
        float3 grad_p = dResult * (-2.0f * (dp.p - center));
        dp = DifferentialPair<float3>(dp.p, dp.d + grad_p);
    }

    BoundingBox getBoundingBox() {
        float margin = radius * 1.1f;
        BoundingBox bbox;
        bbox.minCorner = center - float3(margin, margin, margin);
        bbox.maxCorner = center + float3(margin, margin, margin);
        return bbox;
    }
};

// ---------------------------------------------------------------------------
// Convenience constructors (Slang structs do not support constructor overloads
// the way C++ classes do; use free functions instead).
// ---------------------------------------------------------------------------

UnitSphere makeUnitSphere() {
    UnitSphere s;
    s.radius = 1.0f;
    s.center = float3(0.0f, 0.0f, 0.0f);
    return s;
}

UnitSphere makeSphere(float radius) {
    UnitSphere s;
    s.radius = radius;
    s.center = float3(0.0f, 0.0f, 0.0f);
    return s;
}

UnitSphere makeSphere(float radius, float3 center) {
    UnitSphere s;
    s.radius = radius;
    s.center = center;
    return s;
}
