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
//
//   - radius and center are shape parameters, not AD variables.
//     UnitSphere does NOT implement IDifferentiable; its fields are constants
//     in the AD graph. If you later need gradient-based shape fitting
//     (e.g. optimising sphere parameters to match a point cloud), implement
//     IDifferentiable and remove the no_diff qualifiers.
//
//   - gradient() is not a method — call gradient(sphere, p) from
//     IImplicitFunction. The [BackwardDerivativeOf(eval)] below is what
//     Slang's AD system uses internally.
//
// AUTODIFF DESIGN:
//   eval() is [Differentiable]. eval_bwd() is the hand-coded backward pass
//   registered via [BackwardDerivativeOf(eval)]. This replaces the auto-
//   generated backward pass — same mathematical result, no symbolic overhead.
//
//   Cross-check: auto-diff of `r² - dot(d,d)` where `d = p - center` gives
//   ∂f/∂p = -2*(p - center), which matches eval_bwd. To verify, temporarily
//   remove [BackwardDerivativeOf(eval)] and compare gradient() output.
//
//   VERIFY WITH SLANGC: [BackwardDerivativeOf] signature for a struct method
//   where the struct is NOT IDifferentiable — `this` should not appear as a
//   DifferentialPair. Adjust if slangc rejects this form.

module UnitSphere;

__exported import IImplicitFunction;   // re-exports gradient<T> and BoundingBox to callers

public struct UnitSphere : IImplicitFunction {
    public float  radius;   // shape parameter — not an AD variable (struct is not IDifferentiable)
    public float3 center;   // shape parameter — not an AD variable

    // f(p) = r² − dot(p − center, p − center)
    // [NoDiffThis]: UnitSphere is not IDifferentiable; radius/center have no derivatives.
    // Suppresses compiler warnings about non-differentiable struct members in a
    // [Differentiable] method. Use [Differentiable] struct + remove [NoDiffThis] if
    // you later need gradient-based shape fitting (optimising radius/center).
    [Differentiable]
    [NoDiffThis]
    public float eval(float3 p) {
        float3 d = p - center;
        return radius * radius - dot(d, d);
    }

    // Hand-coded backward pass for eval().
    // ∂f/∂p = −2·(p − center)
    // Slang uses this instead of the auto-generated pass when
    // bwd_diff(_evalForDiff<UnitSphere>) is called.
    [BackwardDerivativeOf(eval)]
    [NoDiffThis]
    public void eval_bwd(inout DifferentialPair<float3> dp, float dResult) {
        float3 grad_p = dResult * (-2.0f * (dp.p - center));
        dp = DifferentialPair<float3>(dp.p, dp.d + grad_p);
    }

    public BoundingBox getBoundingBox() {
        float margin = radius * 1.1f;
        BoundingBox bbox;
        bbox.minCorner = center - float3(margin, margin, margin);
        bbox.maxCorner = center + float3(margin, margin, margin);
        return bbox;
    }
};

// Convenience constructors — Slang structs do not support constructor overloads.
public UnitSphere makeUnitSphere() {
    UnitSphere s;
    s.radius = 1.0f;
    s.center = float3(0.0f, 0.0f, 0.0f);
    return s;
}

public UnitSphere makeSphere(float radius) {
    UnitSphere s;
    s.radius = radius;
    s.center = float3(0.0f, 0.0f, 0.0f);
    return s;
}

public UnitSphere makeSphere(float radius, float3 center) {
    UnitSphere s;
    s.radius = radius;
    s.center = center;
    return s;
}
