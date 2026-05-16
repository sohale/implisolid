// CrispSubtract — CSG subtraction, Slang realisation.
//
// Computes A minus B using the F-REP formula:
//   f(p) = min(f_A(p), −f_B(p))
//        = intersection of A with the complement of B
//
// f(p) > 0  inside A and outside B  (the subtracted region)
// f(p) ≤ 0  elsewhere
//
// Gradient — piecewise, follows the active constraint:
//   f_A(p) < −f_B(p)  →  ∇f = ∇f_A          (A surface is tighter)
//   f_A(p) > −f_B(p)  →  ∇f = −∇f_B         (complement-of-B surface is tighter)
//   f_A(p) = −f_B(p)  →  on the sharp edge; we choose ∇f_A (arbitrary, follows A branch)
//
// C++ reference: crisp_subtract.hpp
//   eval_implicit: f = min(f1, -f2)                        ← same
//   eval_gradient: pre-negate grad2, then select by condition  ← same math, different code style
//
// WHY manual [BackwardDerivativeOf]:
//   min() is not differentiable at f_A = -f_B (the sharp edge).
//   Slang's auto-generated backward pass for min() uses a subgradient that
//   may be incorrect for our piecewise convention. The hand-coded backward
//   pass implements the exact piecewise rule used in the C++ reference.
//
// GENERIC DESIGN:
//   CrispSubtract<A, B> is monomorphised at compile time for each (A, B) pair.
//   eval_bwd calls gradient(shape_a, p) and gradient(shape_b, p) — these
//   themselves use bwd_diff, which is valid to call from within a
//   [BackwardDerivativeOf] function (it's a compile-time transformation).

module CrispSubtract;

__exported import IImplicitFunction;  // re-exports gradient<T>, BoundingBox

public struct CrispSubtract<A : IImplicitFunction, B : IImplicitFunction> : IImplicitFunction {
    public A shape_a;
    public B shape_b;

    // f(p) = min(f_A(p), −f_B(p))
    [Differentiable]
    [NoDiffThis]
    public float eval(float3 p) {
        float f_a = shape_a.eval(p);
        float f_b = shape_b.eval(p);
        return min(f_a, -f_b);
    }

    // Gradient of min(f_A, −f_B).
    // Implements the same piecewise rule as crisp_subtract.hpp::eval_gradient.
    [BackwardDerivativeOf(eval)]
    [NoDiffThis]
    public void eval_bwd(inout DifferentialPair<float3> dp, float dResult) {
        float f_a = shape_a.eval(dp.p);
        float f_b = shape_b.eval(dp.p);
        float3 grad;
        if (f_a < -f_b) {
            // A is the active constraint.
            grad = dResult * gradient(shape_a, dp.p);
        } else {
            // Complement of B is active.  ∇(−f_B) = −∇f_B.
            grad = -dResult * gradient(shape_b, dp.p);
        }
        dp = DifferentialPair<float3>(dp.p, dp.d + grad);
    }

    public BoundingBox getBoundingBox() {
        // A − B ⊆ A; bounding the result by A's bounding box is conservative.
        return shape_a.getBoundingBox();
    }
}

// Convenience constructor — type params A, B are inferred from arguments.
public CrispSubtract<A, B> makeCrispSubtract<A : IImplicitFunction, B : IImplicitFunction>(A a, B b) {
    CrispSubtract<A, B> cs;
    cs.shape_a = a;
    cs.shape_b = b;
    return cs;
}
