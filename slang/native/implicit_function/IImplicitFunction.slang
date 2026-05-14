// IImplicitFunction — core interface for all implicit solids.
//
// Slang realisation of mp5_implicit::implicit_function.
// Reference: js_iteration_2/implicit_function/implicit_function.hpp
//
// KEY DESIGN DECISIONS vs the C++ reference:
//
//   1. Single-point interface (not batched).
//      C++ eval_implicit / eval_gradient take vectorized_vect (a batch of points)
//      for cache efficiency. Here each method handles one point; GPU thread groups
//      provide the equivalent parallelism. One thread = one point.
//
//   2. gradient() is NOT in the interface. It is derived via Slang's autodiff from
//      eval(). Implementing types can provide a hand-coded backward pass via
//      [BackwardDerivativeOf(eval)] for performance, but consumers always call the
//      free function gradient<T>(shape, p) — never a method.
//      This is the standard Slang/ML idiom: require the function, derive its derivative.
//
//   3. eval() is marked [Differentiable]. This attribute must appear on the interface
//      declaration (not only on the implementing struct) so that bwd_diff() works
//      through a generic T : IImplicitFunction parameter.
//
//   4. __exported import types — BoundingBox must be visible to any file that
//      imports IImplicitFunction; without __exported, callers would need to
//      separately import types themselves.
//
// USING THE INTERFACE IN A GENERIC POLYGONISER:
//
//   // Correct: generic over T, monomorphized per type — no virtual dispatch on GPU.
//   void doWork<T : IImplicitFunction>(T shape, float3 p) {
//       float  f    = shape.eval(p);
//       float3 grad = gradient(shape, p);   // free function below
//   }
//
// Precision: currently fixed at float32. Can be generalised to
//   interface IImplicitFunction<T : __BuiltinFloatingPointType>
//   when multiple precisions are needed.

module IImplicitFunction;

__exported import types;   // re-exports BoundingBox to all importers of this module

interface IImplicitFunction {
    // Implicit function value at a single point p.
    // f > 0: inside, f = 0: surface, f < 0: outside  (F-REP convention).
    // Must be [Differentiable] so bwd_diff() works through generic parameters.
    [Differentiable]
    float eval(float3 p);

    // Bounding box for the solid in world space.
    // Used by the polygoniser to set grid extent.
    BoundingBox getBoundingBox();
};

// ---------------------------------------------------------------------------
// _evalForDiff<T> — free wrapper so bwd_diff can differentiate eval().
//
// bwd_diff() cannot take a reference to a non-static member function —
// `bwd_diff(shape.eval)` is not valid Slang. Wrapping eval() in a free
// function with an explicit `no_diff T shape` parameter is the correct idiom:
// shape is passed through unchanged; only `p` is in the AD graph.
// ---------------------------------------------------------------------------
[Differentiable]
float _evalForDiff<T : IImplicitFunction>(no_diff T shape, float3 p) {
    return shape.eval(p);
}

// ---------------------------------------------------------------------------
// gradient<T> — derive ∇f via autodiff. Works for any IImplicitFunction.
//
// Implementing types may supply [BackwardDerivativeOf(eval)] to replace the
// auto-generated backward pass with a hand-coded one (same result, faster).
// bwd_diff(_evalForDiff<T>) dispatches through shape.eval(p), which will use
// that hand-coded backward pass when T provides one.
//
// Returns ∇f(p). Not unit-length; magnitude varies. Normalize for rendering.
// Direction convention: points toward the solid interior (inward normal).
// ---------------------------------------------------------------------------
float3 gradient<T : IImplicitFunction>(T shape, float3 p) {
    var dp = diffPair(p, float3(0.0f, 0.0f, 0.0f));
    bwd_diff(_evalForDiff<T>)(shape, dp, 1.0f);
    return dp.d;
}
