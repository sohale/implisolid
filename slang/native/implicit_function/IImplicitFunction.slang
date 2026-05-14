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
//   4. __exported import types — re-exports BoundingBox to all importers of this
//      module so they don't need a separate `import types`.
//
//   5. Slang module visibility: all exported symbols must be `public`.
//      Struct fields default to `internal` — see types.slang.
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

public interface IImplicitFunction {
    // Implicit function value at a single point p.
    // f > 0: inside, f = 0: surface, f < 0: outside  (F-REP convention).
    // [Differentiable] must be on the interface declaration (not only the struct)
    // so that bwd_diff() works when T is a generic parameter.
    [Differentiable]
    float eval(float3 p);

    // Bounding box in world space. Used by the polygoniser to set grid extent.
    BoundingBox getBoundingBox();
};

// ---------------------------------------------------------------------------
// _evalForDiff<T> — free wrapper so bwd_diff can differentiate eval().
//
// bwd_diff() cannot take a reference to a non-static member function;
// `bwd_diff(shape.eval)` is not valid Slang. This wrapper with an explicit
// `no_diff T shape` parameter is the correct idiom: shape passes through
// unchanged; only `p` participates in the AD graph.
// ---------------------------------------------------------------------------
[Differentiable]
public float _evalForDiff<T : IImplicitFunction>(no_diff T shape, float3 p) {
    return shape.eval(p);
}

// ---------------------------------------------------------------------------
// gradient<T> — derive ∇f via autodiff. Works for any IImplicitFunction.
//
// Implementing types may supply [BackwardDerivativeOf(eval)] to replace the
// auto-generated backward pass with a hand-coded one (same result, faster).
// bwd_diff(_evalForDiff<T>) dispatches through shape.eval(p) and will use
// T's hand-coded backward derivative when one is provided.
//
// Returns ∇f(p). Not unit-length. Normalize for rendering normals.
// Direction: points toward the solid interior (inward normal), matching C++ ref.
// ---------------------------------------------------------------------------
public float3 gradient<T : IImplicitFunction>(T shape, float3 p) {
    var dp = diffPair(p, float3(0.0f, 0.0f, 0.0f));
    bwd_diff(_evalForDiff<T>)(shape, dp, 1.0f);
    return dp.d;
}
