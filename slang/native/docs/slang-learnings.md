# Slang Learnings — ImpliSolid Native Slang Port

Accumulated from first compilation session (May 2026).
Slangc version: **2025.24.3-1** (`slangc -v`)
Slangc binary: `/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc`

See also: `../scripts/sanity-test-implicit.sh` — runnable sanity check that exercises everything here.

---

## 1. Module system

### Declaring a module
Put `module name;` as the first declaration in the file. Module name should match the filename (without `.slang`).

```slang
// types.slang
module types;
public struct BoundingBox { ... };
```

### Importing
```slang
import types;            // makes types symbols available in THIS file only
__exported import types; // re-exports types symbols to anyone who imports THIS module
```

**Critical rule:** If module A imports module B, consumers of A do NOT automatically see B's symbols unless A uses `__exported import B`. Unlike C++ `#include` — closer to Rust's `pub use`.

### Compile with `-I` for module search path
```bash
slangc -I /path/to/module/dir  entry_shader.slang  -target spirv  -o out.spv
```
The `-I` path is where slangc resolves `import name` → `name.slang`.

---

## 2. Access control — the biggest gotcha

**Slang module members default to `internal` (not `public`).** Unlike HLSL/GLSL where everything in a struct is accessible. Every symbol used outside the module needs explicit `public`:

```slang
// WRONG — fields invisible to importers:
public struct BoundingBox {
    float3 minCorner;   // internal — caller sees the type but cannot read fields
    float3 maxCorner;
};

// CORRECT:
public struct BoundingBox {
    public float3 minCorner;
    public float3 maxCorner;
};
```

Applies to: struct fields, interface declarations, free functions, everything.
Error when wrong: `error 30600: 'fieldName' is not accessible from the current context.`

---

## 3. Autodiff — key patterns

### `[Differentiable]` must be on the INTERFACE declaration

```slang
public interface IImplicitFunction {
    [Differentiable]     // ← must be here, on the interface method
    float eval(float3 p);
};
```

If only on the struct method (not the interface), `bwd_diff()` fails when called through a generic `T : IImplicitFunction` parameter.

### `bwd_diff(shape.eval)` is INVALID — use a free wrapper

```slang
// WRONG — error 30098: non-static function reference not allowed
bwd_diff(shape.eval)(dp, 1.0f);

// CORRECT — wrap eval in a free function:
[Differentiable]
public float _evalForDiff<T : IImplicitFunction>(no_diff T shape, float3 p) {
    return shape.eval(p);
}

public float3 gradient<T : IImplicitFunction>(T shape, float3 p) {
    var dp = diffPair(p, float3(0.0f, 0.0f, 0.0f));
    bwd_diff(_evalForDiff<T>)(shape, dp, 1.0f);
    return dp.d;
}
```

`no_diff T shape` tells the AD system: shape is a constant; only `p` has derivatives.
`_evalForDiff<T>` is monomorphized per concrete type — no virtual dispatch.

### Hand-coded backward pass: `[BackwardDerivativeOf(eval)]`

Registers a hand-coded backward pass so Slang uses it instead of auto-generating one:

```slang
[BackwardDerivativeOf(eval)]
[NoDiffThis]
public void eval_bwd(inout DifferentialPair<float3> dp, float dResult) {
    float3 grad_p = dResult * (-2.0f * (dp.p - center));  // for unit sphere: ∂f/∂p
    dp = DifferentialPair<float3>(dp.p, dp.d + grad_p);   // accumulate, don't assign
}
```

- `dp.p` = primal value (original `p`)
- `dp.d` = accumulated differential (gradient accumulator — **add to it**, do not assign)
- `dResult` = upstream gradient (1.0 when this is the root of the AD call)

Verify wiring in WGSL output: `s_bwd_prop_UnitSphere_eval_0` should call `UnitSphere_eval_bwd_0`.

### `[NoDiffThis]` — for structs that are not `IDifferentiable`

When a struct has `[Differentiable]` methods but is not itself `IDifferentiable`, Slang warns for every access to `this.field`:

```
warning 31159: There is no derivative calculated for member 'center' because
the parent struct is not differentiable. Consider using [NoDiffThis]...
```

Fix: add `[NoDiffThis]` to both forward and backward methods:

```slang
[Differentiable]
[NoDiffThis]
public float eval(float3 p) { ... }

[BackwardDerivativeOf(eval)]
[NoDiffThis]
public void eval_bwd(inout DifferentialPair<float3> dp, float dResult) { ... }
```

**When to remove `[NoDiffThis]`:** If you ever need gradient-based shape fitting (optimising `radius`/`center` to match data), make the struct implement `IDifferentiable` and remove `[NoDiffThis]`. The polygoniser use case never needs this — it only differentiates f w.r.t. the query point `p`.

---

## 4. Generics — mandatory for GPU, no virtual dispatch

Always write polygonisers as generics:

```slang
// CORRECT — monomorphized per concrete type, zero overhead
void polygonise<T : IImplicitFunction>(T shape, float3 gridMin, float3 gridMax) { ... }
```

The generic is specialized at compile time for each concrete type. The GPU never sees a vtable.

---

## 5. Convenience constructors

Slang structs do not support C++-style constructor overloads. Use free functions:

```slang
public UnitSphere makeUnitSphere() { ... }
public UnitSphere makeSphere(float radius) { ... }
public UnitSphere makeSphere(float radius, float3 center) { ... }
```

---

## 6. Compilation targets and flags

```bash
SLANGC=/dataneura/gpu-experimentations/experiments/20_slang_shaders/slang/build/RelWithDebInfo/bin/slangc

# SPIRV — most universal sanity check:
$SLANGC -I <module_dir> shader.slang -target spirv -o out.spv

# WGSL — eventual browser/WebGPU target:
$SLANGC -I <module_dir> shader.slang -target wgsl  -o out.wgsl

# HLSL — needs explicit entry when not inferred:
$SLANGC -I <module_dir> shader.slang -target hlsl -entry main -stage compute -o out.hlsl

# Check version:
$SLANGC -v

# List all available targets:
$SLANGC -h target
```

SPIRV and WGSL both work without `-entry`/`-stage` when `[shader("compute")]` is in source.
HLSL requires `-entry <name> -stage compute` explicitly.

Library modules (no `[shader]` entry point) cannot be compiled to binary targets directly —
always compile through a test shader that has an entry point.

---

## 7. Checklist for each new shape primitive

- [ ] `public struct ShapeName : IImplicitFunction`
- [ ] All fields `public float ...` (no `no_diff` on fields needed — struct is not `IDifferentiable`)
- [ ] `[Differentiable] [NoDiffThis] public float eval(float3 p)`
- [ ] `[BackwardDerivativeOf(eval)] [NoDiffThis] public void eval_bwd(inout DifferentialPair<float3> dp, float dResult)`
- [ ] `public BoundingBox getBoundingBox()`
- [ ] Free functions for construction (`makeShapeName(...)`)
- [ ] `__exported import IImplicitFunction` (not bare `import`) so callers get `gradient<T>`

### CSG ops pattern

CSG ops are generic structs over two shape types:

```slang
// Crisp union: f = max(f_a, f_b)  (R-function convention, f > 0 = inside)
public struct CrispUnion<A : IImplicitFunction, B : IImplicitFunction> : IImplicitFunction {
    public A left;
    public B right;

    [Differentiable][NoDiffThis]
    public float eval(float3 p) {
        return max(left.eval(p), right.eval(p));
    }
    // eval_bwd: gradient of max — passes through to whichever operand is larger
    // BoundingBox: union of both bboxes
}
// Crisp intersection: min(f_a, f_b)
// Crisp subtract:     min(f_a, -f_b)
```

Note: The backward pass of `max(a, b)` is a sub-gradient — gradient flows to whichever operand is larger. Slang's autodiff handles this automatically for `max()`; you can omit `[BackwardDerivativeOf]` and let it auto-generate.

---

## 8. Open questions / to verify after slangc updates

- `[BackwardDerivativeOf]` exact signature when struct is not `IDifferentiable` — worked in 2025.24.3; re-check after upgrades
- `no_diff` on struct fields (vs. function parameters) — compiler accepted it silently; `[NoDiffThis]` is the documented path
- `bwd_diff(_evalForDiff<T>)` with generic `T` — worked but not explicitly in Slang spec; watch for regressions
- WGSL buffer alignment: `@align(16)` on `f32` fields — verify correct for WebGPU buffer layout on the JS side
