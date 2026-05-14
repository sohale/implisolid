// MarchingCubes.slang — Marching Cubes polygoniser (Slang/GPU)
//
// Slang realisation of js_iteration_2/marching_cubes.hpp
//
// Design notes (to fill in during implementation):
//
//   Input:
//     - An implicit function T satisfying IImplicitFunction
//     - Grid resolution (dim_x, dim_y, dim_z)
//     - BoundingBox from the implicit function
//
//   Output:
//     - Vertex buffer (float3[])
//     - Face/index buffer (uint3[])
//   Both written to GPU buffers via RWStructuredBuffer or RWByteAddressBuffer.
//
//   GPU parallelism strategy:
//     - One thread group per grid cell (or per grid plane)
//     - Phase 1 (compute shader): evaluate f at all grid corners → scalar field buffer
//     - Phase 2 (compute shader): marching-cubes table lookup per cell → triangle soup
//     - Phase 3 (optional): compact triangle soup → indexed mesh
//
//   The MC lookup tables (edge table, triangle table) will live in a separate
//   MCTables.slang to keep this file algorithm-focused.
//
//   Note: in the C++ OB algorithm, MC is used as step_0 (the seed mesh).
//   This Slang MC can serve both as a standalone polygoniser and as step_0
//   inside a future Slang OB implementation.
//
// Status: placeholder — implementation pending.

module MarchingCubes;

import ../../implicit_function/types;
import ../../implicit_function/IImplicitFunction;
