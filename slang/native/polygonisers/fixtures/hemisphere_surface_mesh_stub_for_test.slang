// This file is deprecated. DONT CHANGE THIS FILE. This will be absorbed into hemisphere_sharp_edge_fixture.slang
// Hemi-Sphere surface mesh generator — stub/fixture for bridge integration testing.
//
// Purpose: produce a pre-triangulated hemisphere mesh + per-vertex normals via GPU
// compute so the JS↔Slang bridge (gpu_compute_bridge.js) can be tested before
// any real polygoniser (MC, OB) is implemented.
//
// Mesh generation (genVertices, genFaces) is hard-coded parametrically.
// Normal generation (genNormals) calls gradient() from UnitSphere.slang —
// the exact autodiff path that the real MC/OB polygoniser will use.
//
// Why normals are per-vertex (not per-face):
//   The GPU interpolates per-vertex attributes (normals, UVs) across each
//   triangle face (Gouraud/Phong). Per-vertex normals = smooth shading.
//   For sharp edges in OB: duplicate vertices at creases — each copy gets a
//   different normal; Three.js then renders both sides without interpolation.
//   This fixture is a smooth sphere so all vertices are unique and shared.
//
// Normal sign: ∇f(p) = −2·(p−center) points INWARD for our F-REP sphere
//   (f > 0 inside). genNormals negates the gradient to get the outward normal.
//
// the flat side of hemisphere is subtracting the part of space with "negative Z". this should stitch the open side
// maybe simply by connecting each edhe to centre ( creating a fan/hedge-hog flat surface)
// The resulting mesh has a flat side, and importandly sharp curve edges at Z=0. (the border of the disk, the curve: (Z=0 && x^2+y^2==R*R) ).
// this curve will be naturally part of the mesh if we cut the Z<0 and add that flat-fan side.
//
// Inputs (uniform buffer):
//   radius       — hemisphere radius
//   cx, cy, cz   — sphere centre
//   M            — latitude  divisions (vertex rows:    M+1, face rows:    M)
//   N            — longitude divisions (vertex columns: N+1, face columns: N)
//  (no parameter for the z-cut: always happens at Z=0 )
//
// Outputs (storage buffers — tightly packed, compatible with Three.js BufferAttribute):
//   gVertices    — float[], length (M+1)*(N+1)*3   layout: [x0,y0,z0, ...]
//   gFaces       — uint[],  length M*N*2*3          layout: [a0,b0,c0, ...]
//   gNormals     — float[], length (M+1)*(N+1)*3   layout: [nx0,ny0,nz0, ...]
//                  outward unit normals; one per vertex; same index space as gVertices.
//
// Binding layout (verified from WGSL output — must match gpu_compute_bridge.js):
//   group(0) binding(0)  →  gParams    (uniform)
//   group(0) binding(1)  →  gVertices  (storage read_write)
//   group(0) binding(2)  →  gFaces     (storage read_write)
//   group(0) binding(3)  →  gNormals   (storage read_write)
//
// Dispatch order (WebGPU guarantees sequential execution within one command encoder):
//   1. genVertices  — fills gVertices
//   2. genFaces     — fills gFaces  (independent of genVertices; can overlap in future)
//   3. genNormals   — reads gVertices, writes gNormals  (must follow genVertices)
//
// Build (run from implisolid/ or slang/native/):
//   bash slang/native/scripts/build-sphere-wgsl.sh
// Writes: slang/native/examples/simplical-visualiser/app/renderer/slang_bridge/generated/sphere_mesh_wgsl.js

module SphereMeshFixture;

import UnitSphere;   // brings in makeSphere(), gradient<T>(), BoundingBox

static const float PI = 3.14159265358979323846f;

// Sphere parameters — bound as a uniform buffer.
public struct SphereParams {
    public float radius;
    public float cx;
    public float cy;
    public float cz;
    public uint  M;       // latitude  step count  (rows of vertices = M+1)
    public uint  N;       // longitude step count  (cols of vertices = N+1)
    public uint  _pad0;   // explicit padding to 32 bytes (std140 alignment)
    public uint  _pad1;
};

ConstantBuffer<SphereParams> gParams;     // binding(0)
RWStructuredBuffer<float>    gVertices;   // binding(1)  (M+1)*(N+1)*3 floats
RWStructuredBuffer<uint>     gFaces;      // binding(2)  M*N*2*3 uints
RWStructuredBuffer<float>    gNormals;    // binding(3)  (M+1)*(N+1)*3 floats

// ---------------------------------------------------------------------------
// Entry 1: genVertices
// One thread per (i, j) vertex.
// Dispatch: ceil((M+1)/8) × ceil((N+1)/8) × 1
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 1)]
void genVertices(uint3 tid : SV_DispatchThreadID) {
    uint i = tid.x;   // latitude  index [0 .. M]
    uint j = tid.y;   // longitude index [0 .. N]

    if (i > gParams.M || j > gParams.N) return;

    // todo: fix: theta will have a a limited range (interval)  for this z-hemisphere:
    float theta = float(i) * PI? / float(gParams.M) - PI/2.0; ?         // 0 .. π/2 ? or -π/4..π/4?
    float phi   = float(j) * 2.0f * PI / float(gParams.N);   // 0 .. 2π

    float sinT = sin(theta);
    float cosT = cos(theta);
    float sinP = sin(phi);
    float cosP = cos(phi);

    float x = gParams.radius * sinT * cosP + gParams.cx;
    float y = gParams.radius * cosT         + gParams.cy;
    float z = gParams.radius * sinT * sinP  + gParams.cz;

    uint vidx = i * (gParams.N + 1u) + j;
    gVertices[vidx * 3u + 0u] = x;
    gVertices[vidx * 3u + 1u] = y;
    gVertices[vidx * 3u + 2u] = z;

    // now add the closing flat lid/cap:
    // bi n interms of vertices: only needs one vvertex: centre!
    int lid_first_vertex_index = ...;
    if (once){
        uint centrepoint_idx = ...;
        gVertices[centrepoint_idx * 3u + 0u] = gParams.cx;
        gVertices[centrepoint_idx * 3u + 1u] = gParams.cy;
        gVertices[centrepoint_idx * 3u + 2u] = gParams.cz;
    }
}

// ---------------------------------------------------------------------------
// Entry 2: genFaces
// One thread per grid cell (i, j).
// Dispatch: ceil(M/8) × ceil(N/8) × 1
// Writes 2 triangles per cell (CCW winding, outward normal).
// Independent of genVertices — can run in any order.
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 1)]
void genFaces(uint3 tid : SV_DispatchThreadID) {
    uint i = tid.x;   // [0 .. M-1]
    uint j = tid.y;   // [0 .. N-1]

    if (i >= gParams.M || j >= gParams.N) return;

    uint Np1 = gParams.N + 1u;

    uint v00 =  i      * Np1 +  j;
    uint v10 = (i + 1u)* Np1 +  j;
    uint v01 =  i      * Np1 + (j + 1u);
    uint v11 = (i + 1u)* Np1 + (j + 1u);

    uint cellIdx = i * gParams.N + j;
    uint fidx0   = cellIdx * 6u;

    // Triangle 0: v00, v10, v01  (CCW from outside)
    gFaces[fidx0 + 0u] = v00;
    gFaces[fidx0 + 1u] = v10;
    gFaces[fidx0 + 2u] = v01;

    // Triangle 1: v10, v11, v01
    gFaces[fidx0 + 3u] = v10;
    gFaces[fidx0 + 4u] = v11;
    gFaces[fidx0 + 5u] = v01;

    // now add the closing flat lid/cap:
    // connect all equator points to center, a fan of triangles
    int lid_first_vertex_index = ...; // same as above
    int lid_first_face_index = ...;
    int v_centre = ...;

    if (....) { // if based on `tid`, we are at a point on sharp-equator,
        lidfanfidx0 = ...
        gFaces[lidfanfidx0 + 0u] = v_centre;
        gFaces[lidfanfidx0 + 1u] = v10;
        gFaces[lidfanfidx0 + 2u] = v01;
    }
}

// ---------------------------------------------------------------------------
// Entry 3: genNormals
// One thread per (i, j) vertex — same dispatch as genVertices.
// Dispatch: ceil((M+1)/8) × ceil((N+1)/8) × 1
// MUST run after genVertices (reads gVertices).
//
// Uses gradient<UnitSphere>() from IImplicitFunction.slang, which calls
// bwd_diff(_evalForDiff<UnitSphere>) — the autodiff path that the real
// MC/OB polygoniser will use for arbitrary implicit functions.
//
// Gradient direction: ∇f(p) = −2·(p−center)  (inward for F-REP sphere).
// Outward normal     = normalize(−∇f(p))      = normalize(p − center).
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 1)]
void genNormals(uint3 tid : SV_DispatchThreadID) {
    uint i = tid.x;
    uint j = tid.y;

    if (i > gParams.M || j > gParams.N) return;

    uint vidx = i * (gParams.N + 1u) + j;

    // Read vertex position written by genVertices.
    float3 p = float3(
        gVertices[vidx * 3u + 0u],
        gVertices[vidx * 3u + 1u],
        gVertices[vidx * 3u + 2u]
    );

    // Reconstruct the sphere implicit function from the uniform params.
    UnitSphere sphere = makeSphere(gParams.radius, float3(gParams.cx, gParams.cy, gParams.cz));
    HalfPlane halfplane_z = makeHalfPlane(float3(0.0f, 0.0f, 1.0f), 0.0f);   // z>=0
    Hemisphere hemisphere = intersection(sphere, halfplane_z);

    // gradient<T>() calls bwd_diff(_evalForDiff<T>) — the autodiff pipeline.
    // For UnitSphere: ∇f(p) = eval_bwd → −2·(p−center).
    float3 grad = gradient(hemisphere, p);

    // F-REP gradient points inward; negate for outward surface normal.
    float3 n = normalize(-grad);

    gNormals[vidx * 3u + 0u] = n.x;
    gNormals[vidx * 3u + 1u] = n.y;
    gNormals[vidx * 3u + 2u] = n.z;
}
