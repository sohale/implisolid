// This is a fake/mock/mifixture polygoniser. Used for test of brige interfaceing between ThreeJS and `implisolid/slang/native` sub-project.
// polygnoser stub/mock/fixture for testing (progressive) updates in polygonizer-generation / mesh generation (progressive/reactive) loop
// This generates a surface mesh for a 3D sphere, hard-coded: for test purposes. For testing the bridge

// Sphere surface mesh generator — stub/fixture for bridge integration testing.
//
// Purpose: produce a triangulated UV-sphere mesh via GPU compute so the
// JS↔Slang bridge (gpu_compute_bridge.js) can be tested before any real
// polygoniser (MC, OB) is implemented.
//
// Inputs (uniform buffer):
//   radius       — sphere radius
//   cx, cy, cz   — sphere centre
//   M            — latitude  divisions (vertex rows:    M+1, face rows:    M)
//   N            — longitude divisions (vertex columns: N+1, face columns: N)
//
// Outputs (storage buffers — tightly packed, compatible with Three.js BufferAttribute):
//   gVertices    — float[], length (M+1)*(N+1)*3   layout: [x0,y0,z0, x1,y1,z1, ...]
//   gFaces       — uint[],  length M*N*2*3          layout: [a0,b0,c0, a1,b1,c1, ...]
//                  Two triangles per grid cell; M*N cells; 2*M*N triangles total.
//
// IMPORTANT: Uses RWStructuredBuffer<float> / RWStructuredBuffer<uint> (not float3/uint3)
// to avoid the vec3<f32> 16-byte alignment padding that WGSL storage buffers impose.
// Three.js BufferAttribute expects tightly-packed [x,y,z, x,y,z, ...] with stride 3.
//
// Build (run from implisolid/ or slang/native/):
//   bash slang/native/scripts/build-sphere-wgsl.sh
// This compiles to WGSL and writes slang/native/examples/simplical-visualiser/app/renderer/slang_bridge/generated/sphere_mesh_wgsl.js
//
// Cross-reference: see IImplicitFunction.slang — the UnitSphere struct evaluates f(p)
// and gradient. This fixture skips evaluation entirely and hard-codes the sphere
// surface parametrically, for fast bridge testing only.

module SphereMeshFixture;

static const float PI = 3.14159265358979323846f;

// Sphere parameters — bound as a uniform buffer.
// Declared as a plain struct; Slang maps ConstantBuffer<T> → var<uniform> in WGSL.
public struct SphereParams {
    public float radius;
    public float cx;
    public float cy;
    public float cz;
    public uint  M;   // latitude  step count  (rows of vertices = M+1)
    public uint  N;   // longitude step count  (cols of vertices = N+1)
    public uint  _pad0;  // explicit padding to 32 bytes
    public uint  _pad1;
};

// Binding layout (assigned by Slang, verified from WGSL output):
//   group(0) binding(0)  →  gParams    (uniform)
//   group(0) binding(1)  →  gVertices  (storage read_write)
//   group(0) binding(2)  →  gFaces     (storage read_write)
ConstantBuffer<SphereParams> gParams;
RWStructuredBuffer<float>    gVertices;  // (M+1)*(N+1)*3 floats
RWStructuredBuffer<uint>     gFaces;     // M*N*2*3 uints

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

    float theta = float(i) * PI / float(gParams.M);          // 0 .. π
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
}

// ---------------------------------------------------------------------------
// Entry 2: genFaces
// One thread per grid cell (i, j).
// Dispatch: ceil(M/8) × ceil(N/8) × 1
// Writes 2 triangles per cell (CCW winding, outward normal).
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 1)]
void genFaces(uint3 tid : SV_DispatchThreadID) {
    uint i = tid.x;   // [0 .. M-1]
    uint j = tid.y;   // [0 .. N-1]

    if (i >= gParams.M || j >= gParams.N) return;

    uint Np1 = gParams.N + 1u;

    // Corners of this quad cell in the (i,j) grid
    uint v00 =  i      * Np1 +  j;
    uint v10 = (i + 1u)* Np1 +  j;
    uint v01 =  i      * Np1 + (j + 1u);
    uint v11 = (i + 1u)* Np1 + (j + 1u);

    // Two triangles per cell; cell index → face buffer offset
    uint cellIdx = i * gParams.N + j;
    uint fidx0 = cellIdx * 6u;   // each triangle = 3 uints, 2 per cell = 6 per cell

    // Triangle 0: bottom-left, top-left, bottom-right  (CCW from outside)
    gFaces[fidx0 + 0u] = v00;
    gFaces[fidx0 + 1u] = v10;
    gFaces[fidx0 + 2u] = v01;

    // Triangle 1: top-left, top-right, bottom-right
    gFaces[fidx0 + 3u] = v10;
    gFaces[fidx0 + 4u] = v11;
    gFaces[fidx0 + 5u] = v01;
}
