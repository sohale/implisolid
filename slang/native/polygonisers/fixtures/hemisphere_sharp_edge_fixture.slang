// Hemisphere surface mesh — fixture for testing sharp-edge rendering and normals.
//
// Shape: upper hemisphere (z ≥ 0) of a sphere of radius R centred at origin.
// F-REP: CrispSubtract(UnitSphere, ZNegativeHalfSpace)
//   = min(R² − dot(p,p),  p.z)  — positive inside the hemisphere
//
// Mesh:  parametric, hard-coded in two sections:
//   A) Curved cap  — UV-sphere patch, theta ∈ [0, π/2], Z-up
//   B) Flat disk   — fan triangles from centre to equator circle at z = 0
//
// SHARP EDGE at z = 0 equator (the circle x²+y²=R², z=0):
//   Vertices on the equator are DUPLICATED:
//     - Copy in section A (cap side)  → sphere-gradient normal (radially outward in XY)
//     - Copy in section B (disk side) → flat-disk normal (0, 0, −1) = downward
//   This discontinuity in normal direction IS the sharp edge:
//   lighting changes abruptly as the eye crosses the equator boundary.
//
// Z-up coordinate convention (different from sphere fixture which uses Y-up):
//   x = R · sin(θ) · cos(φ)
//   y = R · sin(θ) · sin(φ)
//   z = R · cos(θ)             θ ∈ [0, π/2],  φ ∈ [0, 2π]
//   At θ=0: (0, 0, R)  — north pole
//   At θ=π/2: (R·cos(φ), R·sin(φ), 0)  — equator at z = 0
//
// Vertex layout (linear index space):
//   [0 .. A_count−1]          Section A: cap vertices     (M+1)·(N+1) entries
//     A[i][j] = i·(N+1)+j     theta row i ∈ [0,M], phi col j ∈ [0,N]
//   [A_count .. A_count+N]     Section B: disk rim copies  N+1 entries
//     B[j]  = A_count + j      same positions as A[M][j], normal (0,0,−1)
//   [A_count+N+1]              Section C: disk centre      1 entry
//     C_idx = A_count + N + 1  position (cx,cy,cz), normal (0,0,−1)
//   Total vertices: (M+2)·(N+1) + 1
//
// Face layout (triangle index space, each triangle = 3 uints):
//   [0 .. M·N·2−1]             Cap triangles   — same quad-split as sphere fixture
//   [M·N·2 .. M·N·2+N−1]       Disk fan        — N triangles, C→Bj→Bj+1
//   Total triangles: M·N·2 + N  = N·(2M+1)
//
// Normals (per-vertex, outward):
//   Section A  — sphere gradient: normalize(−∇f_sphere(p)) = normalize(p − center)
//                Called via gradient<UnitSphere>() from IImplicitFunction.slang.
//   Sections B+C — halfspace gradient: normalize(−∇f_halfspace(p)) = (0,0,−1)
//                Called via gradient<HalfSpace>() — constant, avoids hard-coding.
//
// Bindings (verified from WGSL output — must match gpu_compute_bridge.js):
//   binding(0)  gParams    uniform
//   binding(1)  gVertices  storage read_write
//   binding(2)  gFaces     storage read_write
//   binding(3)  gNormals   storage read_write
//
// Entry dispatch order (all independent except genDiskGeom can follow genCapGeom):
//   genCapGeom   ceil((M+1)/8) × ceil((N+1)/8) × 1  — workgroup (8,8,1)
//   genDiskGeom  ceil((N+2)/64) × 1 × 1              — workgroup (64,1,1)
//   genCapFaces  ceil(M/8)     × ceil(N/8)     × 1  — workgroup (8,8,1)
//   genDiskFaces ceil(N/64)    × 1             × 1  — workgroup (64,1,1)
//
// Build:
//   bash slang/native/scripts/build-sphere-wgsl.sh
//   (script also compiles this fixture and emits hemisphere_mesh_wgsl.js)

module HemiSphereMeshFixture;

import UnitSphere;  // makeSphere(), gradient<UnitSphere>()
import HalfSpace;   // makeZPositiveHalfSpace(), gradient<HalfSpace>()

// Same std140 layout as SphereParams — JS-side uniform packing is identical.
// M here = latitude rows for the CAP only (theta ∈ [0, π/2]).
public struct HemiParams {
    public float radius;
    public float cx;
    public float cy;
    public float cz;
    public uint  M;      // cap latitude divisions (M+1 rows, theta=[0,π/2])
    public uint  N;      // longitude divisions    (N+1 cols, phi=[0,2π])
    public uint  _pad0;
    public uint  _pad1;
};

ConstantBuffer<HemiParams> gParams;     // binding(0)
RWStructuredBuffer<float>  gVertices;   // binding(1)  vertex positions
RWStructuredBuffer<uint>   gFaces;      // binding(2)  triangle indices
RWStructuredBuffer<float>  gNormals;    // binding(3)  vertex normals

static const float PI = 3.14159265358979323846f;

// ---------------------------------------------------------------------------
// Entry 1: genCapGeom
// Writes section A: cap vertex positions and normals.
// One thread per (i, j); dispatch ceil((M+1)/8) × ceil((N+1)/8) × 1.
// A[i][j] index = i*(N+1)+j.
// Normal = outward sphere normal = normalize(−gradient(sphere, p)).
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 1)]
void genCapGeom(uint3 tid : SV_DispatchThreadID) {
    uint i = tid.x;   // theta index [0 .. M]
    uint j = tid.y;   // phi   index [0 .. N]
    if (i > gParams.M || j > gParams.N) return;

    float theta = float(i) * (PI / 2.0f) / float(gParams.M);  // 0 .. π/2
    float phi   = float(j) * 2.0f * PI  / float(gParams.N);   // 0 .. 2π

    float sinT = sin(theta);  float cosT = cos(theta);
    float sinP = sin(phi);    float cosP = cos(phi);

    float3 center = float3(gParams.cx, gParams.cy, gParams.cz);

    // Z-up: x=R·sinT·cosP, y=R·sinT·sinP, z=R·cosT
    float3 p = float3(
        gParams.radius * sinT * cosP + center.x,
        gParams.radius * sinT * sinP + center.y,
        gParams.radius * cosT        + center.z
    );

    uint vidx = i * (gParams.N + 1u) + j;
    gVertices[vidx * 3u + 0u] = p.x;
    gVertices[vidx * 3u + 1u] = p.y;
    gVertices[vidx * 3u + 2u] = p.z;

    // Sphere outward normal via autodiff gradient (same pipeline as future OB/MC).
    // ∇f_sphere = −2·(p−center)  (inward); negate for outward.
    UnitSphere sphere = makeSphere(gParams.radius, center);
    float3 n = normalize(-gradient(sphere, p));
    gNormals[vidx * 3u + 0u] = n.x;
    gNormals[vidx * 3u + 1u] = n.y;
    gNormals[vidx * 3u + 2u] = n.z;
}

// ---------------------------------------------------------------------------
// Entry 2: genDiskGeom
// Writes sections B and C: disk rim copies and disk centre.
// One thread per rim vertex (j ∈ [0,N]) plus the centre (thread N+1).
// Dispatch: ceil((N+2)/64) × 1 × 1.
//
// Positions: B[j] = same position as cap row M at column j.
//            C    = (cx, cy, cz).
// Normals: HalfSpace gradient at z=0 plane = (0,0,1); negate → (0,0,−1).
//          Calling gradient() exercises the HalfSpace autodiff path.
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(64, 1, 1)]
void genDiskGeom(uint3 tid : SV_DispatchThreadID) {
    uint t = tid.x;
    uint Np1    = gParams.N + 1u;
    uint A_count = (gParams.M + 1u) * Np1;
    uint C_idx  = A_count + Np1;   // centre vertex index

    // Disk normal via HalfSpace gradient — same (0,0,−1) for all disk vertices.
    HalfSpace zpos = makeZPositiveHalfSpace();
    float3 disk_normal = normalize(-gradient(zpos, float3(0.0f)));

    if (t < Np1) {
        // Section B: rim copy j = t.
        // Position = cap row M, column t.
        float phi   = float(t) * 2.0f * PI / float(gParams.N);
        float3 center = float3(gParams.cx, gParams.cy, gParams.cz);
        float3 p = float3(
            gParams.radius * cos(phi) + center.x,
            gParams.radius * sin(phi) + center.y,
            center.z                               // z = cz (equator at z=cz)
        );
        uint vidx = A_count + t;
        gVertices[vidx * 3u + 0u] = p.x;
        gVertices[vidx * 3u + 1u] = p.y;
        gVertices[vidx * 3u + 2u] = p.z;
        gNormals[vidx * 3u + 0u] = disk_normal.x;
        gNormals[vidx * 3u + 1u] = disk_normal.y;
        gNormals[vidx * 3u + 2u] = disk_normal.z;
    } else if (t == Np1) {
        // Section C: disk centre.
        gVertices[C_idx * 3u + 0u] = gParams.cx;
        gVertices[C_idx * 3u + 1u] = gParams.cy;
        gVertices[C_idx * 3u + 2u] = gParams.cz;
        gNormals[C_idx * 3u + 0u] = disk_normal.x;
        gNormals[C_idx * 3u + 1u] = disk_normal.y;
        gNormals[C_idx * 3u + 2u] = disk_normal.z;
    }
}

// ---------------------------------------------------------------------------
// Entry 3: genCapFaces
// Writes cap triangles (same quad-split as sphere fixture).
// One thread per grid cell (i, j); dispatch ceil(M/8) × ceil(N/8) × 1.
// Cell (i,j): 2 triangles from 4 corners of the cap quad.
// Face buffer range: [0 .. M·N·6).
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(8, 8, 1)]
void genCapFaces(uint3 tid : SV_DispatchThreadID) {
    uint i = tid.x;   // [0 .. M−1]
    uint j = tid.y;   // [0 .. N−1]
    if (i >= gParams.M || j >= gParams.N) return;

    uint Np1 = gParams.N + 1u;
    uint v00 =  i      * Np1 +  j;
    uint v10 = (i + 1u)* Np1 +  j;
    uint v01 =  i      * Np1 + (j + 1u);
    uint v11 = (i + 1u)* Np1 + (j + 1u);

    uint fidx0 = (i * gParams.N + j) * 6u;

    gFaces[fidx0 + 0u] = v00;  gFaces[fidx0 + 1u] = v10;  gFaces[fidx0 + 2u] = v01;
    gFaces[fidx0 + 3u] = v10;  gFaces[fidx0 + 4u] = v11;  gFaces[fidx0 + 5u] = v01;
}

// ---------------------------------------------------------------------------
// Entry 4: genDiskFaces
// Writes disk fan triangles.
// One thread per fan segment j ∈ [0, N−1]; dispatch ceil(N/64) × 1 × 1.
// Triangle j: (C_idx, B[j], B[j+1])  — CCW seen from below (normal = −Z).
//   C_idx = (M+1)·(N+1) + (N+1)
//   B[j]  = (M+1)·(N+1) + j
// Face buffer range: [M·N·6 .. (M·N·2+N)·3).
// ---------------------------------------------------------------------------
[shader("compute")]
[numthreads(64, 1, 1)]
void genDiskFaces(uint3 tid : SV_DispatchThreadID) {
    uint j = tid.x;   // fan segment [0 .. N−1]
    if (j >= gParams.N) return;

    uint Np1    = gParams.N + 1u;
    uint A_count = (gParams.M + 1u) * Np1;
    uint B_base = A_count;          // start of section B
    uint C_idx  = A_count + Np1;   // section C (centre)

    uint cap_tri_count = gParams.M * gParams.N * 2u;
    uint fidx0 = (cap_tri_count + j) * 3u;

    // Fan triangle: centre → B[j] → B[j+1]
    // B[j+1] when j = N−1 gives B[N], which is the seam-duplicate of B[0].
    gFaces[fidx0 + 0u] = C_idx;
    gFaces[fidx0 + 1u] = B_base + j;
    gFaces[fidx0 + 2u] = B_base + j + 1u;
}
