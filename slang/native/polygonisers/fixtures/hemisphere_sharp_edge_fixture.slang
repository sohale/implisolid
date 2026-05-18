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
    // unverified: cap latitude rows  (theta=[0,π/2] → M+1 vertex rows)
    // unverified: longitude columns  (phi=[0,2π]   → N+1 vertex cols)
    public uint  M;      // cap latitude divisions (M+1 rows, theta=[0,π/2])
    public uint  N;      // longitude divisions    (N+1 cols, phi=[0,2π])
    public uint  _pad0;
    public uint  _pad1;
};

ConstantBuffer<HemiParams> gParams;     // binding(0)
RWStructuredBuffer<float>  gVertices;   // binding(1)  vertex positions
RWStructuredBuffer<uint>   gFaces;      // binding(2)  triangle indices
RWStructuredBuffer<float>  gNormals;    // binding(3)  vertex normals
RWStructuredBuffer<uint>   gSharpEdges; // binding(4)  sharp directed-edge indices
RWStructuredBuffer<uint>   gSharpCorners; // binding(5)  sharp corner indices

// The output counts are inferred in this fixture, but I will add a function for them too.
// Mostly about counts. Cannot be inferred in generic way.

// Sharp-feature counts as a nested struct.
// Inner struct must be a separate named type (Slang does not allow inline anonymous structs).
public struct SharpCounts {
    public uint numSharpEdges;
    public uint numSharpCorners;
};

// Output metadata produced by getOutputMetadata.
// All fields are uint so the WGSL layout is a flat u32 array (no padding):
//   [0]=numVerts  [1]=numFaces  [2]=numNormals  [3]=shape.numSharpEdges  [4]=shape.numSharpCorners
public struct ShapeSharpFeatures {
    public uint numVerts;
    public uint numFaces;
    // Should be the same as numVerts
    public uint numNormals;   // temporarily. Currently, still the same as numVerts for now; will diverge with sharp-vertex duplication
    public SharpCounts shape;
};

public struct DispatchSizes {
    public uint dispVX;
    public uint dispVY;
    public uint dispFX;
    public uint dispFY;
}
/*
RWStructuredBuffer<DispatchSizes> gDispatchSizes; // binding(7)? // not implemented yet,
*/
// Dispatch sizes should be computed here, NOT on JS side, because they are part of the algorithm's logic, not just a detail of how to call it. The bridge will remain generic.
// They must be calculated here bsased on the complex polygonisatoion algorithms, on the Slang side itself
//  and returned to JS via gOutputMetadata.
// This is important for correctness and maintainability, as well as for supporting different shapes and algorithms in the future.
// Dispatch sizes were computed on the JS side (see gpu_compute_bridge.js).
// But need to be calculated here. The logic belogs to here, on the Slang side itself.
// binding(6): one-element structured buffer; written by getOutputMetadata, read back to JS.

RWStructuredBuffer<ShapeSharpFeatures> gOutputMetadata; // binding(6)



[shader("compute")]
[numthreads(64, 1, 1)]
void getOutputMetadata(uint3 tid : SV_DispatchThreadID) {
    uint t = tid.x;
    if (t >= 1) return; // no parallelisation!

    // These can be iinferred in this fixture, but in general, are non-trivial numbers.
    // Currently the code is hard-coded in IMPLISOLID/slang/native/examples/simplical-visualiser/app/renderer/slang_bridge
    /*
    // --- Dispatch sizes ---
    const dispVX = Math.ceil((M + 1) / 8)
    const dispVY = Math.ceil((N + 1) / 8)
    const dispFX = Math.ceil(M / 8)
    const dispFY = Math.ceil(N / 8)

    ...
    const A_count   = (M + 1) * (N + 1)
    const C_idx     = A_count + (N + 1)
    const totalVerts = C_idx + 1
    const capFaces  = M * N * 2
    const totalFaces = capFaces + N

    const vertBytes = totalVerts * 3 * 4
    const faceBytes = totalFaces * 3 * 4
    const normBytes = totalVerts * 3 * 4

    const nVerts     = (M + 1) * (N + 1)
    const nFaces     = M * N * 2
    const vertBytes  = nVerts * 3 * 4     // float32
    const faceBytes  = nFaces * 3 * 4     // uint32
    const normBytes  = nVerts * 3 * 4     // float32, same count as vertices

    from comments:
    // vertices: dispatch ceil((M+1)/8) × ceil((N+1)/8) × 1.
    // faces: dispatch ceil(M/8) × ceil(N/8) × 1.
    //      + dispatch ceil(N/64) × 1 × 1.
    Normals are not dispatched here, in this fixture, yet.
    In strategic solution, normals will be generated in the same polygniisation process/pipeline/chain.


    This is for callers, not for the algorithms.
    In the fixture, simlar things are calculated based on the input paramets (`gParams`, which is part of initial input).
    Since the algorithm is progressive, these will be avilable only after each rout of improvements.
    The inprovments may or may not have observable results.
    Thesea re produced only each time we have obrvable retuslts to give back.
    */
    const uint M = gParams.M;
    const uint N = gParams.N;
    // totalVerts = (M+1)*(N+1) [cap A] + (N+1) [disk rim B] + 1 [disk centre C]
    uint totalVerts = (M + 1u) * (N + 1u) + (N + 1u) + 1u;
    uint capFaces   = M * N * 2u;
    gOutputMetadata[0].numVerts              = totalVerts;
    gOutputMetadata[0].numFaces              = capFaces + N;
    gOutputMetadata[0].numNormals            = totalVerts;   // will differ when sharp-vertex duplication is added
    gOutputMetadata[0].shape.numSharpEdges   = N;            // one per equator segment
    gOutputMetadata[0].shape.numSharpCorners = 3u;           // placeholder: vertices 0 and 1
}


static const float PI = 3.14159265358979323846f;


// geometric shape-specific:  (expired/deprecated comments?)
// Shape: upper hemisphere (z ≥ 0) of radius R centred at (cx,cy,cz).
// Z-up coordinate convention.
// Vertex layout in the returned buffers:
//   Section A  [0 .. (M+1)*(N+1))    curved cap — sphere-gradient normals
//   Section B  [A .. A+N]            disk rim copies — flat normal (0,0,-1)
//   Section C  [A+N+1]               disk centre    — flat normal (0,0,-1)
//   totalVerts = (M+2)*(N+1)+1
//
// Face layout:
//   Cap triangles  [0 .. M*N*2)      two per grid cell
//   Disk fan       [M*N*2 .. M*N*2+N) fan from centre to rim
//   totalFaces = M*N*2 + N
 // geometric shape-specific: (expired/deprecated comments?) ^.

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
// note: either genCapFaces or genDiskFaces is redundant. combine & absorb
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
// note: either genCapFaces or genDiskFaces is redundant. combine & absorb
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

    // Here, (by construction), we know which one is the sharp edge:
    // Directed-edge encoding: e = faceIndex * 3 + localEdgeIndex:
    //
    // It is between (e1,e2) = (B_base + j , B_base + j + 1u)
    // (uint, uint) : sharp_edge =  (B_base + j , B_base + j + 1u)
    // sharp_edge_idx = ( B_base + j, 1); // 1: always the second vertex in this, connecting 2->3
    /*
    uint sharp_edge_idx = (B_base + j)*3 +  1u;
    gSharpEdges[j] = sharp_edge_idx; // mark this edge as the one that participates in a sharp edge
    */
    //
    //
    //   so that JS can decode: faceIdx = e / 3,  v0 = faces[faceIdx*3 + e%3],  v1 = faces[faceIdx*3 + (e%3+1)%3]
    // Disk fan face j has index (cap_tri_count + j) — cap_tri_count already declared above.
    // Its local edge 1 connects gFaces[fidx0+1]=B_base+j to gFaces[fidx0+2]=B_base+j+1.
    uint fi = cap_tri_count + j;
    // Mark this edge as the one that participates in a sharp edge:
    gSharpEdges[j] = fi * 3u + 1u;   // directed equator edge j

    // and now for sharp corners:
    // Placeholder sharp corners: first two equator vertices (vertex indices B_base+0, B_base+1).
    if (j == 0) {
        // lol, N>3
        // for now, some arbitrary vertice is considered as sharp ones in this fixture
        gSharpCorners[0] = B_base + 0;
        gSharpCorners[1] = B_base + 1;
        gSharpCorners[2] = C_idx;
    }
}



[shader("compute")]
[numthreads(64, 1, 1)]
void genCapDiscSharpEdges(uint3 tid : SV_DispatchThreadID) {
    uint j = tid.x;   // cap as disc segment [0 .. N−1]
    if (j >= gParams.N) return;

    // gSharpEdges[j] = index of the edge that participates in a sharp edge: ridge, etc
    // The value of `gSharpEdges[j]` is actually the index of a `gFaces` entry.
    // It encodes via a single integer:
    // The triangle index of the triangle that straddles the sharp edge between cap and disk at equator segment j.
    // if e = gSharpEdges[j],
    // the edge vertices are as follows:
    //   (fi, ei) = (e / 3, e % 3)
    //   v0 = gFaces[fi*3 + ei]
    //   v1 = gFaces[fi*3 + (ei+1) % 3]
    // It reuses the directionality implicitly encoded in the face (winding-number).
    // Let's call it "directed edge" represention. Ironically, it is more compact, shorte, but represent more infomtion: the edge + the edge dir.
    // directed edge ( but not directed in sense iof graph-direciton, but winding-direction ; (that is defined for) of the edge)
    // Brilliant, isn't it? -- Sun 17 May 2026 11:37:45 BST

    // gSharpEdges[j] = ...; // mark edge between cap and disk at equator segment j
    // The values of gSharpEdges[] are already generated in genDiskFaces, but the calculations may move here. I just did that beucae it was  hard to reuse the relevant part of that code.

}
