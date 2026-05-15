// This is a fake/mock/mifixture polygoniser. Used for test of brige interfaceing between ThreeJS and `implisolid/slang/native` sub-project.
// polygnoser stub/mock/fixture for testing (progressive) updates in polygonizer-generation / mesh generation (progressive/reactive) loop
// This generates a surface mesh for a 3D sphere, hard-coded: for test purposes. For testing the bridge


// Inputs:
// Receives R, C : radius and center (x,y,z), and produces vertices, and a way to connect them. in a simple way. Use M,N integer args as number of latitude, logitude steps for vertices. (a grid on surface of wpshere, a spherical-polar is fine, it's jsut to test I/O: I want to test bridge that can transfer vertices and faced faces (array of triple reals, and another aray of triple integers) )
// Outputs:
// data (blocks? for now, the whole block) (verts,faces): vertices and faced faces of surface triangulated mesh, ie,  an array of triple reals, and another aray of triple integers (indives of vertices forming each facet, ie each triangle)
