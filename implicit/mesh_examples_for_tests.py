import numpy as np


def testcase_square():
    v = np.array([[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0] ])*10 / 4.
    f = np.array([[0, 1, 2], [0, 2, 3]])
    return v, f


def testcase_cube():
    faces_xyz = []
    for d in range(3):
        dims = np.array([(d+1) % 3, (d+2) %3, d], dtype=np.int)  # new dimensions that x, y, z are mapped into
        for s in range(2):
            side_verts = []
            for x in range(2):
                for y in range(2):
                    xyz = np.array([x, y, s], dtype=np.int)[dims]
                    side_verts.append(tuple(xyz))
            #reorder = np.array([0,1,3,2,0])
            #reorder = [[0, 1], [1, 3], [3, 2], [2, 0]]
            reorder = [[0, 1, 3], [3, 2, 0]]
            verts_reorders = map( lambda face: [side_verts[face[0]], side_verts[face[1]], side_verts[face[2]] ] , reorder)
            faces_xyz += verts_reorders
    #print faces_xyz
    vertdict = {}
    vert_index_dict = {}
    pure_verts = []
    faces = []
    for i in range(len(faces_xyz)):
        f = faces_xyz[i]
        face1 = []
        for v in f:
            key = str(v)
            if key in vertdict:
                idx = vert_index_dict[key]
            else:
                pure_verts.append(v)
                idx = len(pure_verts) - 1
                vert_index_dict[key] = idx
            vertdict[key] = v
            face1.append(idx)
            del idx
        faces.append(face1)


    #for v in  vertdict:
    #    print v
    #print vertdict
    #print vert_index_dict
    #print pure_verts
    for i in range(len(faces_xyz)):
        f = faces_xyz[i]
        for v in f:
            xyz = vertdict[str(v)]
    verts = pure_verts
    #print
    #print faces
    v =  np.array(verts)
    f =  np.array(faces, dtype=np.int64)
    return v, f
