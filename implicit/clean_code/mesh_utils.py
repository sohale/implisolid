import numpy as np

VERBOSE = False


def centroids(verts, faces):
    # print("faces: ", faces.shape)
    # print(verts[faces[:],:].shape)
    c = np.mean(verts[faces[:], :], axis=1)
    # print(c.shape)
    return c

# The only useful method in our method. This needs to be done once and can be reused.
def make_neighbour_faces_of_vertex(faces):
    """ neighbour_faces_of_vertex is a list. index=vertex, v1,v2,v3 """
    vertex_count = np.max(np.max(faces)) + 1
    neighbour_faces_of_vertex = {}  # np.zeros( (vertex_count,3) , dtype=np.type) - 1
    for fi in range(faces.shape[0]):
        for vi in range(3):
            v1 = faces[fi, vi]
            if v1 not in neighbour_faces_of_vertex:
                neighbour_faces_of_vertex[v1] = [fi]
            else:
                assert fi not in neighbour_faces_of_vertex[v1]
                neighbour_faces_of_vertex[v1].append(fi)
    # print(neighbour_faces_of_vertex)
    # print(map( lambda k: len(neighbour_faces_of_vertex[k]), neighbour_faces_of_vertex   )) #for k,v in neighbour_faces_of_vertex:
    #    print(len(v), end="")
    return neighbour_faces_of_vertex

#import vectorized

if __name__ == '__main__':
    pass
