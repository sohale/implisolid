#include <iostream>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <math.h>
#include <cassert>
#include <map>
// #include <vector>
#include <string>
#include <tuple>
#include <fstream>

using namespace std;
using namespace mp5_implicit;

#define ASSERTS 1
//#define VERBOSE  1
typedef float REAL;
typedef struct {
   REAL x, y, z;
} XYZ;


REAL ABS(REAL x){
  if(x<0)
    return -x;
  return x;
}

typedef boost::multi_array<REAL, 2> verts_t;
typedef boost::multi_array<int, 2> faces_t;
typedef std::vector<int> vector_int;
typedef std::vector<std::vector<int>> neighbour;
typedef pair<verts_t, faces_t> vf_t;


REAL norm_2(REAL x, REAL y, REAL z){
  REAL norm = sqrt(x*x + y*y + z*z);
  return norm;
}


REAL compute_average_edge_length(const faces_t& faces, const verts_t& verts){
  int nfaces = faces.shape()[0];
  REAL edge_length;
  for (int j=0; j<nfaces; j++){
    edge_length += norm_2((verts[faces[j][0]][0] - verts[faces[j][1]][0], verts[faces[j][0]][1] - verts[faces[j][1]][1], verts[faces[j][0]][2] - verts[faces[j][1]][2]));
    edge_length += norm_2((verts[faces[j][0]][0] - verts[faces[j][2]][0], verts[faces[j][0]][1] - verts[faces[j][2]][1], verts[faces[j][0]][2] - verts[faces[j][2]][2]));
    edge_length += norm_2((verts[faces[j][2]][0] - verts[faces[j][1]][0], verts[faces[j][2]][1] - verts[faces[j][1]][1], verts[faces[j][2]][2] - verts[faces[j][1]][2]));
  }
  return edge_length/(3.*nfaces);
}


void compute_centroids(const faces_t& faces, const verts_t& verts, verts_t& centroids){
  int nt = faces.shape()[0];
  for (int j=0; j<nt; j++){
    int f0 = faces[j][0];
    int f1 = faces[j][1];
    int f2 = faces[j][2];
    for (int di=0; di<3; di++){
        centroids[j][di] = (verts[f0][di] + verts[f1][di] + verts[f2][di])/3.;

    }
  }
}





std::vector< std::vector<int>> make_neighbour_faces_of_vertex(const verts_t& verts, const faces_t& faces){
  int nt = faces.shape()[0];
  int vt = verts.shape()[0];
  std::vector< std::vector<int>> neighbour_faces_of_vertex;
  for (int fi=0; fi< vt; fi++){
    neighbour_faces_of_vertex.push_back(std::vector<int>());
  }
  for (int fi=0; fi< nt; fi++){
    for (int vi=0; vi<3; vi++){
      int v1 = faces[fi][vi];
      neighbour_faces_of_vertex[v1].push_back(fi);
    }
  }

  return neighbour_faces_of_vertex;
}


void compute_centroid_gradient(const verts_t& centroids, verts_t& centroid_normals_normalized, implicit_function* gradou){

  gradou->eval_gradient(centroids, &centroid_normals_normalized);
    for(int i = 0; i < centroid_normals_normalized.shape()[0]; i++){
      REAL norm = norm_2(centroid_normals_normalized[i][0], centroid_normals_normalized[i][1], centroid_normals_normalized[i][2]);
      assert(norm!=0.);
      for(int j = 0; j < 3; j++){
        centroid_normals_normalized[i][j]=centroid_normals_normalized[i][j]/norm;
        assert(centroid_normals_normalized[i][j] <= 1.);
        assert(centroid_normals_normalized[i][j] >= -1.);
      }
    }
}


void centroids_projection(implicit_function* object, verts_t& verts, const faces_t& faces){

  REAL average_edge;
  average_edge = compute_average_edge_length(verts, faces);

  verts_t centroids;
  compute_centroids(faces, verts, centroids);

  set_centers_on_surface(object, &centroids, average_edge);

  std::vector< std::vector<int>> vertex_neighbours_list;
  vertex_neighbours_list = make_neighbour_faces_of_vertex(verts, faces);

  verts_t centroid_gradients;
  compute_centroid_gradients(centroids, centroid_gradients, object);

  vertex_apply_qem(&verts, faces, centroids, vertex_neighbours_list, centroid_gradients);


}
