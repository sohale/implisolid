// Authors: Marc, Solene, Sohail

#pragma once


#include <iostream>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <cmath>
#include <cassert>
#include <map>
// #include <vector>
#include <string>
#include <tuple>
#include <fstream>



#include "../js_iteration_2/implicit_vectorised_algorithms.hpp"

using namespace std;
using namespace mp5_implicit;

using mp5_implicit::compute_centroid_gradient;
using mp5_implicit::compute_centroids;
using mp5_implicit::vectorised_algorithms::norm_2_squared;
using mp5_implicit::vectorised_algorithms::norm_2;

#include "../js_iteration_2/faces_verts_algorithms.hpp"


typedef boost::multi_array<REAL, 2> verts_t;
typedef boost::multi_array<int, 2> faces_t;
typedef std::vector<int> vector_int;
typedef std::vector<std::vector<int>> neighbour;
typedef pair<verts_t, faces_t> vf_t;





void make_edge_lookup(faces_t faces, faces_t& edges_of_faces, faces_t& faces_of_edges){
  int nfaces = faces.shape()[0];
  cout << "nfaces is : " << nfaces << endl;
  assert(nfaces % 2 == 0);
  int num_edges = nfaces*3./2.;

  long modulo = long(num_edges);
  long lookup_array_size = modulo*num_edges + num_edges;
  map<int, int> eulookup;
  map<int, int>::iterator iter;
  int edge_counter = 0;
  for (int fi=0; fi<nfaces; fi++){
    for (int vj=0; vj<3; vj++){
      assert(fi<nfaces);
      bool new_edge = false;
      int v2j = (vj+1)%3;
      int e1 = faces[fi][vj];
      int e2 = faces[fi][v2j];
      int eu_pair_int;
      if (e2 > e1){
        eu_pair_int = int(e1 + e2*modulo);
      }
      else{
        eu_pair_int = int(e2 + e1*modulo);
      }

      iter = eulookup.find(eu_pair_int);
      if (iter== eulookup.end()){
        new_edge = true;
      }
      if (new_edge){
        int e_id = edge_counter;
        edges_of_faces[fi][vj] = e_id;
        faces_of_edges[e_id][0] = fi;
        assert (vj!= v2j);

        eulookup.insert(pair<int,int>(eu_pair_int,e_id));
        edge_counter ++;
        assert(edge_counter <= num_edges);
      }
      else{
        int e_id = iter->second;
        assert (e_id >= 0);
        edges_of_faces[fi][vj] = e_id;
        faces_of_edges[e_id][1] = fi;
        assert (vj!= v2j);
        int other_fi = faces_of_edges[e_id][0];

      }
    }
  }

}

void build_faces_of_faces(faces_t& edges_of_faces, faces_t& faces_of_edges, faces_t& faces_of_faces){
  for(int face = 0; face < edges_of_faces.shape()[0]; face++){
    for(int edge = 0; edge < 3; edge++){
      if(faces_of_edges[edges_of_faces[face][edge]][0]!=face)
        faces_of_faces[face][edge]=faces_of_edges[edges_of_faces[face][edge]][0];
      else{
        assert(faces_of_edges[edges_of_faces[face][edge]][1] != face);
        faces_of_faces[face][edge]=faces_of_edges[edges_of_faces[face][edge]][1];}
    }
  }
  for(int face = 0; face < faces_of_faces.shape()[0]; face ++){
    for(int faces = 0; faces<3; faces++){
      assert(face==faces_of_faces[faces_of_faces[face][faces]][0] ||
          face==faces_of_faces[faces_of_faces[face][faces]][1] ||
          face==faces_of_faces[faces_of_faces[face][faces]][2]);
    }
  }
}

inline REAL kij(int i, int j, const verts_t& centroids, const verts_t& centroid_normals_normalized){
  assert (i!=j);
  REAL pi_x = centroids[i][0];
  REAL pi_y = centroids[i][1];
  REAL pi_z = centroids[i][2];
  REAL pj_x = centroids[j][0];
  REAL pj_y = centroids[j][1];
  REAL pj_z = centroids[j][2];
  REAL mi_x = centroid_normals_normalized[i][0];
  REAL mi_y = centroid_normals_normalized[i][1];
  REAL mi_z = centroid_normals_normalized[i][2];
  REAL mj_x = centroid_normals_normalized[j][0];
  REAL mj_y = centroid_normals_normalized[j][1];
  REAL mj_z = centroid_normals_normalized[j][2];

  REAL mimj = mi_x*mj_x + mi_y*mj_y + mi_z*mj_z;

  if(mimj > 1.0){
    mimj = 1.0;
  }
  if(mimj < -1.0){
    mimj = -1.0;
  }

  REAL pipj = norm_2(pi_x - pj_x, pi_y - pj_y, pi_z - pj_z);
  if (pipj == 0){
    return 0;
  }
  REAL kij = REAL(std::acos(REAL(mimj)))/pipj;
  return kij;
}

REAL wi(int i, const faces_t& faces_of_faces, const verts_t& centroids, const verts_t& centroid_normals_normalized, float c){
    //clog << "wi1" << std::endl;  // called 2696 times

  REAL ki = 0;
  for (int j_faces=0; j_faces<3; j_faces ++){
    ki += kij(i, faces_of_faces[i][j_faces], centroids, centroid_normals_normalized);
  }
  REAL wi = 1.0 + c*ki;
  return wi;

}

void vertex_resampling_VV1(
        verts_t& new_verts,
        const std::vector< std::vector<int>>& faceslist_neighbours_of_vertex, const faces_t& faces_of_faces,
        const verts_t& centroids, const verts_t& centroid_normals_normalized, const float c
    ) {
    clog << "VV1" << std::endl;
    // exit(1);
    int nfaces = centroids.shape()[0];

    boost::array<int, 2> wi_total_array_shape = {nfaces, 1 };
    boost::multi_array<REAL, 1> wi_total_array(wi_total_array_shape);

    for (int i_faces=0; i_faces<nfaces; i_faces++){
        REAL w = wi(i_faces, faces_of_faces, centroids, centroid_normals_normalized, c);
        wi_total_array[i_faces] = w;
    }
    for (int i=0; i< new_verts.shape()[0]; i++){
        const std::vector<int> & umbrella_faces = faceslist_neighbours_of_vertex[i];
        REAL sum_w = 0;
        for (int j=0; j< umbrella_faces.size(); j++){
            sum_w += wi_total_array[umbrella_faces[j]];
        }

        // todo: remove w, aoid std::vector<REAL>. e.g. reuse wi_total_array for normalisation.
        std::vector<REAL> w;
        for (int j=0; j< umbrella_faces.size(); j++){
            w.push_back(wi_total_array[umbrella_faces[j]] / sum_w);
        }
        REAL new_verts_x = 0;
        REAL new_verts_y = 0;
        REAL new_verts_z = 0;
        for (int j=0; j< umbrella_faces.size(); j++){
            new_verts_x += w[j]*centroids[umbrella_faces[j]][0];
            new_verts_y += w[j]*centroids[umbrella_faces[j]][1];
            new_verts_z += w[j]*centroids[umbrella_faces[j]][2];
        }
        new_verts[i][0] = new_verts_x;
        new_verts[i][1] = new_verts_y;
        new_verts[i][2] = new_verts_z;
    }
}

//
/*
  Applies the vertex relaxation given constant c.
  It is the main function for part I (vertex relaxation, i.e. vertex resampling)
    1- Forms the graph internal structures (FoF, FoE, EoF).
    2- calculates the centroids and their gradients.
    3- applies the resampling.
    4- The result is stored in new_verts argument.
*/
void process2_vertex_resampling_relaxation_v1(
        // outputs
        verts_t& new_verts,
        // input
        const faces_t& faces, verts_t& verts,
        // output
        verts_t& centroids, // overwritten! It's output-only
        // inputs
        mp5_implicit::implicit_function* object, float c
    )
{

    int nfaces = faces.shape()[0];
    assert(nfaces % 2 == 0);
    int num_edges = nfaces*3./2.;
    boost::array<int, 2> edges_of_faces_shape = { nfaces, 3 };
    boost::array<int, 2> faces_of_edges_shape = { num_edges, 2 };
    boost::array<int, 2> faces_of_faces_shape = { nfaces, 3 };

    boost::multi_array<int, 2>  edges_of_faces(edges_of_faces_shape);
    boost::multi_array<int, 2>  faces_of_edges(faces_of_edges_shape);
    boost::multi_array<int, 2>  faces_of_faces(edges_of_faces_shape);

    compute_centroids(faces, verts, centroids);

    boost::array<int, 2> centroid_normals_normalized_shape = { nfaces, 3 };
    vectorized_vect  centroid_normals_normalized(centroid_normals_normalized_shape);

    // compute_centroid_gradient(centroids, centroid_normals_normalized, object);
    object->eval_gradient(centroids, &centroid_normals_normalized);
    vectorised_algorithms::normalize_1111(centroid_normals_normalized);

    std::vector< std::vector<int>> faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(verts, faces);

    make_edge_lookup(faces, edges_of_faces, faces_of_edges);

    build_faces_of_faces(edges_of_faces, faces_of_edges, faces_of_faces);

    vertex_resampling_VV1(new_verts, faceslist_neighbours_of_vertex, faces_of_faces,
        centroids, centroid_normals_normalized, c
    );
}
