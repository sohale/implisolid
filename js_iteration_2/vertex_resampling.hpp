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
// #include <tuple>
#include <fstream>


#include "../js_iteration_2/implicit_vectorised_algorithms.hpp"

using namespace std;
using namespace mp5_implicit;

using mp5_implicit::compute_centroid_gradient;
using mp5_implicit::compute_centroids;
using mp5_implicit::vectorised_algorithms::norm_2_squared;
using mp5_implicit::vectorised_algorithms::norm_2;

#include "../js_iteration_2/faces_verts_algorithms.hpp"


// typedef boost::multi_array<REAL, 2> verts_t;
// typedef boost::multi_array<int, 2> faces_t;
// typedef std::vector<int> vector_int;
// typedef std::vector<std::vector<int>> neighbour;
// typedef pair<verts_t, faces_t> vf_t;


#include "../js_iteration_2/mesh_algorithms.hpp"

using mp5_implicit::build_faces_of_faces;
using mp5_implicit::make_edge_lookup;

// for debugging purposes
#include "pointset_set.hpp"


inline REAL kij(int i, int j, const vectorized_vect& centroids, const vectorized_vect& centroid_normals_normalized){
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
  REAL kij = std::acos(mimj) / pipj;
  return kij;
}

REAL wi(int i, const faces_of_xxx_type& faces_of_faces, const vectorized_vect& centroids, const vectorized_vect& centroid_normals_normalized, float c){
    //clog << "wi1" << std::endl;  // called 2696 times

  assert(faces_of_faces.shape()[1] == 3);  // turns dyanmic (runtime) into constexpr (literal 3)

  REAL ki = 0;
  for (int j_faces=0; j_faces<3; j_faces ++){
    ki += kij(i, faces_of_faces[i][j_faces], centroids, centroid_normals_normalized);
  }
  REAL wi = 1.0 + c * ki;
  return wi;

}

void vertex_resampling_VV1(
        vectorized_vect& new_verts,
        const std::vector< std::vector<faceindex_type>>& faceslist_neighbours_of_vertex,
        const faces_of_xxx_type& faces_of_faces,
        const vectorized_vect& centroids,
        const vectorized_vect& centroid_normals_normalized,
        const float c
    ) {
    // clog << "VV1: c=" << c << std::endl;
    // exit(1);
    int nfaces = centroids.shape()[0];

    assert(faces_of_faces.shape()[1] == 3);

    boost::array<int, 2> wi_total_array_shape = {nfaces, 1 };
    boost::multi_array<REAL, 1> wi_total_array(wi_total_array_shape);

    for (int i_faces=0; i_faces<nfaces; i_faces++){
        REAL w = wi(i_faces, faces_of_faces, centroids, centroid_normals_normalized, c);
        wi_total_array[i_faces] = w;
    }
    for (int i=0; i< new_verts.shape()[0]; i++){

        // todo: DONT USE  std::vector<?>
        const std::vector<faceindex_type> & umbrella_faces = faceslist_neighbours_of_vertex[i];
        REAL sum_w = 0;
        for (int j=0; j< umbrella_faces.size(); j++){
            sum_w += wi_total_array[umbrella_faces[j]];
        }

        // todo: DONT USE  std::vector<?>
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
        vectorized_vect& new_verts,
        // input
        const vectorized_faces& faces, vectorized_vect& verts,
        // output
        vectorized_vect& centroids, // overwritten! It's output-only
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


    edges_of_xxx_type  edges_of_faces(edges_of_faces_shape);
    faces_of_xxx_type  faces_of_edges(faces_of_edges_shape);
    faces_of_xxx_type  faces_of_faces(edges_of_faces_shape);

    if (STORE_POINTSETS)
    {
    vectorized_vect ps1 = verts;
    point_set_set.emplace(std::make_pair(std::string("pre_resampling_vertices"), ps1));
    }

    compute_centroids(faces, verts, centroids);

    // boost::array<int, 2> centroid_normals_normalized_shape = { nfaces, 3 };
    vectorized_vect  centroid_normals_normalized(boost::extents[nfaces][3]);

    compute_centroid_gradient(centroids, centroid_normals_normalized, object);
    object->eval_gradient(centroids, &centroid_normals_normalized);
    vectorised_algorithms::normalize_1111(centroid_normals_normalized);

    std::vector< std::vector<faceindex_type>> faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(faces, verts.shape()[0]);

    make_edge_lookup(faces, edges_of_faces, faces_of_edges);

    build_faces_of_faces(edges_of_faces, faces_of_edges, faces_of_faces);

    if (STORE_POINTSETS)
    {
    // vectorized_vect ps1 = centroids;
    // point_set_set.emplace(std::make_pair(std::string("pre_p_centroid"), ps1));
    // vectorized_vect ps1 = new_verts;
    // point_set_set.emplace(std::make_pair(std::string("pre_p_centroid"), ps1));
    }

    vertex_resampling_VV1(new_verts, faceslist_neighbours_of_vertex, faces_of_faces,
        centroids, centroid_normals_normalized, c
    );

    if (STORE_POINTSETS)
    {
    vectorized_vect ps2 = new_verts;
    point_set_set.emplace(std::make_pair(std::string("post_resampling_vertices"), ps2));
    }

    // vectorized_vect ps2 = centroids;
    /*
    for (int i=0;i<ps2.shape()[0]; ++i) {
        ps2[i][0] += (rand01()*2.0-1.0)*0.2;
        ps2[i][1] += (rand01()*2.0-1.0)*0.2;
        ps2[i][2] += (rand01()*2.0-1.0)*0.2;
    }
    */
    // point_set_set.emplace(std::make_pair(std::string("post_p_centroids"), ps2));

}
