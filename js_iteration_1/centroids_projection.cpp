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

REAL compute_average_edge_length(const faces_t& faces, const verts_t& verts){
  int nfaces = faces.shape()[0];
  REAL edge_length;
  for (int j=0; j<nfaces; j++){
    edge_length += norm_2(verts[faces[j][0]][0] - verts[faces[j][1]][0], verts[faces[j][0]][1] - verts[faces[j][1]][1], verts[faces[j][0]][2] - verts[faces[j][1]][2]);
    edge_length += norm_2(verts[faces[j][0]][0] - verts[faces[j][2]][0], verts[faces[j][0]][1] - verts[faces[j][2]][1], verts[faces[j][0]][2] - verts[faces[j][2]][2]);
    edge_length += norm_2(verts[faces[j][2]][0] - verts[faces[j][1]][0], verts[faces[j][2]][1] - verts[faces[j][1]][1], verts[faces[j][2]][2] - verts[faces[j][1]][2]);
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


void bisection(mp5_implicit::implicit_function* object, verts_t& res_x_arr, verts_t& x1_arr, verts_t& x2_arr, REAL ROOT_TOLERANCE){
  // initilization step
  int n = x1_arr.shape()[0];

  vectorized_scalar v1;
  vectorized_scalar v2;

  object->eval_implicit(x1_arr, &v1);
  object->eval_implicit(x2_arr, &v2);

  // boost::array<int, 1> active_indices_shape = {n};
  // boost::array<int, 1> active_indices= { active_indices_shape };
  //
  // for (int i=0; i<n; i++){
  //   active_indices[i] = i;
  // }

  std::vector<int> active_indices;

  for (int i=0; i<n; i++){
    active_indices.push_back(i);
  }

  int active_count = n;
  int solved_count = 0;

  boost::array<int, 3> x_mid_shape = {n,3};
  boost::multi_array<REAL, 2> x_mid(x_mid_shape);

  vectorized_scalar v_mid;
  vectorized_scalar abs_;

  std::vector<int> indices_boundary;
  std::vector<int> indices_outside;
  std::vector<int> indices_inside;
  std::vector<int> indices_eitherside;
  std::vector<int> which_zeroed;

  int iteration = 1;
 // loop
  while (true){

    for (int i=0; i<active_count; i++){
      x_mid[i][0] = (x1_arr[i][0]+ x2_arr[i][0])/2.;
      x_mid[i][1] = (x1_arr[i][1]+ x2_arr[i][1])/2.;
      x_mid[i][2] = (x1_arr[i][2]+ x2_arr[i][2])/2.;
    }

    object->eval_implicit(x_mid, &v_mid);

    for (int i=0; i<active_count; i++){
      abs_[i] = ABS(v_mid[i]);
      if(abs_[i] <= ROOT_TOLERANCE){
        indices_boundary.push_back(i);
      }
      else{
        indices_eitherside.push_back(i);
      }
      if(v_mid[i] < - ROOT_TOLERANCE){
        indices_outside.push_back(i);
      }
      if(v_mid[i] > ROOT_TOLERANCE){
        indices_inside.push_back(i);
      }

    }

    for (int i=0; i<indices_boundary.size(); i++){
      which_zeroed.push_back(active_indices[indices_boundary[i]]);
    }

    int found_count = indices_boundary.size();
    solved_count += found_count;

    for (int i=0; i<which_zeroed.size(); i++){
      res_x_arr[which_zeroed[i]][0] = x_mid[indices_boundary[i]][0];
      res_x_arr[which_zeroed[i]][1] = x_mid[indices_boundary[i]][1];
      res_x_arr[which_zeroed[i]][2] = x_mid[indices_boundary[i]][2];
    }

    for (int i=0; i<indices_inside.size(); i++){
      v2[indices_inside[i]] = v_mid[indices_inside[i]];
      x2_arr[indices_inside[i]][0] = x_mid[indices_inside[i]][0];
      x2_arr[indices_inside[i]][1] = x_mid[indices_inside[i]][1];
      x2_arr[indices_inside[i]][2] = x_mid[indices_inside[i]][2];
    }

    for (int i=0; i<indices_outside.size(); i++){
      v1[indices_outside[i]] = v_mid[indices_outside[i]];
      x1_arr[indices_outside[i]][0] = x_mid[indices_outside[i]][0];
      x1_arr[indices_outside[i]][1] = x_mid[indices_outside[i]][1];
      x1_arr[indices_outside[i]][2] = x_mid[indices_outside[i]][2];
    }

    //next round
    active_indices.clear();

    for (int i=0; i<indices_eitherside.size(); i++){
      active_indices.push_back(indices_eitherside[i]);
    }

    active_count = active_count - found_count;

    iteration += 1;

    for(int i=0; i<active_count; i++){
      v1[i] = v1[indices_eitherside[i]];

      v2[i] = v2[indices_eitherside[i]];

      x1_arr[i][0] = x1_arr[indices_eitherside[i]][0];
      x1_arr[i][1] = x1_arr[indices_eitherside[i]][1];
      x1_arr[i][2] = x1_arr[indices_eitherside[i]][2];

      x2_arr[i][0] = x2_arr[indices_eitherside[i]][0];
      x2_arr[i][1] = x2_arr[indices_eitherside[i]][1];
      x2_arr[i][2] = x2_arr[indices_eitherside[i]][2];
    }

    if (active_indices.size() == 0){
      break;
    }


    indices_boundary.clear();
    indices_outside.clear();
    indices_inside.clear();
    indices_eitherside.clear();
    which_zeroed.clear();
  }

}



void  set_centers_on_surface(mp5_implicit::implicit_function* object, verts_t& centroids,const REAL average_edge){
  REAL min_gradient_len = 0.000001;
  REAL ROOT_TOLERANCE = 0.0001;
  int max_iter = 20;

  REAL max_dist = average_edge;



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


void centroids_projection(mp5_implicit::implicit_function* object, std::vector<REAL>& result_verts, const std::vector<int>& result_faces){

  boost::array<int, 2> verts_shape = { (int)result_verts.size()/3 , 3 };
  boost::multi_array<REAL, 2> verts(verts_shape);
  boost::array<int, 2> faces_shape = { (int)result_faces.size()/3 , 3 };
  boost::multi_array<int, 2> faces(faces_shape);

  int output_verts=0;
  auto i = result_verts.begin();
  auto e = result_verts.end();
  for(; i!=e; i++, output_verts++){
      verts[output_verts][0] = (*i);
      i++;
      verts[output_verts][1] = (*i);
      i++;
      verts[output_verts][2] = (*i);
  }

  int output_faces=0;
  auto i_f = result_faces.begin();
  auto e_f = result_faces.end();
  for(; i_f!=e_f; i_f++, output_faces++){
      faces[output_faces][0] = (*i_f);
      i_f++;
      faces[output_faces][1] = (*i_f);
      i_f++;
      faces[output_faces][2] = (*i_f);
  }

  REAL average_edge;
  average_edge = compute_average_edge_length(faces,  verts);

  verts_t centroids;
  compute_centroids(faces, verts, centroids);

  set_centers_on_surface(object, centroids, average_edge);

  std::vector< std::vector<int>> vertex_neighbours_list;
  vertex_neighbours_list = make_neighbour_faces_of_vertex(verts, faces);

  verts_t centroid_gradients;
  compute_centroid_gradient(centroids, centroid_gradients, object);

  vertex_apply_qem(&verts, faces, centroids, vertex_neighbours_list, centroid_gradients);


}
