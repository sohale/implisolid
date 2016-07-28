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


void  set_centers_on_surface(mp5_implicit::implicit_function* object, verts_t& centroids, const REAL average_edge){
  REAL min_gradient_len = 0.000001;
  int max_iter = 20;

  REAL max_dist = average_edge;

  vectorized_scalar fc_a;
  vectorized_scalar signs_c;

  int n = fc_a.shape()[0];

  object->eval_implicit(centroids, &fc_a);

  verts_t g_a;
  compute_centroid_gradient(centroids, g_a, object);

  boost::array<int, 3> dx0_c_grad_shape = {n,3};
  boost::multi_array<REAL, 2> dx0_c_grad(dx0_c_grad_shape);

  for (int i=0; i<fc_a.shape()[0]; i++){
    if (fc_a[i] > ROOT_TOLERANCE ){
      signs_c[i] = +1.;
    }
    else if(fc_a[i] < -ROOT_TOLERANCE){
      signs_c[i] = -1.;
    }
    else{
      signs_c[i] = 0.;
    }

    dx0_c_grad[i][0] = g_a[i][0]*signs_c[i];
    dx0_c_grad[i][1] = g_a[i][1]*signs_c[i];
    dx0_c_grad[i][2] = g_a[i][2]*signs_c[i];
  }

  REAL step_size = max_dist;

  std::vector<int> alpha_list;

  while(step_size > 0.001){
    step_size = step_size*0.5;
    int max_step;
    max_step = min(max_iter, int(floor(max_dist/ABS(step_size)+0.001)));

    for (int i=1; i< max_step+1; i+=2){
      REAL alpha = float(i)*step_size;
      alpha_list.push_back(alpha/average_edge);
      alpha_list.push_back(-alpha/average_edge);
    }
  }

  //THE algorithm
  boost::array<int, 3> best_result_x_shape = {n,3};
  boost::multi_array<REAL, 2> best_result_x(best_result_x_shape);

  boost::array<int, 3> x1_half_shape = {n,3};
  boost::multi_array<REAL, 2> x1_half(x1_half_shape);

  boost::array<int, 3> xa4_shape = {n,3};
  boost::multi_array<REAL, 2> xa4(xa4_shape);

  vectorized_scalar f_a;
  vectorized_scalar signs_a;

  std::vector<int> success0;
  std::vector<int> success;

  std::vector<int> active_indices;

  for (int i=0; i<n; i++){
    active_indices.push_back(i);
  }

  int active_count = n;
  std::vector<int> still_nonsuccess_indices = active_indices;

  std::vector<int> new_success_indices;
  std::vector<int> already_success;

  for (int i=0; i<n; i++){
    already_success.push_back(0);
  }


  int counter = -1;

  for (int i=0; i< alpha_list.size(); i++){
    counter += 1;
    for (int j=0; j<n; j++){
      x1_half[j][0] = centroids[j][0] + (max_dist*alpha_list[i])*dx0_c_grad[j][0];
      x1_half[j][1] = centroids[j][1] + (max_dist*alpha_list[i])*dx0_c_grad[j][1];
      x1_half[j][2] = centroids[j][2] + (max_dist*alpha_list[i])*dx0_c_grad[j][2];
    }

    active_indices = still_nonsuccess_indices;

    for (int j=0; j<active_indices.size(); j++){
      xa4[j][0] = x1_half[active_indices[j]][0];
      xa4[j][1] = x1_half[active_indices[j]][1];
      xa4[j][2] = x1_half[active_indices[j]][2];
    }

    object->eval_implicit(xa4, &f_a);

    for (int j=0; j<f_a.shape()[0]; j++){
      if (f_a[j] > ROOT_TOLERANCE ){
        signs_a[j] = +1.;
      }
      else if(f_a[j] < -ROOT_TOLERANCE){
        signs_a[j] = -1.;
      }
      else{
        signs_a[j] = 0.;
      }

      success.push_back(0);

      if(signs_a[j] * signs_c[active_indices[j]] <=0){
        success0.push_back(1);
      }
      else{
        success0.push_back(0);
      }

    }

    for (int j=0; j< active_indices.size(); j++){
      success[active_indices[j]] = success0[j];
    }

    still_nonsuccess_indices.clear();

    for (int j=0; j<n; j++){
      if (success[j] == 1 and already_success[j] == 0){
        new_success_indices.push_back(j);
      }
      else{
        still_nonsuccess_indices.push_back(j);
      }
    }

    for (int j=0; j< new_success_indices.size(); j++){
      best_result_x[new_success_indices[j]][0] = x1_half[new_success_indices[j]][0];
      best_result_x[new_success_indices[j]][1] = x1_half[new_success_indices[j]][1];
      best_result_x[new_success_indices[j]][2] = x1_half[new_success_indices[j]][2];
    }

    for (int j=0; j<n; j++){
      if(success[j] == 1 or already_success[j] == 1){
        already_success[j] = 1;
      }
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


void vertex_apply_qem(verts_t* verts, faces_t& faces, verts_t& centroids, std::vector< std::vector<int>> vertex_neighbours_list, verts_t& centroid_gradients){

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
