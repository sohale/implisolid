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

// vectorized bisection
void bisection(mp5_implicit::implicit_function* object, verts_t& res_x_arr, verts_t& x1_arr, verts_t& x2_arr, REAL ROOT_TOLERANCE){
  // initilization step
  int n = x1_arr.shape()[0];

  vectorized_scalar v1;
  vectorized_scalar v2;

  object->eval_implicit(x1_arr, &v1);
  object->eval_implicit(x2_arr, &v2);

  boost::multi_array<int, 1> active_indices;

  for (int i=0; i< n; i++){
    active_indices[i] = i;
  }

  int active_count = n;
  int solved_count = 0;

  boost::array<int, 3> x_mid_shape = {n,3};
  boost::multi_array<REAL, 2> x_mid(x_mid_shape);

  vectorized_scalar v_mid;
  vectorized_scalar abs_;

  boost::multi_array<int, 1> indices_boundary;
  boost::multi_array<int, 1> indices_outside;
  boost::multi_array<int, 1> indices_inside;
  boost::multi_array<int, 1> indices_eitherside;
  boost::multi_array<int, 1> which_zeroed;

  int iteration = 1;
 // loop
  while (true){

    for (int i=0; i<active_count; i++){
      x_mid[i][0] = (x1_arr[i][0]+ x2_arr[i][0])/2.;
      x_mid[i][1] = (x1_arr[i][1]+ x2_arr[i][1])/2.;
      x_mid[i][2] = (x1_arr[i][2]+ x2_arr[i][2])/2.;
    }

    object->eval_implicit(x_mid, &v_mid);

    int i_b = 0;
    int i_e = 0;
    int i_i = 0;
    int i_o = 0;
    for (int i=0; i<active_count; i++){
      abs_[i] = ABS(v_mid[i]);
      if(abs_[i] <= ROOT_TOLERANCE){
        indices_boundary[i_b] = i;
        i_b ++;
      }
      else{
        indices_eitherside[i_e] = i;
        i_e ++;
      }
      if(v_mid[i] < - ROOT_TOLERANCE){
        indices_outside[i_o] = i;
        i_o ++;
      }
      if(v_mid[i] > ROOT_TOLERANCE){
        indices_inside[i_i] = i;
        i_i ++;
      }

    }

    for (int i=0; i<indices_boundary.shape()[0]; i++){
      which_zeroed[i] = (active_indices[indices_boundary[i]]);
    }

    int found_count = indices_boundary.shape()[0];
    solved_count += found_count;

    for (int i=0; i<which_zeroed.shape()[0]; i++){
      res_x_arr[which_zeroed[i]][0] = x_mid[indices_boundary[i]][0];
      res_x_arr[which_zeroed[i]][1] = x_mid[indices_boundary[i]][1];
      res_x_arr[which_zeroed[i]][2] = x_mid[indices_boundary[i]][2];
    }

    for (int i=0; i<indices_inside.shape()[0]; i++){
      v2[indices_inside[i]] = v_mid[indices_inside[i]];
      x2_arr[indices_inside[i]][0] = x_mid[indices_inside[i]][0];
      x2_arr[indices_inside[i]][1] = x_mid[indices_inside[i]][1];
      x2_arr[indices_inside[i]][2] = x_mid[indices_inside[i]][2];
    }

    for (int i=0; i<indices_outside.shape()[0]; i++){
      v1[indices_outside[i]] = v_mid[indices_outside[i]];
      x1_arr[indices_outside[i]][0] = x_mid[indices_outside[i]][0];
      x1_arr[indices_outside[i]][1] = x_mid[indices_outside[i]][1];
      x1_arr[indices_outside[i]][2] = x_mid[indices_outside[i]][2];
    }

    //next round

    for (int i=0; i<indices_eitherside.shape()[0]; i++){
      active_indices[i] = indices_eitherside[i];
    }

    indices_boundary.resize(boost::extents[i_b]);
    indices_eitherside.resize(boost::extents[i_e]);
    indices_outside.resize(boost::extents[i_o]);
    indices_inside.resize(boost::extents[i_i]);
    which_zeroed.resize(boost::extents[i_e]);
    active_indices.resize(boost::extents[i_e]);

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

    if (active_indices.shape()[0] == 0){
      break;
    }

  }

}

// main function
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

  boost::array<int, 3> vector_shape = {n,3};
  boost::multi_array<REAL, 2> dx0_c_grad(vector_shape);

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

  vectorized_scalar alpha_list;

  while(step_size > 0.001){
    step_size = step_size*0.5;
    int max_step;
    max_step = min(max_iter, int(floor(max_dist/ABS(step_size)+0.001)));

    for (int i=1; i< max_step+1; i+=2){
      REAL alpha = float(i)*step_size;
      alpha_list[i] = alpha/average_edge;
      alpha_list[i+1] = -alpha/average_edge;
    }
  }

  //THE algorithm
  boost::array<int, 3> best_result_x_shape = {n,3};
  boost::multi_array<REAL, 2> best_result_x(best_result_x_shape);

  boost::array<int, 3> x1_half_shape = {n,3};
  boost::multi_array<REAL, 2> x1_half(x1_half_shape);

  boost::array<int, 3> xa4_shape = {n,3};
  boost::multi_array<REAL, 2> xa4(xa4_shape);

  boost::array<int, 1> bool_shape  = {n};

  vectorized_scalar f_a;
  vectorized_scalar signs_a;

  boost::multi_array<bool_t, 1> success0;
  boost::multi_array<bool_t, 1> success;
  boost::multi_array<int, 1> active_indices;
  boost::multi_array<int, 1> still_nonsuccess_indices;


  for (int i=0; i<n; i++){
    active_indices[i] = i;
  }

  int active_count = n;

  still_nonsuccess_indices = active_indices;

  boost::multi_array<bool_t, 1> already_success;
  boost::multi_array<bool_t, 1> new_success_indices;


  for (int i=0; i<n; i++){
    already_success[i] = b_false;
  }


  int counter = -1;

  for (int i=0; i< alpha_list.shape()[0]; i++){
    counter += 1;
    for (int j=0; j<n; j++){
      x1_half[j][0] = centroids[j][0] + (max_dist*alpha_list[i])*dx0_c_grad[j][0];
      x1_half[j][1] = centroids[j][1] + (max_dist*alpha_list[i])*dx0_c_grad[j][1];
      x1_half[j][2] = centroids[j][2] + (max_dist*alpha_list[i])*dx0_c_grad[j][2];
    }

    active_indices = still_nonsuccess_indices;

    for (int j=0; j<active_indices.shape()[0]; j++){
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

      success[i] = b_false;

      if(signs_a[j] * signs_c[active_indices[j]] <=0){
        success0[i] = b_true;
      }
      else{
        success0[i] = b_false;
      }

    }

    for (int j=0; j< active_indices.shape()[0]; j++){
      success[active_indices[j]] = success0[j];
    }

    still_nonsuccess_indices.resize(boost::extents[0]);

    int n_s = 0;
    int s_n_s = 0;
    for (int j=0; j<n; j++){
      if (success[j] == 1 and already_success[j] == 0){
        new_success_indices[n_s] = j;
        n_s ++;
      }
      else{
        still_nonsuccess_indices[s_n_s] = j;
        s_n_s ++;
      }
    }

    for (int j=0; j< new_success_indices.shape()[0]; j++){
      best_result_x[new_success_indices[j]][0] = x1_half[new_success_indices[j]][0];
      best_result_x[new_success_indices[j]][1] = x1_half[new_success_indices[j]][1];
      best_result_x[new_success_indices[j]][2] = x1_half[new_success_indices[j]][2];
    }

    for (int j=0; j<n; j++){
      if(success[j] == 1 || already_success[j] == 1){
        already_success[j] = 1;
      }
    }

    if (still_nonsuccess_indices.shape()[0] == 0){
      break;
    }

  }

  for (int i=0; i<still_nonsuccess_indices.shape()[0]; i++){
    best_result_x[still_nonsuccess_indices[i]][0] = centroids[still_nonsuccess_indices[i]][0];
    best_result_x[still_nonsuccess_indices[i]][1] = centroids[still_nonsuccess_indices[i]][1];
    best_result_x[still_nonsuccess_indices[i]][2] = centroids[still_nonsuccess_indices[i]][2];
  }

  boost::array<int, 3> xa1_shape = {n,3};
  boost::multi_array<REAL, 2> xa1(xa1_shape);
  boost::multi_array<REAL, 2> xa2(xa1_shape);

  vectorized_scalar f1;
  vectorized_scalar f2;
  f1 = fc_a;

  xa1 = centroids;
  xa2 = best_result_x;

  object->eval_implicit(xa2, &f2);


  boost::multi_array<bool_t, 1> zeros2_bool(bool_shape);
  boost::multi_array<bool_t, 1> zeros1_bool(bool_shape);
  boost::multi_array<bool_t, 1> zeros1or2(bool_shape);
  boost::multi_array<int, 1> relevants_bool;

  int r_b;
  for (int i=0; i<n; i++){

    if (ABS(f2[i])<= ROOT_TOLERANCE){
      zeros2_bool[i] = b_true;
    }
    else{
      zeros2_bool[i] = b_false;
    }

    if (ABS(f1[i])<= ROOT_TOLERANCE){
      zeros1_bool[i] = b_true;
      best_result_x[i][0] = centroids[i][0];
      best_result_x[i][1] = centroids[i][1];
      best_result_x[i][2] = centroids[i][2];
    }
    else{
      zeros1_bool[i] = b_false;
    }

    if (zeros2_bool[i] == b_true || zeros1_bool[i] == b_true ){
      zeros1or2[i] = b_true;
    }
    else{
      zeros1or2[i] = b_false;
    }

    if (already_success[i] == b_true && zeros1or2[i] == b_true ){
      relevants_bool[r_b] = i;
      r_b ++;
    }

  }

  int m = relevants_bool.shape()[0];

  boost::array<int, 3> x1_relevant_shape = {m,3};
  boost::multi_array<REAL, 2> x1_relevant(x1_relevant_shape);
  boost::multi_array<REAL, 2> x2_relevant(x1_relevant_shape);

  vectorized_scalar f1_relevants;
  vectorized_scalar f2_relevants;

  for (int i=0; i<m; i++){
    x1_relevant[i][0] = centroids[relevants_bool[i]][0];
    x1_relevant[i][1] = centroids[relevants_bool[i]][1];
    x1_relevant[i][2] = centroids[relevants_bool[i]][2];

    x2_relevant[i][0] = best_result_x[relevants_bool[i]][0];
    x2_relevant[i][1] = best_result_x[relevants_bool[i]][1];
    x2_relevant[i][2] = best_result_x[relevants_bool[i]][2];
  }

  object->eval_implicit(x1_relevant, &f1_relevants);
  object->eval_implicit(x2_relevant, &f2_relevants);

  verts_t temp;
  for (int i=0; i<m; i++){
    int k = 0;
    if (f2_relevants[i] < -ROOT_TOLERANCE){
      temp[k][0] = x2_relevant[i][0];
      temp[k][1] = x2_relevant[i][1];
      temp[k][2] = x2_relevant[i][2];
      k += 1;
      x2_relevant[i][0] = x1_relevant[i][0];
      x2_relevant[i][1] = x1_relevant[i][1];
      x2_relevant[i][2] = x1_relevant[i][2];
      x1_relevant[i][0] = temp[k][0];
      x1_relevant[i][1] = temp[k][1];
      x1_relevant[i][2] = temp[k][2];
    }


  }

  boost::multi_array<REAL, 2> x_bisect(x1_relevant_shape);
  // calling the vectorized bisection
  bisection(object, x_bisect, x1_relevant, x2_relevant, ROOT_TOLERANCE);

  //changing the values of the centroids
  for (int i=0; i<m; i++){
    centroids[relevants_bool[i]][0] = x_bisect[i][0];
    centroids[relevants_bool[i]][1] = x_bisect[i][1];
    centroids[relevants_bool[i]][2] = x_bisect[i][2];
  }

  for (int i=0; i<n; i++){
    if (zeros1or2[i] == 1){
      centroids[i][0] = best_result_x[i][0];
      centroids[i][1] = best_result_x[i][1];
      centroids[i][2] = best_result_x[i][2];
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
