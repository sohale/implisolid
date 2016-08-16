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

  // implicit function of the two arrays
  boost::array<int, 1> v1_shape = {n};
  // vectorized_scalar v1(v1_shape);
  // vectorized_scalar v2(v1_shape);

  // object->eval_implicit(x1_arr, &v1);
  // object->eval_implicit(x2_arr, &v2);

  boost::multi_array<int, 1> active_indices(v1_shape);

  for (int i=0; i< n; i++){
    active_indices[i] = i;
  }

  int active_count = n;
  int solved_count = 0;

  boost::array<int, 2> x_mid_shape = {n,3};
  boost::multi_array<REAL, 2> x_mid(x_mid_shape);

  vectorized_scalar v_mid(v1_shape); // implicit function for x_mid
  vectorized_scalar abs_(v1_shape); // absolute value of the implicit function

  // array of indices
  boost::multi_array<int, 1> indices_boundary(v1_shape);
  boost::multi_array<int, 1> indices_outside(v1_shape);
  boost::multi_array<int, 1> indices_inside(v1_shape);
  boost::multi_array<int, 1> indices_eitherside(v1_shape);
  boost::multi_array<int, 1> which_zeroed(v1_shape);

  int iteration = 1;

 // loop
  while (true){
    //mean of (x1,x2)
    for (int i=0; i<active_count; i++){
      x_mid[i][0] = (x1_arr[i][0]+ x2_arr[i][0])/2.;
      x_mid[i][1] = (x1_arr[i][1]+ x2_arr[i][1])/2.;
      x_mid[i][2] = (x1_arr[i][2]+ x2_arr[i][2])/2.;
    }

    object->eval_implicit(x_mid, &v_mid);


    // imcrementing the size of the indices arrays
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

    for (int i=0; i<i_b; i++){
      which_zeroed[i] = (active_indices[indices_boundary[i]]);
    }

    int found_count = i_b;
    solved_count += found_count;

    for (int i=0; i<i_b; i++){
      res_x_arr[which_zeroed[i]][0] = x_mid[indices_boundary[i]][0];
      res_x_arr[which_zeroed[i]][1] = x_mid[indices_boundary[i]][1];
      res_x_arr[which_zeroed[i]][2] = x_mid[indices_boundary[i]][2];

    }


    // changing the values of x2 and x1
    for (int i=0; i<i_i; i++){
      // v2[indices_inside[i]] = v_mid[indices_inside[i]];
      x2_arr[indices_inside[i]][0] = x_mid[indices_inside[i]][0];
      x2_arr[indices_inside[i]][1] = x_mid[indices_inside[i]][1];
      x2_arr[indices_inside[i]][2] = x_mid[indices_inside[i]][2];
    }

    for (int i=0; i<i_o; i++){
      // v1[indices_outside[i]] = v_mid[indices_outside[i]];
      x1_arr[indices_outside[i]][0] = x_mid[indices_outside[i]][0];
      x1_arr[indices_outside[i]][1] = x_mid[indices_outside[i]][1];
      x1_arr[indices_outside[i]][2] = x_mid[indices_outside[i]][2];
    }

    //next round

    for (int i=0; i<i_e; i++){
      active_indices[i] = active_indices[indices_eitherside[i]];
    }

    indices_boundary.resize(boost::extents[i_e]);
    indices_eitherside.resize(boost::extents[i_e]);
    indices_outside.resize(boost::extents[i_e]);
    indices_inside.resize(boost::extents[i_e]);
    which_zeroed.resize(boost::extents[i_e]);
    active_indices.resize(boost::extents[i_e]);

    active_count = active_count - found_count;

    iteration += 1;

    for(int i=0; i<active_count; i++){
      // v1[i] = v1[indices_eitherside[i]];
      // v2[i] = v2[indices_eitherside[i]];

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
  // intilization, objects creation
  REAL min_gradient_len = 0.000001;
  int max_iter = 20;

  REAL max_dist = average_edge;
  int n = centroids.shape()[0];
  boost::array<int, 1> scalar_shape = {n};

  vectorized_scalar fc_a(scalar_shape);
  vectorized_scalar signs_c(scalar_shape);

  object->eval_implicit(centroids, &fc_a);

  boost::array<int, 2> g_a_shape = { n , 3 };
  boost::multi_array<REAL, 2> g_a(g_a_shape);
  // applying the centroid gradient
  compute_centroid_gradient(centroids, g_a, object);

  boost::array<int, 2> vector_shape = {n,3};

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

  vectorized_scalar alpha_list(scalar_shape);

  int iter = 0;
  while(step_size > 0.001){

    step_size = step_size*0.5;
    int max_step;
    max_step = min(max_iter, int(floor(max_dist/ABS(step_size)+0.001)));

    for (int i=1; i< max_step+1; i+=2){
      REAL alpha = float(i)*step_size;
      alpha_list[i + iter] = alpha/average_edge;
      alpha_list[i + iter +1] = -alpha/average_edge;
    }
    iter += max_step;
  }

  alpha_list.resize(boost::extents[iter]);

  //THE algorithm

 //array definition
  boost::multi_array<REAL, 2> best_result_x(vector_shape);
  boost::multi_array<REAL, 2> x1_half(vector_shape);
  boost::multi_array<REAL, 2> xa4(vector_shape);

  vectorized_scalar f_a(scalar_shape);
  vectorized_scalar signs_a(scalar_shape);

  // boolean
  boost::multi_array<bool_t, 1> success0(scalar_shape);
  boost::multi_array<bool_t, 1> success(scalar_shape);

  // indices arrays
  boost::multi_array<int, 1> active_indices(scalar_shape);
  boost::multi_array<int, 1> still_nonsuccess_indices(scalar_shape);


  for (int i=0; i<n; i++){
    active_indices[i] = i;
  }

  int active_count = n;

  still_nonsuccess_indices = active_indices;

  boost::multi_array<bool_t, 1> already_success(scalar_shape);
  boost::multi_array<bool_t, 1> new_success_indices(scalar_shape);


  for (int i=0; i<n; i++){
    already_success[i] = b_false;
  }

  int counter = -1;
  int s_n_s = 0;
  // main part of the algor
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

    for (int j=0; j<active_indices.shape()[0]; j++){
      if (f_a[j] > ROOT_TOLERANCE ){
        signs_a[j] = +1.;
      }
      else if(f_a[j] < -ROOT_TOLERANCE){
        signs_a[j] = -1.;
      }
      else{
        signs_a[j] = 0.;
      }

      success[j] = b_false;

      if(signs_a[j] * signs_c[active_indices[j]] <=0){
        success0[j] = b_true;
      }
      else{
        success0[j] = b_false;
      }

    }

    for (int j=0; j< active_indices.shape()[0]; j++){
      success[active_indices[j]] = success0[j];
    }


    int n_s = 0;
    s_n_s = 0;
    for (int j=0; j<active_indices.shape()[0]; j++){
      if (success[j] == b_true && already_success[j] == b_false){
        new_success_indices[n_s] = j;
        n_s ++;
      }
      else{
        still_nonsuccess_indices[s_n_s] = j;
        s_n_s ++;
      }
    }

    for (int j=0; j< n_s; j++){
      best_result_x[new_success_indices[j]][0] = x1_half[new_success_indices[j]][0];
      best_result_x[new_success_indices[j]][1] = x1_half[new_success_indices[j]][1];
      best_result_x[new_success_indices[j]][2] = x1_half[new_success_indices[j]][2];
    }

    for (int j=0; j<n; j++){
      if(success[j] == b_true){
        already_success[j] = b_true;
      }
    }

    if (s_n_s == 0){
      break;
    }

    active_indices.resize(boost::extents[s_n_s]);
    still_nonsuccess_indices.resize(boost::extents[s_n_s]);
  }

  for (int i=0; i<s_n_s; i++){
    best_result_x[still_nonsuccess_indices[i]][0] = centroids[still_nonsuccess_indices[i]][0];
    best_result_x[still_nonsuccess_indices[i]][1] = centroids[still_nonsuccess_indices[i]][1];
    best_result_x[still_nonsuccess_indices[i]][2] = centroids[still_nonsuccess_indices[i]][2];
  }

  boost::multi_array<REAL, 2> xa1(vector_shape);
  boost::multi_array<REAL, 2> xa2(vector_shape);

  vectorized_scalar f1(scalar_shape);
  vectorized_scalar f2(scalar_shape);
  f1 = fc_a;

  xa1 = centroids;
  xa2 = best_result_x;

  object->eval_implicit(xa2, &f2);

  boost::multi_array<bool_t, 1> zeros2_bool(scalar_shape);
  boost::multi_array<bool_t, 1> zeros1_bool(scalar_shape);
  boost::multi_array<bool_t, 1> zeros1or2(scalar_shape);
  boost::multi_array<int, 1> relevants_bool(scalar_shape);

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

    if (already_success[i] == b_true && zeros1or2[i] == b_false ){
      relevants_bool[r_b] = i;
      r_b ++;
    }

  }
  relevants_bool.resize(boost::extents[r_b]);

  int m = relevants_bool.shape()[0];

  boost::array<int, 2> x1_relevant_shape = {m,3};
  boost::array<int, 1> f1_relevant_shape = {m};
  boost::multi_array<REAL, 2> x1_relevant(x1_relevant_shape);
  boost::multi_array<REAL, 2> x2_relevant(x1_relevant_shape);

  vectorized_scalar f2_relevants(f1_relevant_shape);

  for (int i=0; i<m; i++){
    x1_relevant[i][0] = centroids[relevants_bool[i]][0];
    x1_relevant[i][1] = centroids[relevants_bool[i]][1];
    x1_relevant[i][2] = centroids[relevants_bool[i]][2];

    x2_relevant[i][0] = best_result_x[relevants_bool[i]][0];
    x2_relevant[i][1] = best_result_x[relevants_bool[i]][1];
    x2_relevant[i][2] = best_result_x[relevants_bool[i]][2];
  }

  object->eval_implicit(x2_relevant, &f2_relevants);

  REAL temp0;
  REAL temp1;
  REAL temp2;

  for (int i=0; i<m; i++){
    if (f2_relevants[i] < -ROOT_TOLERANCE){
      temp0 = x2_relevant[i][0];
      temp1 = x2_relevant[i][1];
      temp2 = x2_relevant[i][2];
      x2_relevant[i][0] = x1_relevant[i][0];
      x2_relevant[i][1] = x1_relevant[i][1];
      x2_relevant[i][2] = x1_relevant[i][2];
      x1_relevant[i][0] = temp0;
      x1_relevant[i][1] = temp1;
      x1_relevant[i][2] = temp2;
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

//get the matrix A and b used in vertex_apply_qem
void get_A_b(const std::vector<int> nai, const verts_t& centroids, const verts_t& centroid_gradients, verts_t* A, vectorized_scalar* b){

    int m = nai.size();

    boost::array<int, 2> center_array_shape = {m,3};
    boost::multi_array<REAL, 2> center_array(center_array_shape);
    boost::multi_array<REAL, 2> normals(center_array_shape);

    for (int i=0; i<m; i++){

      normals[i][0] = centroid_gradients[nai[i]][0];
      normals[i][1] = centroid_gradients[nai[i]][1];
      normals[i][2] = centroid_gradients[nai[i]][2];

      center_array[i][0] = centroids[nai[i]][0];
      center_array[i][1] = centroids[nai[i]][1];
      center_array[i][2] = centroids[nai[i]][2];
    }

    for (int j=0; j<m; j++){

      REAL a00 = normals[j][0]*normals[j][0];
      REAL a01 = normals[j][0]*normals[j][1];
      REAL a02 = normals[j][0]*normals[j][2];
      REAL a11 = normals[j][1]*normals[j][1];
      REAL a12 = normals[j][1]*normals[j][2];
      REAL a22 = normals[j][2]*normals[j][2];

      (*A)[0][0] += a00;
      (*A)[0][1] += a01;
      (*A)[0][2] += a02;
      (*A)[1][1] += a11;
      (*A)[1][2] += a12;
      (*A)[2][2] += a22;

      (*b)[0] += a00*center_array[j][0] + a01*center_array[j][1] + a02*center_array[j][2];
      (*b)[1] += a01*center_array[j][0] + a11*center_array[j][1] + a12*center_array[j][2];
      (*b)[2] += a02*center_array[j][0] + a12*center_array[j][1] + a22*center_array[j][2];

    }

    (*A)[1][0] = (*A)[0][1];
    (*A)[2][0] = (*A)[0][2];
    (*A)[2][1] = (*A)[1][2];

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


void vertex_apply_qem(verts_t* verts, faces_t& faces, verts_t& centroids, std::vector< std::vector<int>>& vertex_neighbours_list, verts_t& centroid_gradients){

  int nverts = verts->shape()[0];


  boost::array<int, 2> A_shape = { 3 , 3 };

  verts_t A(A_shape);
  verts_t u(A_shape);
  verts_t s(A_shape);
  verts_t v(A_shape);

  boost::array<int, 2> b_shape = { 3 , 1 };
  vectorized_scalar b(b_shape);
  vectorized_scalar y(b_shape);
  vectorized_scalar utb(b_shape);
  vectorized_scalar new_x(b_shape);

  for (int vi=0; vi<nverts; vi++){

    std::vector<int> nlist;
    for (int i=0; i< vertex_neighbours_list[vi].size(); i++){
      nlist.push_back(vertex_neighbours_list[vi][i]);
    }

    get_A_b(nlist, centroids, centroid_gradients, &A, &b);
    SVD(A, u, s, v); // the SVD

    REAL tau = 680.;
    int rank = 0;
    if (s[1][1]/s[0][0] < 1./tau){
      s[1][1] = 0.;
    }
    else{
      rank ++;
    }

    if (s[2][2]/s[0][0] < 1./tau){
      s[2][2] = 0.;
    }
    else{
      rank ++;
    }

    y[0] = v[0][0]*(*verts)[vi][0] + v[0][1]*(*verts)[vi][1] + v[0][2]*(*verts)[vi][2];
    y[1] = v[1][0]*(*verts)[vi][0] + v[1][1]*(*verts)[vi][1] + v[1][2]*(*verts)[vi][2];
    y[2] = v[2][0]*(*verts)[vi][0] + v[2][1]*(*verts)[vi][1] + v[2][2]*(*verts)[vi][2];

    utb[0] = - u[0][0]*b[0] - u[1][0]*b[1] - u[2][0]*b[2];
    utb[1] = - u[0][1]*b[0] - u[1][1]*b[1] - u[2][1]*b[2];
    utb[2] = - u[0][2]*b[0] - u[1][2]*b[1] - u[2][2]*b[2];

    for (int i=0; i<rank; i++){
      y[i] = utb[i]/s[i][i];
    }

    new_x[0] = v[0][0]*y[0] + v[1][0]*y[1] + v[2][0]*y[2];
    new_x[1] = v[0][1]*y[0] + v[1][1]*y[1] + v[2][1]*y[2];
    new_x[2] = v[0][2]*y[0] + v[1][2]*y[1] + v[2][2]*y[2];


    (*verts)[vi][0] = new_x[0];
    (*verts)[vi][1] = new_x[1];
    (*verts)[vi][2] = new_x[2];
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

  boost::array<int, 2> centroids_shape = { (int)result_faces.size()/3 , 3 };
  boost::multi_array<REAL, 2> centroids(centroids_shape);

  compute_centroids(faces, verts, centroids);

  set_centers_on_surface(object, centroids, average_edge);
cout << "set_centers done " << endl;
  std::vector< std::vector<int>> vertex_neighbours_list;
  vertex_neighbours_list = make_neighbour_faces_of_vertex(verts, faces);

  boost::multi_array<REAL, 2> centroid_gradients(centroids_shape);

  compute_centroid_gradient(centroids, centroid_gradients, object);
cout << "centroids gradient computed" << endl;
  vertex_apply_qem(&verts, faces, centroids, vertex_neighbours_list, centroid_gradients);

  for (int i=0; i<verts.shape()[0]; i++) {
    result_verts[i*3+0] = verts[i][0];
    result_verts[i*3+1] = verts[i][1];
    result_verts[i*3+2] = verts[i][2];
  }

}
