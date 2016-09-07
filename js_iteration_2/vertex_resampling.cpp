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

#include "implicit_function.hpp"

#include "implicit_vectorised_algorithms.hpp"

using namespace std;
//using namespace mp5_implicit;

using mp5_implicit::compute_centroid_gradient;
using mp5_implicit::compute_centroids;
using mp5_implicit::vectorised_algorithms::norm_2_squared;
using mp5_implicit::vectorised_algorithms::norm_2;

#include "faces_verts_algorithms.hpp"


/*
typedef float REAL;
typedef struct {
   REAL x, y, z;
} XYZ;

REAL ABS(REAL x){
  if(x<0)
    return -x;
  return x;
}
*/


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

REAL kij(int i, int j, verts_t& centroids, verts_t& centroid_normals_normalized){
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
  REAL kij = REAL(acos(REAL(mimj)))/pipj;
  return kij;
}

REAL wi(int i, faces_t& faces_of_faces, verts_t& centroids, verts_t& centroid_normals_normalized, float c){
  REAL ki = 0;
  for (int j_faces=0; j_faces<3; j_faces ++){
    ki += kij(i, faces_of_faces[i][j_faces], centroids, centroid_normals_normalized);
  }
  REAL wi = 1.0 + c*ki;
  return wi;

}

void vertex_resampling(verts_t& new_verts, std::vector< std::vector<int>>& faceslist_neighbours_of_vertex, faces_t& faces_of_faces,
verts_t& centroids, verts_t& centroid_normals_normalized, float c){
  int nfaces = centroids.shape()[0];

  boost::array<int, 2> wi_total_array_shape = {nfaces, 1 };
  boost::multi_array<REAL, 1> wi_total_array(wi_total_array_shape);

  for (int i_faces=0; i_faces<nfaces; i_faces++){
    REAL w = wi(i_faces, faces_of_faces, centroids, centroid_normals_normalized, c);
    wi_total_array[i_faces] = w;
  }
  for (int i=0; i< new_verts.shape()[0]; i++){
    std::vector<int> umbrella_faces = faceslist_neighbours_of_vertex[i];
    std::vector<REAL> w;
    REAL sum_w = 0;

    for (int j=0; j< umbrella_faces.size(); j++){
      sum_w += wi_total_array[umbrella_faces[j]];
    }

    for (int j=0; j< umbrella_faces.size(); j++){
       w.push_back(wi_total_array[umbrella_faces[j]]/sum_w);
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

// main function
void process2_vertex_resampling_relaxation(
        // outputs
        verts_t& new_verts,
        // input
        const faces_t& faces, verts_t& verts,
        // output
        verts_t& centroids,
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

    compute_centroid_gradient(centroids, centroid_normals_normalized, object);

    std::vector< std::vector<int>> faceslist_neighbours_of_vertex = make_neighbour_faces_of_vertex(verts, faces);

    make_edge_lookup(faces, edges_of_faces, faces_of_edges);

    build_faces_of_faces(edges_of_faces, faces_of_edges, faces_of_faces);

    vertex_resampling(new_verts, faceslist_neighbours_of_vertex, faces_of_faces,
        centroids, centroid_normals_normalized, c
    );
}

void writing_test_file(
        string  output_file_name,
        std::vector<REAL>&result_verts, std::vector<int>&result_faces,
        vectorized_vect & new_verts, vectorized_vect & centroids
    )
{

    ofstream f_out(output_file_name);

    f_out << "0ld vertex :" << endl;
    for (int i=0; i< result_verts.size()/3.; i++){
        f_out << result_verts[3*i];
        f_out << " ";
        f_out << result_verts[3*i+1];
        f_out << " ";
        f_out << result_verts[3*i+2];
        f_out <<  "\n";
    }
    f_out << endl;

    /*
    boost::array<int, 2> verts_shape = { (int)result_verts.size()/3 , 3 };
    vectorized_vect  verts(verts_shape);
    boost::array<int, 2> faces_shape = { (int)result_faces.size()/3 , 3 };
    boost::multi_array<int, 2> faces(faces_shape);
    vectorized_vect  new_verts (verts_shape);
    vectorized_vect  centroids (faces_shape);
    process2_vertex_resampling_relaxation(new_verts, faces, verts, centroids, object, c);
    */
    //new_verts, centroids

    for (int i=0; i<new_verts.shape()[0]; i++){
        result_verts[i*3+0] = new_verts[i][0];
        result_verts[i*3+1] = new_verts[i][1];
        result_verts[i*3+2] = new_verts[i][2];
    }

    f_out << "n3w vertex :" << endl;
    for (int i=0; i< result_verts.size()/3.; i++){
        f_out << result_verts[3*i];
        f_out << " ";
        f_out << result_verts[3*i+1];
        f_out << " ";
        f_out << result_verts[3*i+2];
        f_out <<  "\n";
    }

    f_out << "faces :" << endl;
    for (int i=0; i< result_faces.size()/3.; i++){
        f_out << result_faces[3*i];
        f_out << " ";
        f_out << result_faces[3*i+1];
        f_out << " ";
        f_out << result_faces[3*i+2];
        f_out <<  "\n";
    }

    f_out << "centroids:" << endl;
    for (int i=0; i< centroids.shape()[0]; i++){
        f_out << centroids[i][0];
        f_out << " ";
        f_out << centroids[i][1];
        f_out << " ";
        f_out << centroids[i][2];
        f_out <<  "\n";
    }
    f_out.close();
}


REAL rand01() {
    return static_cast<REAL>(rand())/ (static_cast<REAL>( RAND_MAX) + 1 );
}
void vertex_resampling(mp5_implicit::implicit_function* object,  float c,
      std::vector<REAL>&result_verts, std::vector<int>& result_faces
      )
{

    /*
    for (int i=0; i< result_verts.size()/3.; i++){
        result_verts[i*3+0] += (rand01()-0.5 )* 0.3;
    }
    return;
    */
    /*
    boost::array<int, 2> verts_shape = { (int)result_verts.size()/3 , 3 };
    vectorized_vect  verts(verts_shape);

    boost::array<int, 2> faces_shape = { (int)result_faces.size()/3 , 3 };
    boost::multi_array<int, 2> faces(faces_shape);

    //vectorized_vect  centroids (faces_shape);
    //vectorized_vect  new_verts (verts_shape);
    */
    boost::array<int, 2> verts_shape = { (int)result_verts.size()/3 , 3 };
    vectorized_vect  verts(verts_shape);
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


    vectorized_vect  new_verts (verts_shape);
    vectorized_vect  centroids (faces_shape);
    process2_vertex_resampling_relaxation(new_verts, faces, verts, centroids, object, c);

    bool writing_test_file_ = false;
    if (writing_test_file_) {

        writing_test_file(
            "/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_algo_cpp.txt",
            result_verts, result_faces,
            new_verts, centroids
        );

    } else {
        REAL d = 0;
        for (int i=0; i<verts.shape()[0]; i++) {
            // std::clog << "result_verts  = new_verts: " <<result_verts[i*3+0] << " = " << new_verts[i][0] << std::endl;
          for(int j=0;j<3;j++) {
              d += std::abs( result_verts[i*3+j] - new_verts[i][j]);
          }
            result_verts[i*3+0] = new_verts[i][0];
            result_verts[i*3+1] = new_verts[i][1];
            result_verts[i*3+2] = new_verts[i][2];
        }
        std::clog << "<d> = " <<  d/((REAL)(verts.shape()[0]))/3 << std::endl;

   }


}
