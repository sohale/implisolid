
#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include "../js_iteration_2/primitives.cpp"
#include "vertex_resampling.cpp"
#include "centroids_projection.cpp"
#include <fstream>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"

using namespace std;
using namespace mp5_implicit;

const bool VERBOSE = false;
const bool REPORT_STATS = false;
bool writing_test_file = false;

typedef unsigned short int dim_t; //small integers for example the size of one side of the grid
typedef float REAL;

typedef  boost::multi_array<REAL, 1>  array1d;
typedef boost::array<array1d::index, 1>  array_shape_t;
typedef array1d::index  index_t;

typedef index_t index3_t; //Range of the element type has to be large enough, larger than (size^3)*3.
typedef boost::multi_array<index3_t, 1>   array1d_e3;
typedef std::map<index3_t,int>  e3map_t;


struct callback_t { void call (void*) const { } callback_t(){} };

#include "mcc2_marching_cubes_MS.hpp"

extern "C" {
    void build_geometry(int resolution, REAL mc_size, REAL time);
    int get_v_size();
    int get_f_size();
    void get_f(int*, int);
    void get_v(REAL*, int);
    void finish_geometry();
    void* get_f_ptr();
    void* get_v_ptr();
};


typedef struct {
    bool active = 0;
    MarchingCubes* mc = 0;
} state_t;

state_t _state;

void vertex_resampling(implicit_function* object, REAL f_argument,  float c, MarchingCubes& mc){

      boost::array<int, 2> verts_shape = { (int)mc.result_verts.size()/3 , 3 };
      boost::multi_array<REAL, 2> verts(verts_shape);

      boost::array<int, 2> faces_shape = { (int)mc.result_faces.size()/3 , 3 };
      boost::multi_array<int, 2> faces(faces_shape);

      boost::multi_array<REAL, 2> centroids (faces_shape);
      boost::multi_array<REAL, 2> new_verts (verts_shape);

      int output_verts=0;
      auto i = mc.result_verts.begin();
      auto e = mc.result_verts.end();
      for(; i!=e; i++, output_verts++){
        verts[output_verts][0] = (*i);
        i++;
        verts[output_verts][1] = (*i);
        i++;
        verts[output_verts][2] = (*i);
      }

      int output_faces=0;
      auto i_f = mc.result_faces.begin();
      auto e_f = mc.result_faces.end();
      for(; i_f!=e_f; i_f++, output_faces++){
        faces[output_faces][0] = (*i_f);
        i_f++;
        faces[output_faces][1] = (*i_f);
        i_f++;
        faces[output_faces][2] = (*i_f);
      }


      if (writing_test_file){

      ofstream f_out("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_algo_cpp.txt");

      f_out << "0ld vertex :" << endl;
      for (int i=0; i< mc.result_verts.size()/3.; i++){
          f_out << mc.result_verts[3*i];
          f_out << " ";
          f_out << mc.result_verts[3*i+1];
          f_out << " ";
          f_out << mc.result_verts[3*i+2];
          f_out <<  "\n";
      }
      f_out << endl;

      process2_vertex_resampling_relaxation(new_verts, faces, verts, centroids, object, f_argument, c);

      for (int i=0; i<verts.shape()[0]; i++){
        mc.result_verts[i*3+0] = new_verts[i][0];
        mc.result_verts[i*3+1] = new_verts[i][1];
        mc.result_verts[i*3+2] = new_verts[i][2];

      }

      f_out << "n3w vertex :" << endl;
      for (int i=0; i< mc.result_verts.size()/3.; i++){
          f_out << mc.result_verts[3*i];
          f_out << " ";
          f_out << mc.result_verts[3*i+1];
          f_out << " ";
          f_out << mc.result_verts[3*i+2];
          f_out <<  "\n";
      }

      f_out << "faces :" << endl;
      for (int i=0; i< mc.result_faces.size()/3.; i++){
          f_out << mc.result_faces[3*i];
          f_out << " ";
          f_out << mc.result_faces[3*i+1];
          f_out << " ";
          f_out << mc.result_faces[3*i+2];
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

    else {
    process2_vertex_resampling_relaxation(new_verts, faces, verts, centroids, object, f_argument, c);

    for (int i=0; i<verts.shape()[0]; i++){
      mc.result_verts[i*3+0] = new_verts[i][0];
      mc.result_verts[i*3+1] = new_verts[i][1];
      mc.result_verts[i*3+2] = new_verts[i][2];

    }

   }

}

void check_state() {
    if(!_state.active) loger << "Error: not active.";
}
void check_state_null() {
    if(_state.active)
        loger << "Error: should not be active.";
}

void build_geometry(int resolution, REAL mc_size, REAL time){

    check_state_null();


    bool enableUvs = true;
    bool enableColors = true;

    string name = "cube";
    _state.mc = new MarchingCubes(resolution, mc_size, enableUvs, enableColors);

    _state.mc -> isolation = 0.0;
      // before we had some amazing meatballs! merde a celui qui le lira !


      //********this should become an input of build geometry (and so be set in the html file)*******
      REAL grid_real_size = 1.5 ;


      // f_argument is made to always be between 0. and 1.
      REAL f_argument = 0.5;

    implicit_function * object;

    if (name == "double_mushroom"){
      double_mushroom mushroom(1.4, 0.3, 0.3, 2, 0.1,0.1,0.1); //3.3
      object = &mushroom;
    }
    else if (name == "egg"){
      egg segg(0.3, 0.4, 0.5, 0.2, 0.1, 0.3);
      object = &segg; // super egg !
    }
    else if (name == "torus"){
      torus Tor(4.,0.2,0.2,0.2);
      object = &Tor; // super egg !
    }
    else if (name == "honey_comb"){
      honey_comb honey_comb(0.5, 30., 1.);
      object = &honey_comb; // super egg !
    }
    else if (name == "sphere"){
      unit_sphere sphere(f_argument, 0.2, 0.1, 0.3);
      object = &sphere;
    }
    else if (name == "cube"){
      cube cube(0.6, 0.6, 0.6);
      boost::array<int, 2> direction_shape = { 1, 3 };
      boost::multi_array<REAL, 2> direction(direction_shape);
      direction[0][0] = 0.;
      direction[0][1] = 0.0;
      direction[0][2] = 0.3;

    //  cube.rotate(2., direction);
      object = &cube;
    }
    else if (name == "super_bowl"){
      super_bowl super_bowl(0.5, 0.2, 0.3, 0.); //0.5
      object = &super_bowl;
    }

    else if (name == "scone"){
      scone scone(0.8, 0.1, 0.5, 0., 0., 0.5);
      object = &scone;
    }

    else if (name == "lego"){
      legoland legoland(0.8, 0.1, 0.5);
      object = &legoland;
    }
    else if (name == "dice"){
      dice dichu(0.8, 0.1, 0.5);
      object = &dichu;
    }
    else if (name == "heart"){
      heart heartou(0.8, 0.1, 0.5);
      object = &heartou;
    }
    else if (name == "scylinder"){
      boost::array<int, 2> direction_shape = { 1, 3 };
      boost::multi_array<REAL, 2> direction(direction_shape);
      direction[0][0] = 0.;
      direction[0][1] = 0.0;
      direction[0][2] = 0.3;
      REAL w[3];
      w[0] = 0;
      w[1] = 1;
      w[2] = 0;
      scylinder scylinder(w , 0.4, 0.4, 0.8, 0., 0.0, 0.0); //0.7
  //    scylinder.rotate(2., direction);
      object = &scylinder;
    }
    else if (name == "egg_cylinder"){

      REAL w[3];
      w[0] = 0;
      w[1] = 1;
      w[2] = 0;
      egg segg(0.6, 0.5, 0.5);
      scylinder scylinder(w, 0.2, 0.1, 0.2);
  //    cube cube(0.4, 0.4, 0.4);
      CrispIntersection crispou(segg, scylinder);
      object = &crispou;
    }
    else if (name == "egg_transform"){
      boost::array<int, 2> direction_shape = { 1, 3 };
      boost::multi_array<REAL, 2> direction(direction_shape);
      direction[0][0] = 0.;
      direction[0][1] = 0.3;
      direction[0][2] = 0.3;

      unit_sphere segg(0.5);
      // segg.move(direction);
      // segg.resize(12.);
      segg.rotate(2., direction);
      object = &segg; // super egg !
    }
    else {
      cout << "Error! You must enter a valid name! So I made a sphere!" << endl;
      unit_sphere sphere(f_argument);
      object = &sphere;
     }


    _state.mc->create_shape(object, grid_real_size);

    _state.mc->seal_exterior();

    const callback_t renderCallback;
    _state.mc->render_geometry(renderCallback);

    if(REPORT_STATS){
    int mapctr = 0;
    for (auto& kv_pair: _state.mc->result_e3map){
        if(0)
            loger << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
      }
    }

    float c=2000.;
    for (int i=0; i<3; i++){
     vertex_resampling(object, f_argument, c, *(_state.mc));
    }

    centroids_projection(object, _state.mc->result_verts, _state.mc->result_faces);

    if(VERBOSE){
        loger << resolution << " " << time << std::endl;
        loger << _state.mc << std::endl;
    }
    _state.active = true;

    check_state();
}
int get_f_size() {
    check_state();
    return _state.mc->result_faces.size()/3;
}
int get_v_size(){
    check_state();
    return _state.mc->result_verts.size()/3;
}
void get_v(REAL* v_out, int vcount){
    check_state();

    // Vertices
    int ctr = 0;
    for(std::vector<REAL>::iterator it=_state.mc->result_verts.begin(); it < _state.mc->result_verts.end(); it+=3 ){
        for(int di=0; di<3; di++){
            v_out[ctr] = *( it + di );
            ctr++;
        }
    }

    if(vcount*3 != ctr)  loger << "sizes dont match: " << (float)ctr/3. << " " << vcount << std::endl;
}

void get_f(int* f_out, int fcount){
    check_state();

    int ctr = 0;
    for(std::vector<int>::iterator it=_state.mc->result_faces.begin(); it < _state.mc->result_faces.end(); it+=3 ){
        for(int di=0; di<3; di++){
            f_out[ctr] = *( it + di );

            ctr++;
        }
    }
    if(fcount*3 != ctr)  loger << "sizes dont match: " << (float)ctr/3. << " " << fcount << std::endl;

};

void* get_v_ptr(){
    check_state();
    return (void*)(_state.mc->result_verts.data());
}

void* get_f_ptr(){
    check_state();
    return (void*)(_state.mc->result_faces.data());
}


void finish_geometry() {
    check_state();
    if(_state.mc == 0){
        loger << "Error: finish_geometry() before producing the shape()" << std::endl;
    }
    if(!_state.active){

    }
    else{
    }
    delete _state.mc;
    _state.active = false;
    _state.mc = 0;
};


#include "timer.hpp"

int main() {
  //***** This part was used to creat test files *****//
  //
  // if (writing_test_file){
  //
  // int resolution = 28;
  // REAL time = 0.2;
  // check_state_null();
  // bool enableUvs = true;
  // bool enableColors = true;
  //
  // _state.mc = new MarchingCubes(resolution, 1.0, enableUvs, enableColors);
  //
  // _state.mc -> isolation = 0.0;
  //   // before we had some amazing meatballs! merde a celui qui le lira!
  //   REAL grid_real_size= 10;
  // string name = "sphere";
  // _state.mc->create_shape(name,grid_real_size,f_argument);
  //
  // _state.mc->seal_exterior();
  //
  // const callback_t renderCallback;
  // _state.mc->render_geometry(renderCallback);
  //
  //   _state.mc->vertex_resampling(name, f_argument);
  //
  //  }
  //   loger << "main();" << std::endl;

    return 0;
}
