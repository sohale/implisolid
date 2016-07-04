
#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include "../js_iteration_2/primitives.cpp"
#include "vertex_resampling.cpp"
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
    void build_geometry(int resolution, REAL time);
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


void check_state() {
    if(!_state.active) std::cout << "Error: not active.";
}
void check_state_null() {
    if(_state.active)
        std::cout << "Error: should not be active.";
}

void build_geometry(int resolution, REAL time){

    check_state_null();


    bool enableUvs = true;
    bool enableColors = true;

    string name = "super_bowl";
    _state.mc = new MarchingCubes(resolution, enableUvs, enableColors);

    _state.mc -> isolation = 0.0;
      // before we had some amazing meatballs! merde a celui qui le lira!
      REAL real_size = 10;
      // f_argument is made to always be between 0. and 1.
      REAL f_argument = 0.5;

    _state.mc->create_shape(name, real_size, f_argument);

    _state.mc->seal_exterior();

    const callback_t renderCallback;
    _state.mc->render_geometry(renderCallback);

    if(REPORT_STATS){
    int mapctr = 0;
    for (auto& kv_pair: _state.mc->result_e3map){
        if(0)
            std::cout << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
      }
    }

    for (int i=0; i<3; i++){
    _state.mc->vertex_resampling(name, f_argument);
    }

    if(VERBOSE){
        std::cout << resolution << " " << time << std::endl;
        std::cout << _state.mc << std::endl;
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

    if(vcount*3 != ctr)  std::cout << "sizes dont match: " << (float)ctr/3. << " " << vcount << std::endl;
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
    if(fcount*3 != ctr)  std::cout << "sizes dont match: " << (float)ctr/3. << " " << fcount << std::endl;

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
        std::cout << "Error: finish_geometry() before producing the shape()" << std::endl;
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

  // if (writing_test_file){
  //
  // int resolution = 28;
  // REAL time = 0.2;
  // check_state_null();
  // bool enableUvs = true;
  // bool enableColors = true;
  //
  // _state.mc = new MarchingCubes(resolution, enableUvs, enableColors);
  //
  // _state.mc -> isolation = 0.0;
  //   // before we had some amazing meatballs! merde a celui qui le lira!
  //   REAL real_size= 10;
  // string name = "sphere";
  // _state.mc->create_shape(name,real_size,f_argument);
  //
  // _state.mc->seal_exterior();
  //
  // const callback_t renderCallback;
  // _state.mc->render_geometry(renderCallback);
  //
  // if(REPORT_STATS){
  // int mapctr = 0;
  // for (auto& kv_pair: _state.mc->result_e3map){
  //     if(0)
  //         std::cout << " [" << kv_pair.first << ':' << kv_pair.second << ']';
  //     mapctr++;
  //   }
  // }
  //
  //   _state.mc->vertex_resampling(name, f_argument);
  //
  //  }
    std::cout << "main();" << std::endl;

    return 0;
}
