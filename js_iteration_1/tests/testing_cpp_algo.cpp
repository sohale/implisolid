/*
   author ; Solene Chauvier & Marc Fraysse
   This program generates output textfile to be compared with the output of the equivalent python program.
*/

#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include "../js_iteration_2/primitives.hpp"
#include "../js_iteration_2/vertex_resampling.hpp"
#include <fstream>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "../mcc2_MS.cpp"

using namespace std;
using namespace mp5_implicit;

const bool VERBOSE = false;
const bool REPORT_STATS = false;
bool writing_test_file = true;

typedef unsigned short int dim_t; //small integers for example the size of one side of the grid
typedef float REAL;

typedef  boost::multi_array<REAL, 1>  array1d;
typedef boost::array<array1d::index, 1>  array_shape_t;
typedef array1d::index  index_t;

typedef index_t index3_t; //Range of the element type has to be large enough, larger than (size^3)*3.
typedef boost::multi_array<index3_t, 1>   array1d_e3;
//typedef std::map<index3_t,int>  e3map_t;


#include "mcc2_marching_cubes_MS.hpp"

#include "timer.hpp"

int main() {
  //***** This part was used to creat test files *****//

  if (writing_test_file){

  int resolution = 28;
  REAL time = 0.2;
  REAL f_argument = 0.5;
  mp5_implicit::implicit_function * object;

  check_state_null();
  bool enableUvs = true;
  bool enableColors = true;

  _state.mc = new MarchingCubes(resolution, 1.0, enableUvs, enableColors);

  _state.mc -> isolation = 0.0;
    // before we had some amazing meatballs! merde a celui qui le lira!
  REAL grid_real_size= 10;

  unit_sphere sphere(f_argument);
  object = &sphere;
  _state.mc->create_shape(object, grid_real_size);


  _state.mc->seal_exterior();

  //const callback_t renderCallback;
  _state.mc->render_geometry(/*renderCallback*/);

  mcc2_MS::vertex_resampling_v1(object, f_argument, c, *(_state.mc));

   }
    std::clog << "main();" << std::endl;

    return 0;
}
