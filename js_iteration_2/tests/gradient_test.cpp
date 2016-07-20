// author ; Solene Chauvier & Marc Fraysse

#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "../../js_iteration_1/mcc2_MS.cpp"


using namespace std;
using namespace mp5_implicit;

int main(){

bool enableUvs = true;
bool enableColors = true;
int resolution = 28;

string name = "egg";
REAL mc_size = 1.0;
_state.mc = new MarchingCubes(resolution, mc_size, enableUvs, enableColors);

_state.mc -> isolation = 0.0;
  // before we had some amazing meatballs! merde a celui qui le lira !

  //********this should become an input of build geometry (and so be set in the html file)*******
  REAL grid_real_size = 10.;

  // f_argument is made to always be between 0. and 1.
  REAL f_argument = 0.5;

implicit_function * object;

// if (name == "double_mushroom"){
//   double_mushroom mushroom(1.4, 0.3, 0.3, 2, 0.1,0.1,0.1); //3.3
//   object = &mushroom;
// }
// else if (name == "egg"){
//   egg segg(0.3, 0.4, 0.5, 0.2, 0.1, 0.3);
//   object = &segg; // super egg !
// }
// else if (name == "sphere"){
//   unit_sphere sphere(f_argument, 0.2, 0.1, 0.3);
//   object = &sphere;
// }
// else if (name == "cube"){
  cube cube(0.4, 0.4, 0.4);
  // boost::array<int, 2> direction_shape = { 1, 3 };
  // boost::multi_array<REAL, 2> direction(direction_shape);
  // direction[0][0] = 0.;
  // direction[0][1] = 0.3;
  // direction[0][2] = 0.3;
  //
  // cube.rotate(2., direction);
  object = &cube;
// }
// else if (name == "super_bowl"){
//   super_bowl super_bowl(0.5, 0.2, 0.3, 0.); //0.5
//   object = &super_bowl;
// }
// else if (name == "scone"){
//   scone scone(0.8, 0.3, 0.3, 2, 0.1, 0.1, -0.1);
//   object = &scone;
// }
// else if (name == "scylinder"){
//   REAL w[3];
//   w[0] = 0;
//   w[1] = 1;
//   w[2] = 0;
//   scylinder scylinder(w , 0.4, 0.4, 0.8, 0., 0.0, 0.0); //0.7
//   object = &scylinder;
// }
// else if (name == "egg_cylinder"){
//
//   REAL w[3];
//   w[0] = 0;
//   w[1] = 1;
//   w[2] = 0;
//   egg segg(0.6, 0.5, 0.5);
//   scylinder scylinder(w, 0.2, 0.1, 0.2);
// //    cube cube(0.4, 0.4, 0.4);
//   CrispIntersection crispou(segg, scylinder);
//   object = &crispou;
// }
// else if (name == "egg_transform"){
//   boost::array<int, 2> direction_shape = { 1, 3 };
//   boost::multi_array<REAL, 2> direction(direction_shape);
//   direction[0][0] = 0.;
//   direction[0][1] = 0.3;
//   direction[0][2] = 0.3;
//
//   unit_sphere segg(0.5);
//   // segg.move(direction);
//   // segg.resize(12.);
//   segg.rotate(2., direction);
//   object = &segg; // super egg !
// }
// else {
//   cout << "Error! You must enter a valid name! So I made a sphere!" << endl;
//   unit_sphere sphere(f_argument);
//   object = &sphere;
//  }


//_state.mc->create_shape(object, grid_real_size);


boost::array<int, 2> grid_shape = { 5*5*5 , 3 };
boost::multi_array<REAL, 2> grid(grid_shape);

boost::array<int, 1> implicit_values_shape = { 5*5*5 };
boost::multi_array<REAL, 1> implicit_values(implicit_values_shape);

boost::array<int, 2> gradou_values_shape = { 5*5*5,3 };
boost::multi_array<REAL, 2> gradou_values(gradou_values_shape);

for (int z = 0; z < 5; z++ ) {
    for (int y = 0; y < 5; y++ ) {
        for (int x = 0; x < 5; x++ ) {
            grid[x + y*5 + z*25][0] = 2.*(REAL)x/(REAL)5. -1.;
            grid[x + y*5 + z*25][1] =2.*(REAL)y/(REAL)5.-1.;
            grid[x + y*5 + z*25][2] = 2.*(REAL)z/(REAL)5. -1.;


        }
    }
}

object->eval_implicit(grid, &implicit_values);

object->eval_gradient(grid, &gradou_values);

ofstream f_out("/home/solene/Desktop/mp5-private/solidmodeler/js_iteration_2/tests/gradou_values.txt");

f_out << "Gradient_function in cpp :" << endl;
for (int i=0; i<125.; i++){
    f_out << gradou_values[i][0];
    f_out << " " ;
    f_out << gradou_values[i][1];
    f_out << " " ;
    f_out << gradou_values[i][2];
    f_out <<  "\n";
}
f_out << endl;

}
