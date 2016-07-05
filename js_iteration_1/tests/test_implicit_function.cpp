// author ; Solene Chauvier & Marc Fraysse

#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "../mcc2_MS.cpp"

using namespace std;
using namespace mp5_implicit;

int main(){

REAL grid_real_size = 10;
// f_argument is made to always be between 0. and 1.
REAL f_argument = 0.5;

string name = "double_mushroom";
implicit_function * object;

if (name == "double_mushroom"){
double_mushroom mushroom(f_argument+3.); //3.3
object = &mushroom;
}
else if (name == "egg"){
egg segg(f_argument);
object = &segg; // super egg !
}
else if (name == "sphere"){
unit_sphere sphere(f_argument);
object = &sphere;
}
else if (name == "cube"){
cube cube(f_argument);
object = &cube;
}
else if (name == "super_bowl"){
super_bowl super_bowl(f_argument); //0.5
object = &super_bowl;
}
else if (name == "scone"){
scone scone(f_argument +2.5);
object = &scone;
}
else if (name == "scylinder"){
scylinder scylinder(f_argument); //0.7
object = &scylinder;
}
else {
cout << "Error! You must enter a valid name! So I made a sphere!" << endl;
unit_sphere sphere(f_argument);
object = &sphere;
}


boost::array<int, 2> grid_shape = { 5*5*5 , 3 };
boost::multi_array<REAL, 2> grid(grid_shape);

boost::array<int, 1> implicit_values_shape = { 5*5*5 };
boost::multi_array<REAL, 1> implicit_values(implicit_values_shape);

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

ofstream f_out("/home/solene/Desktop/mp5-private/solidmodeler/js_iteration_1/tests/comparison_implicit_func.txt");

f_out << "Implicit_function in cpp :" << endl;
for (int i=0; i<125.; i++){
    f_out << implicit_values[i];
    f_out <<  "\n";
}
f_out << endl;

}
