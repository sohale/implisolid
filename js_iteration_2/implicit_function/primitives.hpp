#pragma once

#include <iostream>

#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"

#include "../my_assert.hpp"

namespace implicit_functions {

inline bool assert_implicit_function_io(const vectorized_vect& x, const vectorized_scalar& output){
    //std::clog << x.shape()[1] << " " << x.shape()[0] << " " << output->shape()[0] << std::endl;

    //my_assert(x.shape()[1] == 3, "Size should be N x 3. Not " ); //+ x.shape()[1]);
    my_assert(x.shape()[1] == 3, ""); //, "Size should be N x 3. Not " << x.shape()[1]);
    my_assert(x.shape()[0] == output.shape()[0], "") ; //, "Sizes don't match. Prepare an output using the same size. Sizes: " << x.shape()[0] << " !== " << output.shape()[0] );
    return true;
}

}  // namespace implicit_functions

/*
s.implicit_func(x, f);

Sphere s;
s.eval_f(x, &f);

f = s.eval_f(x);


Sphere s = Sphere(2.);
int i = int(2);   // not (int)2

Sphere s(2.);

int i(2)   // construtor // int(2)

// Copy constructor
// Move Constructor

//ImplicitFunction s = Sphere();
Sphere s = Sphere();
s.implicit_func();
//s.implicit_func();


class TwistedObject::ImplicitObject(){
    ImplicitFunction obj;
    ImplicitFunction twist;

    //user operations
    rotate..
    ...
};
*/

//implicit_func, gradient

//implicit_func, gradient_func
//eval_implicit,


#include "implicit_function.hpp"
#include "transformation.hpp"

#include "linearly_transformed.hpp"

//#include "unit_sphere.hpp"
//
#include "crisp_subtract.hpp"
#include "crisp_union.hpp"
#include "crisp_intersection.hpp"
#include "unit_sphere.hpp"
#include "honey_comb.hpp"
#include "double_mushroom.hpp"
#include "egg.hpp"
#include "cube.hpp"
#include "super_bowl.hpp"
#include "scone.hpp"
#include "scylinder.hpp"
#include "legoland.hpp"
#include "dice.hpp"
#include "heart.hpp"
#include "transformed_union.hpp"
#include "transformed_intersection.hpp"
#include "transformed_subtract.hpp"
#include "torus.hpp"
#include "tetrahedron.hpp"
#include "screw.hpp"
#include "sdf_3d.hpp"
#include "meta_balls_Rydgard.hpp"
#include "javascript_implicit_function.hpp"
#include "extrusion.hpp"
#include "half_plane.hpp"
#include "top_bottom_lid.hpp"
#include "inf_top_bot_bound.hpp"
