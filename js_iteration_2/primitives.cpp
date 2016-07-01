
#include "basic_data_structures.hpp"

#include <iostream>
#include "my_assert.hpp"

inline bool assert_implicit_function_io(const vectorized_vect& x, const vectorized_scalar& output){
    //std::cout << x.shape()[1] << " " << x.shape()[0] << " " << output->shape()[0] << std::endl;

    //my_assert(x.shape()[1] == 3, "Size should be N x 3. Not " ); //+ x.shape()[1]);
    my_assert(x.shape()[1] == 3, ""); //, "Size should be N x 3. Not " << x.shape()[1]);
    my_assert(x.shape()[0] == output.shape()[0], "") ; //, "Sizes don't match. Prepare an output using the same size. Sizes: " << x.shape()[0] << " !== " << output.shape()[0] );
    return true;
}

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

inline REAL squared(REAL x){
    return x*x;
}

inline REAL norm_squared(REAL x, REAL y, REAL z){
    return x*x + y*y + z*z;
}

#include "unit_sphere.cpp"
#include "double_mushroom.cpp"
#include "egg.cpp"
#include "cube.cpp"
#include "super_bowl.cpp"
#include "scone.cpp"
