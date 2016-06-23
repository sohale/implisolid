
#include "basic_data_structures.hpp"

#include <iostream>


inline bool assert_implicit_function_io(const vectorized_vect& x, const vectorized_scalar* output){
    //std::cout << x.shape()[1] << " " << x.shape()[0] << " " << output->shape()[0] << std::endl;

    //my_assert(x.shape()[1] == 3, "Size should be N x 3. Not " ); //+ x.shape()[1]);
    my_assert(x.shape()[1] == 3, ""); //, "Size should be N x 3. Not " << x.shape()[1]);
    my_assert(x.shape()[0] == output->shape()[0], "") ; //, "Sizes don't match. Prepare an output using the same size. Sizes: " << x.shape()[0] << " !== " << output.shape()[0] );
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

class ImplicitFunction {
/*
    void func(const vectorized_vect& x, vectorized_scalar& output) = 0;
    void grad(const vectorized_vect& x, vectorized_vect& output) = 0;
*/
public:
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) = 0;
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) = 0;

protected:
    virtual bool integrity_invariant() = 0;

public:
    virtual ~ImplicitFunction() { };  // ?
};

//trait:
//SignedDistanceImplicitPointwise, PrimitiveBase


inline REAL squared(REAL x){
    return x*x;
}

inline REAL norm_squared(REAL x, REAL y, REAL z){
    return x*x + y*y + z*z;
}

class UnitSphere : public ImplicitFunction {

protected:
    REAL r;

public:
    UnitSphere(REAL r){
        this->r = r;
    }

    //boost::array<int, 2> big_shape = {{ 10000, 3 }};
    //boost::multi_array<REAL, 2> huge_test =  boost::multi_array<REAL, 2>(big_shape);

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output){
        my_assert(assert_implicit_function_io(x, f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);
        //auto i = x.begin();
        auto i_end = x.end();
        int output_ctr=0;
        //const vectorized_vect::iterator
        auto i = x.begin();
        for(; i<i_end; i++, output_ctr++){
            //f_output[output_ctr] = (*x)[0];
            (*f_output)[output_ctr] = r2 - norm_squared((*i)[0], (*i)[1], (*i)[2] );
        }
        /*
        for(auto i = x.begin(); i<i_end; i++, output_ctr++){
            //f_output[output_ctr] = (*x)[0];
            f_output[output_ctr] = r2 - norm_squared((*i)[0], (*i)[1], (*i)[2] );
        }
        */


      /*for (int i = 0; i<points.shape()[0];i++){
        f[i] = 1.0 - (pow(points[i][0],2) + pow(points[i][1],2) + pow(points[i][2],2));
      }*/

    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output){
        //(*output) = x;

        int output_ctr=0;
        auto i = x.begin();
        auto i_end = x.end();
        for(; i<i_end; i++, output_ctr++){
            (*output)[output_ctr][0] = -2. * (*i)[0];
            (*output)[output_ctr][1] = -2. * (*i)[1];
            (*output)[output_ctr][2] = -2. * (*i)[2];
        }
    }
    bool integrity_invariant(){
        return true;
    }
};

