#pragma once

#include "implicit_function.hpp"
#include "basic_data_structures.hpp"

/**
 * File: crisp_subtract.cpp
 * ------------------------
 * Defines the class CrispSubtract which implements the csg subtraction
 * operation between two implicit functions.
 *
 * Implementation notes:
 *
 * As a first step we need to hardcode 2 implicit functions.
 */

class CrispSubtract: public implicit_function {
public:

    /**
     * Function Declarations:
     * CrispSubtract(a, b) --> The constructor of the class.
     * 		a and b are of type implicit_function, and need to be evaluated
     * 		on a grid or vector.
     *
     * eval_implicit --> evaluation of the implicit function.
     * eval_gradient --> evaluation of the implicit gradient.
     * ~CrispSubtract --> Deconstructor of the class.
     */

    CrispSubtract(const implicit_function & a_, const implicit_function & b_)
    : a(a_),b(b_)
    {
        // of what type is a and b fa
    }
    void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output){
        int output_ctr = 0;
        auto i = x.begin();
        auto e = x.end();

        const int nsize = x.shape()[0];

        auto sf = make_shape_1d(nsize);

        vectorized_scalar f1 = vectorized_scalar(sf);  // first function
        vectorized_scalar f2 = vectorized_scalar(sf);  // second function

        a.eval_implicit(x, &f1);
        b.eval_implicit(x, &f2);

        for (; i < e; i++, output_ctr++){
            (*f_output)[output_ctr] = (f1[output_ctr] < -f2[output_ctr]) ? (f1[output_ctr]): -f2[output_ctr];
        }
    }
    void eval_gradient(const vectorized_vect& x, vectorized_scalar* output){

    }

    ~CrispSubtract();   // will this compile with no warning? 

private:
    const implicit_function &a, &b;  // reference a and b, this should be considered again.
};
