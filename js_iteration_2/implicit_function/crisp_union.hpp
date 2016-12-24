#pragma once

#include "implicit_function.hpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"


/**
 * File: crisp_union.hpp
 * ----------------------

 * Defines the class CrispUnion which implements the csg union

 * operation between two implicit functions in accordance to the formula:

 * 	f(x) = max(f1(x),f2(x)). Functions f1 and f2 are implicit functions whose union

 * 	we want to compute.

 */

namespace mp5_implicit{
namespace implicit_functions {


class CrispUnion: public implicit_function {
public:

    /**
     * Function Declarations:
     * CrispUnion(a, b) --> The constructor of the class.
     * 		a and b are of type implicit_function, and need to be evaluated
     * 		on a grid or vector.
     *
     * eval_implicit --> evaluation of the implicit function.
     * eval_gradient --> evaluation of the implicit gradient.
     * ~CrispUnion --> Deconstructor of the class.
     */

    CrispUnion(const implicit_function & a_, const implicit_function & b_)
    : a(a_),b(b_)  // <-- member initialization list
                   // http://stackoverflow.com/questions/7665021/c-member-initialization-list
    {
        // of what type is a and b fa
    }

    void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        // checks that the size of input is the same as that of output
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const vectorized_scalar::index nsize = x.shape()[0];

        auto sf = make_shape_1d(nsize);

        vectorized_scalar f1 = vectorized_scalar(sf);  // first function
        vectorized_scalar f2 = vectorized_scalar(sf);  // second function

        a.eval_implicit(x, &f1);
        b.eval_implicit(x, &f2);

        vectorized_scalar::index output_ctr = 0;
        auto e = x.end();

        for (auto i = x.begin(); i < e; i++, output_ctr++){
            (*f_output)[output_ctr] = (f1[output_ctr] > f2[output_ctr]) ? (f1[output_ctr]): f2[output_ctr];
        }

    }

    void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
        const vectorized_scalar::index nsize = x.shape()[0];

        auto sf = make_shape_1d(nsize);

        vectorized_scalar f1 = vectorized_scalar(sf);  // first function
        vectorized_scalar f2 = vectorized_scalar(sf);  // second function

        auto shape = boost::array<vectorized_vect::index, 2>{{nsize,3}};
        vectorized_vect grad1 = vectorized_vect(shape);
        vectorized_vect grad2 = vectorized_vect(shape);

        a.eval_implicit(x, &f1);
        b.eval_implicit(x, &f2);

        a.eval_gradient(x, &grad1);
        b.eval_gradient(x, &grad2);

        vectorized_scalar::index output_ctr = 0;

        auto e = x.end();

        for (auto i = x.begin(); i < e; i++, output_ctr++){
            (*output)[output_ctr] = (f1[output_ctr] > f2[output_ctr]) ? (grad1[output_ctr]): grad2[output_ctr];
        }

    }

    // ~CrispUnion();   // will this compile with no warning?

private:
    const implicit_function &a, &b;  // reference a and b, this should be considered again.
};

}  // namespace implicit_functions
}  // namespace mp5_implicit
