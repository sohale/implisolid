/* Also See implicit_function.hpp
*/
#pragma once

#include "basic_data_structures_2d.hpp"
//using mp5_implicit::bounding_box_2d;
//using mp5_implicit::vectorized_vect_2d;

namespace mp5_implicit
{

class implicit_function_2d {

public:
    virtual void eval_implicit(const vectorized_vect_2d& x, vectorized_scalar* output) const = 0;
    virtual void eval_gradient(const vectorized_vect_2d& x, vectorized_vect_2d* output) const = 0;

protected:
    virtual bool integrity_invariant() const {return true;};

public:
    virtual ~implicit_function_2d() {};

    virtual bounding_box_2d  get_boundingbox() const = 0;
    // return bounding_box_2d{v1, v2, v3, v4, v5, v6};
};


}
