#pragma once
#include "basic_data_structures.hpp"

namespace mp5_implicit {
class unit_sphere : public implicit_function {

protected:
    REAL r;
    REAL x; REAL y; REAL z;

public:
    unit_sphere(REAL radius){
        //std::cout << "ctor sphere" << std::endl;
        this->r = radius;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;
    }

    unit_sphere(REAL radius, REAL center_x, REAL center_y, REAL center_z){
        //std::cout << "ctor sphere" << std::endl;
        this->r = radius;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;
    }


    virtual void eval_implicit(vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = r2 - norm_squared((*i)[0]-this->x, (*i)[1]-this->y, (*i)[2]-this->z );
        }
    }

    virtual void eval_gradient(vectorized_vect& x, vectorized_vect* output) const {

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*output)[output_ctr][0] = -2. * ((*i)[0]-this->x);
            (*output)[output_ctr][1] = -2. * ((*i)[1]-this->y);
            (*output)[output_ctr][2] = -2. * ((*i)[2]-this->z);
        }
    }
    bool integrity_invariant() const {
      if(this->r < MEAN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = r*1.1;
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}  // namespace
