#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class egg : public implicit_function {

protected:
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;

public:
    egg(REAL radius_x, REAL radius_y, REAL radius_z){
        this->a = radius_x;
        this->b = radius_y;
        this->c = radius_z;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;
    }

    egg(REAL radius_x, REAL radius_y, REAL radius_z, REAL center_x, REAL center_y, REAL center_z){
        this->a = radius_x;
        this->b = radius_y;
        this->c = radius_z;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");



        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = 1 - norm_squared(((*i)[0]-this->x)/this->a, ((*i)[1]-this->y)/this->b, ((*i)[2]-this->z)/this->c);

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        const REAL a2 = pow(this->a,2);
        const REAL b2 = pow(this->b,2);
        const REAL c2 = pow(this->c,2);
        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*output)[output_ctr][0] = -2. * ((*i)[0]-this->x)/a2;
            (*output)[output_ctr][1] = -2. * ((*i)[1]-this->y)/b2;
            (*output)[output_ctr][2] = -2. * ((*i)[2]-this->z)/c2;
        }
    }
    bool integrity_invariant() const {
      if(this->a < MEAN_PRINTABLE_LENGTH || this->b < MEAN_PRINTABLE_LENGTH || this->c < MEAN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }
};

}
