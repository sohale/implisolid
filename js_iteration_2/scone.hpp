#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class scone : public implicit_function {

protected:
    REAL r;
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;

public:
    scone(REAL height, REAL radius_x, REAL radius_y, REAL radius_increase_speed ){
        this->r = height;
        this->a = radius_x;
        this->b = radius_y;
        this->c = 1/radius_increase_speed;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;
    }

    scone(REAL height, REAL radius_x, REAL radius_y, REAL radius_increase_speed, REAL center_x, REAL center_y, REAL center_z ){
        this->r = height;
        this->a = radius_x;
        this->b = radius_y;
        this->c = 1/radius_increase_speed;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
          if((*i)[2]-this->z >= 0.0){
            (*f_output)[output_ctr] = -((*i)[0]-this->x)*((*i)[0]-this->x)/a2 - ((*i)[1]-this->y)*((*i)[1]-this->y)/b2 + ((*i)[2]-this->z)*((*i)[2]-this->z)/c2;

          }
          if((*i)[2]-this->z >= this->r || (*i)[2]-this->z <= 0.0 ){
            (*f_output)[output_ctr] = -1.;
          }

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);
        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        const REAL r2 = squared(this->r);
        for(; i!=e; i++, output_ctr++){

          if((*i)[2]-this->z <= this->r){

            (*output)[output_ctr][0] = -2. * ((*i)[0]-this->x)/a2;
            (*output)[output_ctr][1] = -2. * ((*i)[1]-this->y)/b2;
            (*output)[output_ctr][2] = 2. * ((*i)[2]-this->z)/c2;

          }
          else {
            (*output)[output_ctr][0] = 0.;
            (*output)[output_ctr][1] = 0.;
            (*output)[output_ctr][2] = -1.;
          }

        }
    }
    bool integrity_invariant() const {
      if(this->r < MEAN_PRINTABLE_LENGTH || this->a < MEAN_PRINTABLE_LENGTH || this->b < MEAN_PRINTABLE_LENGTH || this->c < MEAN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(r,a+c*r,b+c*r);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
