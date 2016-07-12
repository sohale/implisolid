#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class scylinder : public implicit_function {

protected:
    REAL r; REAL h;
    REAL x; REAL y; REAL z;

public:
    scylinder(REAL radius, REAL height){
        this->r = radius;
        this->h = height;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;
    }

    scylinder(REAL radius, REAL height, REAL center_x, REAL center_y, REAL center_z){
        this->r = radius;
        this->h = height;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;
    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        const REAL r2 = squared(this->r);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
          (*f_output)[output_ctr] = -((*i)[0]-this->x)*((*i)[0]-this->x) - ((*i)[1]-this->y)*((*i)[1]-this->y) + r2;
          if((*i)[2]-this->z >= this->h/2 || (*i)[2]-this->z <= -this->h/2 ){
            (*f_output)[output_ctr] = -1.;
          }

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        const REAL r2 = squared(this->r);
        for(; i!=e; i++, output_ctr++){

          if((*i)[2]-this->z >= this->h/2 || (*i)[2]-this->z <= -this->h/2){

            (*output)[output_ctr][0] = -2. * ((*i)[0]-this->x);
            (*output)[output_ctr][1] = -2. * ((*i)[1]-this->y);
            (*output)[output_ctr][2] = 0.;

          }
          else {
            (*output)[output_ctr][0] = 0.;
            (*output)[output_ctr][1] = 0.;
            (*output)[output_ctr][2] = -1.;
          }

        }
    }
    bool integrity_invariant() const {
      if(this->r < MEAN_PRINTABLE_LENGTH || this->h < MEAN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(r, h, 0.0);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
