#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class cube : public implicit_function {

protected:
    REAL h; REAL w; REAL d;
    REAL x; REAL y; REAL z;

public:
    cube(REAL height, REAL width, REAL depth){
        this->h = height;
        this->w = width;
        this->d = depth;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;
    }

    cube(REAL height, REAL width, REAL depth, REAL center_x, REAL center_y, REAL center_z){
        this->h = height;
        this->w = width;
        this->d = depth;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output)const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            if (((*i)[0]-this->x<=w/2.) && ((*i)[1]-this->y<=d/2.) && ((*i)[2]-this->z<=h/2.) && ((*i)[0]-this->x >= -w/2.) && ((*i)[1]-this->y >= -d/2.) && ((*i)[2]-this->z >= -h/2.)){
              (*f_output)[output_ctr] = 1.;
            }
            else{

              (*f_output)[output_ctr] = - 1.;
            }

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {


        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i!=e; i++, output_ctr++){
            if ((*i)[0]-this->x < -this->w+0.05) {
                  (*output)[output_ctr][0] = +1.* (*i)[0];
            }
            else if((*i)[0]-this->x > this->w-0.05){
                (*output)[output_ctr][0] = -1. * (*i)[0];
            }
            else if ((*i)[1]-this->y < -this->d-0.05){
                  (*output)[output_ctr][1] = +1.*(*i)[1];
            }
            else if((*i)[1]-this->y < -this->d+0.05){
                (*output)[output_ctr][1] = -1.*(*i)[1];
            }
            else if ((*i)[2]-this->z < -this->h-0.05){
                  (*output)[output_ctr][2] = +1.*(*i)[2];
            }
            else if((*i)[2]-this->z < -this->h+0.05){
                (*output)[output_ctr][2] = -1.*(*i)[2];
            }
            else
              (*output)[output_ctr][0]= 1; // arbitrary value to avoid null vectors may need to be changed

        }
    }
    bool integrity_invariant() const {
      if(this->w < MEAN_PRINTABLE_LENGTH || this->d < MEAN_PRINTABLE_LENGTH || this->h < MEAN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }
};

}
