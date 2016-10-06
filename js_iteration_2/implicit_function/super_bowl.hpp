#pragma once
#include "basic_data_structures.hpp"
#include "basic_functions.hpp"

namespace mp5_implicit {

class super_bowl : public transformable_implicit_function {

protected:
    REAL r;
    REAL x; REAL y; REAL z;

public:
    super_bowl(REAL radius){
        this->r = radius;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];
        for (int i=0; i<12; i++){
          if(i==0 || i==5 || i==10){
            this->transf_matrix[i] = 1;
            this->inv_transf_matrix[i] = 1;
          }
          else{
            this->transf_matrix[i] = 0;
            this->inv_transf_matrix[i] = 0;
          }
        }
    }
    super_bowl(REAL radius, REAL center_x, REAL center_y, REAL center_z){
        this->r = radius;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];
        for (int i=0; i<12; i++){
          if(i==0 || i==5 || i==10){
            this->transf_matrix[i] = 1;
            this->inv_transf_matrix[i] = 1;
          }
          else{
            this->transf_matrix[i] = 0;
            this->inv_transf_matrix[i] = 0;
          }
        }
    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL r2 = this->r;

        int output_ctr=0;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
            if(sqrt(pow((*i)[0]-this->x,2)+pow((*i)[1]-this->y,2)+pow((*i)[2]-this->z,2))<r2/1.4){
            (*f_output)[output_ctr] = -1;
            }
            else if(sqrt(pow((*i)[0]-this->x,2)+pow((*i)[1]-this->y,2)+pow((*i)[2]-this->z,2))<r2){
              (*f_output)[output_ctr] = -(pow((*i)[0]-this->x,2)*r2+pow((*i)[1]-this->y,2)*r2 - ((*i)[2]-this->z+r2/2));
            }
            else if((*i)[2]-this->z+r2/2>r2)
              (*f_output)[output_ctr]=-1;
            else
              (*f_output)[output_ctr]= -1;


        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL r2 = squared(this->r);
        int output_ctr=0;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
          if(sqrt(pow((*i)[0]-this->x,2)+pow((*i)[1]-this->y,2)+pow((*i)[2]-this->z,2))<r2/1.2){
            (*output)[output_ctr][0] = 2 * (*i)[0]-this->x;
            (*output)[output_ctr][1] = 2 * (*i)[1]-this->y;
            (*output)[output_ctr][2] = 2 * (*i)[2]-this->z;
          }
          else{
            (*output)[output_ctr][0] = -2 * (*i)[0]-this->x;
            (*output)[output_ctr][1] = -2 * (*i)[1]-this->y;
            (*output)[output_ctr][2] = 1;
          }
          if((*i)[2]-this->z+r2/2>r2){
            (*output)[output_ctr][0] = 0;
            (*output)[output_ctr][1] = 0;
            (*output)[output_ctr][2] = -1;
          }

          REAL g0 = (*output)[output_ctr][0];
          REAL g1 = (*output)[output_ctr][1];
          REAL g2 = (*output)[output_ctr][2];

          (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
          (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
          (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;
        }
    }
    bool integrity_invariant() const {
      if(this->r < MIN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = r*1.1;
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
