#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class egg : public transformable_implicit_function {

protected:
    REAL a, b, c;
    REAL x0, y0, z0;

public:
    egg(REAL radius_x, REAL radius_y, REAL radius_z){
        this->a = radius_x;
        this->b = radius_y;
        this->c = radius_z;
        this->x0 = 0.;
        this->y0 = 0.;
        this->z0 = 0.;

        // Where is the delete? Added to  --> ~destructor
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
        my_assert(this->integrity_invariant(), "");
    }

    egg(REAL matrix[12]) {
        this->a = 0.5;
        this->b = 0.5;
        this->c = 0.5;

        this->x0 = 0.;
        this->y0 = 0.;
        this->z0 = 0.;

        // How to make sure 12 elements are provided? (how to assert)
        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");
    }

   //Do we test if Inverse works in "bad condision" cases ?
   //Invariance should also check the matrix's condition (i.e. being non-singular).
   // Also the eignenvalues whouls be > min_size.
    egg(REAL radius_x, REAL radius_y, REAL radius_z, REAL center_x, REAL center_y, REAL center_z){
        this->a = radius_x;
        this->b = radius_y;
        this->c = radius_z;
        this->x0 = center_x;
        this->y0 = center_y;
        this->z0 = center_z;

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
        my_assert(this->integrity_invariant(), "");
      }

    ~egg() {
        delete this->transf_matrix;
        this->transf_matrix = 0;
        delete this->inv_transf_matrix;
        this->inv_transf_matrix = 0;
    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {

        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        int output_ctr=0;

        auto e = x_copy.end();
        for(auto i = x_copy.begin(); i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = 1 - norm_squared(((*i)[0]-this->x0)/this->a, ((*i)[1]-this->y0)/this->b, ((*i)[2]-this->z0)/this->c);
        }
    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);
        int output_ctr=0;
        auto e = x_copy.end();
        for(auto i = x_copy.begin(); i < e; i++, output_ctr++) {
            REAL g0 = -2. * ((*i)[0]-this->x0)/a2;
            REAL g1 = -2. * ((*i)[1]-this->y0)/b2;
            REAL g2 = -2. * ((*i)[2]-this->z0)/c2;

            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;
        }
    }

    bool integrity_invariant() const {
      if(this->a < MIN_PRINTABLE_LENGTH || this->b < MIN_PRINTABLE_LENGTH || this->c < MIN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(a,b,c);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
