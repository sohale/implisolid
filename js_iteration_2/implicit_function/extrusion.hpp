#pragma once
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
namespace mp5_implicit {
namespace implicit_functions {

class extrusion : public transformable_implicit_function {

protected:
    //REAL h;
    REAL r;
    //REAL x; REAL y; REAL z;


public:
  /*
    extrusion(REAL height, REAL radius){
      this->h = height;
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
    }*/

    extrusion(REAL matrix12[12]){

      this->r = 0.5;

      this->transf_matrix = new REAL [12];
      this->inv_transf_matrix = new REAL [12];

      for (int i=0; i<12; i++){
          transf_matrix[i] = matrix12[i];
      }

      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
      my_assert(this->integrity_invariant(), "");
    }

    //copied from unit_sphere
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {

        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);
        const REAL r2 = squared(this->r);

        int output_ctr=0;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
            const REAL x = (*i)[0];
            const REAL y = (*i)[1];
            (*f_output)[output_ctr] = r2 - norm_squared(x, y, 0);
        }
    }


    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        //const REAL r = this->r;


        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
            const REAL px = (*i)[0];
            const REAL py = (*i)[1];
            const REAL pz = (*i)[2];

            const REAL gx = -2*px;
            const REAL gy = -2*py;
            const REAL gz = 0;

            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*gx + this->inv_transf_matrix[4]*gy + this->inv_transf_matrix[8]*gz;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*gx + this->inv_transf_matrix[5]*gy + this->inv_transf_matrix[9]*gz;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*gx + this->inv_transf_matrix[6]*gy + this->inv_transf_matrix[10]*gz;
        }
    }
/*
    bool integrity_invariant() const {
      return transformable_implicit_function::integrity_invariant()
    }
*/
/*
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(a,b,c);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }*/
};

}  // namespace implicit_functions
}  // namespace mp5_implicit
