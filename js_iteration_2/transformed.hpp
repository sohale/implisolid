#pragma once
#include "basic_data_structures.hpp"
#include "transformation.hpp"

namespace mp5_implicit {

class transformed : public transformable_implicit_function {

protected:
    REAL* transf_matrix;
    REAL* inv_transf_matrix;

public:
    //transformed(REAL matrix[12], const implicit_function& obj)
    transformed(REAL matrix[12]) {
        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");
    }

    virtual void rotate(const REAL angle, const vectorized_vect axis) const {
      REAL ca = cos(angle);
      REAL sa = sin(angle);
      REAL norm = sqrt(axis[0][0]*axis[0][0] + axis[0][1]*axis[0][1] + axis[0][2]*axis[0][2]);
      REAL a1 = axis[0][0]/norm;
      REAL a2 = axis[0][1]/norm;
      REAL a3 = axis[0][2]/norm;

      REAL rotation[12];
      rotation[0] = ca + a1*a1*(1.-ca);
      rotation[1] = a1*a2*(1.-ca) - a3*sa;
      rotation[2] = a1*a3*(1.-ca) + a2*sa;
      rotation[3] = 0.;
      rotation[4] = a1*a2*(1.-ca) + a3*sa;
      rotation[5] = ca + a2*a2*(1.-ca);
      rotation[6] = a2*a3*(1.-ca) - a1*sa;
      rotation[7] = 0.;
      rotation[8] = a1*a3*(1.-ca) - a2*sa;
      rotation[9] = a2*a3*(1.-ca) + a1*sa;
      rotation[10] = ca + a3*a3*(1.-ca);
      rotation[11] = 0.;

      matrix_matrix_product(this->transf_matrix, rotation);

      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
    }

    virtual void move(const vectorized_vect direction) const{
      this->transf_matrix[3] += direction[0][0];
      this->transf_matrix[7] += direction[0][1];
      this->transf_matrix[11] += direction[0][2];
      invert_matrix(this->transf_matrix, this->inv_transf_matrix);

    }
    virtual void resize(const REAL ratio) const{
      for (int i=0; i<12; i++){
        if(i==3 || i==7 || i==11){
        }
        else{
        this->transf_matrix[i] *= ratio;
        }
      }
      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const = 0;

    /* Makes a copy and applied the matrix. To be called inside the eval_implicit() and eval_gradient() */
    vectorized_vect prepare_inner_vectors(const vectorized_vect& x) const {
        //my_assert(assert_implicit_function_io(x, *f_output), "");
        //my_assert(this->integrity_invariant(), ""); // fixme: has problems
        vectorized_vect x_copy = x;

        Matrix_Vector_Product(this->inv_transf_matrix, x_copy);

        return x_copy;
    }
    void gradient_post_implace_transformation(const vectorized_vect& x) const {
        // todo(sohail):
        // Matrix_Vector_Product(this->inv_transf_matrix, x_copy);
    }

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const = 0;

    /*bool integrity_invariant() const {
        //TODO(sohail): Check the minimum size based on the SVD of the matrix.
        return true;
    }
    */

    //virtual mp5_implicit::bounding_box  get_boundingbox() const {
    //    REAL max_size = norm_squared(a,b,c);
    //    return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    //}
};

}
