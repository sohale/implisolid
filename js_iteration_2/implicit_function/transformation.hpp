#pragma once

#include "basic_data_structures.hpp"
#include "basic_functions.hpp"


class transformable_implicit_function : public mp5_implicit::implicit_function {

public:

    REAL* transf_matrix;
    REAL* inv_transf_matrix;

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

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const = 0;
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const = 0;

protected:
    virtual bool integrity_invariant() const {return true;};

public:
    virtual ~transformable_implicit_function() {};

protected:

    /* Makes a copy and applied the matrix. To be called inside the eval_implicit() and eval_gradient() */
    /*
    vectorized_vect prepare_inner_vectors(const vectorized_vect& x) const {
        //my_assert(assert_implicit_function_io(x, *f_output), "");
        //my_assert(this->integrity_invariant(), ""); // fixme: has problems
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        return x_copy;
    }
    */

};
