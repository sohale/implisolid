#pragma once
#include "basic_data_structures.hpp"
#include "basic_functions.hpp"

#include "transformation.hpp"

namespace mp5_implicit {

class transformed : public transformable_implicit_function {

protected:

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


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const = 0;

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
    void gradient_post_implace_transformation(const vectorized_vect& x) const {
        // todo(sohail):
        // matrix_vector_product(this->inv_transf_matrix, x_copy);
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
