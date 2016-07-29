#pragma once
#include "basic_data_structures.hpp"
/**
 * File:
 * 		torus.hpp
 * Description:
 * 		Defines class torus, an implicit object with its implicit function
 *   	and implicit gradient.
 */

/* part of the namespace mp5_implicit */
namespace mp5_implicit {

class torus : public transformable_implicit_function {

protected:
    REAL rx ,ry, rz;
    // what else should be on protected ?

public:
    torus(REAL matrix12[12]){
        this->r = 4.
        this->rx = 0.1;
        this->ry = 0.1;
        this->rz = 0.1;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            this->transf_matrix[i] = matrix12[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");

    }


    virtual void eval_implicit(const vectorized_vect & x, vectorized_scalar * f_output) const {
        my_assert(assert_implicit_function_io(x, * f_output, "");
        my_assert(this->integrity_invariant(), "");


        // declarations
        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix ,x);


        // this is the main for loop that iterates over the input
        for (auto i = local_x.begin(); i != e; i++, output_ctr++){

            // the actual implicit function is computed here
        }
    }

    virtual void eval_gradient(const vectorized_vect & x, vectorized_scalar * f_output) const {

        // declarations


        // body of eval_gradient

        /* main for-loop for gradient */

        for (auto i = local_x.begin(; i != e; i++, output_ctr++){

            /* I will need the to compute the gradient by hand first to write this function */
        }

    }
}
}
