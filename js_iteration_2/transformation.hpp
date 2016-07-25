#pragma once

#include "basic_data_structures.hpp"


class transformable_implicit_function : public implicit_function {

public:
    virtual void rotate(const REAL angle, const vectorized_vect axis) const = 0;
    virtual void move(const vectorized_vect direction) const = 0;
    virtual void resize(const REAL ratio) const = 0;

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
