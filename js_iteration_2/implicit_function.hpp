#pragma once

#include "basic_data_structures.hpp"
/**
 * Class: implicit_function
 * ----------------------------------------
 * The base class which objects inherit from.
 *
 */
class implicit_function {

public:
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output, REAL grid_real_size) const = 0;
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output, REAL grid_real_size) const = 0;

protected:
    virtual bool integrity_invariant() const {return true;};

public:
    virtual ~implicit_function() {};
};
