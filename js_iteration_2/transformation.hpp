#pragma once

#include "basic_data_structures.hpp"


class transformable_implicit_function : public implicit_function {

public:
    virtual void rotate(const REAL angle, const vectorized_vect axis) const = 0;
    virtual void move(const vectorized_vect direction) const = 0;
    virtual void resize(const REAL ratio) const = 0;
    virtual void eval_implicit(vectorized_vect& x, vectorized_scalar* output) const = 0;
    virtual void eval_gradient(vectorized_vect& x, vectorized_vect* output) const = 0;

protected:
    virtual bool integrity_invariant() const {return true;};

public:
    virtual ~transformable_implicit_function() {};
};
