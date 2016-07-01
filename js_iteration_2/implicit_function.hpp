#ifndef IMPLICIT_FUNCTION_HPP
#define IMPLICIT_FUNCTION_HPP

#include "basic_data_structures.hpp"

class implicit_function {

public:
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const = 0;
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const = 0;

protected:
     bool integrity_invariant(){return true;};

public:
    ~implicit_function() {};  // ?
};

//trait:
//SignedDistanceImplicitPointwise, PrimitiveBase

#endif // IMPLICIT_FUNCTION_HPP
