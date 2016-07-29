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
    //
public:
    // Constructors


    virtual void eval_implicit(const vectorized_vect & x, vectorized_scalar * f_output  )
}
}
