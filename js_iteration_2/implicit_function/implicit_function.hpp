#pragma once

#include "../basic_data_structures.hpp"
//#include "../basic_functions.hpp"

namespace mp5_implicit {


/**
 * Class: implicit_function
 * ----------------------------------------
 * The base class which objects inherit from.
 *
 */

/*************************************************************
            Here is how each object is constructed: (please note that some classes name are written using all lowercase)

            NOTE : Each class takes as their three last arguments values for x for y and for z (REAL) to define the center of the figure.
                   Those are set to 0 if not defined. For now, you need to either define all three or none (constructor overloading)

    -CrispSubtract : takes two implicit_function as input and substract one from the other
    -cube : Takes a height, wigth and depth and computes the according Cube
    -double_mushroom : Takes a radius (horizontal size of the shape)  and a a, b and c value, dividing x², y² and z² in the formula for the hyperboloïde
     (note that those a,b,c must be choosen well for the shape to be closed (see the doc in double_mushroom.hpp))
    -egg : takes a, b, c dividing x², y² and z² in the formula for an ellipsoid. NOTE This should be changed to a 4 by 4 matrix filled with 9 inputs for total control over the shape
    -scone : takes a radius (horizontal size of the cone) and a,b,c dividing x²,y² and z² in the cone formula
    -scylinder : takes radius and height
    -super_bowl : takes radius NOTE : this shape is more of a first understanding of the crisp substract than a really usefull base shape
    -unit_sphere: take radius
 *************************************************************/
class implicit_function {

public:
    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* output) const = 0;
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const = 0;

protected:
    virtual bool integrity_invariant() const {return true;};

public:
    virtual ~implicit_function() {};

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        return mp5_implicit::bounding_box{-0.51, +0.51, -0.51, +0.51, -0.51, +0.51 };
    };

    // todo:
    // virtual vectorized_vect suggest_points_on_surface();
    // virtual vectorized_vect is_signed_distance_function();

};

}  // namespace mp5_implicit
