/**
 * File: crisp_subtract.cpp
 * ------------------------
 * Defines the class CrispSubtract which implements the csg subtraction
 * operation between two implicit functions.
 *
 * Implementation notes:
 *
 * As a first step we need to hardcode 2 implicit functions.
 */

class CrispSubtract: public implicit_function {
public:

    /**
     * Function Declarations:
     * CrispSubtract --> The constructor of the class.
     * eval_implicit --> evaluation of the implicit function of self.
     * eval_gradient --> evaluation of the implicit gradient of self.
     * ~CrispSubtract --> Deconstructor of the class.
     */
    // define the constructor

    CrispSubtract(implicit_function a, implicit_function b);
    void eval_implicit(const vectorized_vect& x, vectorized_scalar* output);
    void eval_gradient(const vectorized_vect& x, vectorized_scalar* output);
    ~CrispSubtract();
};

CrispSubtract::eval_implicit(const vectorized_vect& x, vectorized_scalar* output){};
