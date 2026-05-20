#pragma once
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
namespace mp5_implicit {
namespace implicit_functions {

class heart : public transformable_implicit_function {

protected:
    // REAL a; REAL b; REAL c;
    // REAL cx; REAL cy; REAL cz;

    // See docs/math/implicit_primitives.tex
    // Classical Taubin (1993) heart surface: (TF1, TF2) = ( 9/4, 9/80 ) and TF2b=TF2, TP3=3
    // We use (TF1, TF2) = ( 9/4, 9/200 ) to make the heart more pointy and less flat at the top.
    REAL TF1 = 9./4.;
    REAL TF2 = 9./200.;
    REAL TF2b = 9./100.;
    REAL TF1b = (27./2);
    REAL TF1c = 27./200.;
    size_t TP3 = 3;

    // some helpers for readablilty and a DSL feel
    inline static REAL p2(REAL x) { return x * x; }
    inline static REAL p3(REAL x) { return x * x * x; }
    // inline static REAL p4(REAL x) { return x * x * x * x; }

    inline static void populate_identity_matrix(REAL m[12]) {
        for (int i=0; i<12; i++){
          if(i==0 || i==5 || i==10){
            m[i] = 1;
          }
          else{
            m[i] = 0;
          }
        }
    }

    void init_the_transform_matrix_as_identity() {
      this->transf_matrix = new REAL [12];
      this->inv_transf_matrix = new REAL [12];
      populate_identity_matrix(this->transf_matrix);
      populate_identity_matrix(this->inv_transf_matrix);
    }

    void init_the_transform_matrix_from_given(REAL matrix[12]) {
        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
    }
public:
    // Was never used
    heart(int IGNORE_THIS_CONSTRUCTOR /*REAL radius_x, REAL radius_y, REAL radius_z*/){
        // this->a = radius_x;
        // this->b = radius_y;
        // this->c = radius_z;
        // this->cx = 0.;
        // this->cy = 0.;
        // this->cz = 0.;

        this->init_the_transform_matrix_as_identity();

        my_assert(this->integrity_invariant(), "");
    }

    heart(REAL matrix[12]) {
        // this->a = 6.;
        // this->b = 2.5;
        // this->c = 1.;

        // this->cx = 0.;
        // this->cy = 0.;
        // this->cz = 0.;

        this->init_the_transform_matrix_from_given(matrix);
        my_assert(this->integrity_invariant(), "");
    }

    // Was never used
    heart(/*REAL radius_x, REAL radius_y, REAL radius_z,*/ /*REAL center_x, REAL center_y, REAL center_z*/){
        // this->a = radius_x;
        // this->b = radius_y;
        // this->c = radius_z;
        // this->cx = center_x;
        // this->cy = center_y;
        // this->cz = center_z;

        init_the_transform_matrix_as_identity();
        my_assert(this->integrity_invariant(), "");
      }


    virtual void eval_implicit(const vectorized_vect& X, vectorized_scalar* f_output) const {


        my_assert(this->integrity_invariant(), "");
        vectorized_vect X_copy = X;

        matrix_vector_product(this->inv_transf_matrix, X_copy);
        // const REAL r = this->a*this->a;
        int output_ctr=0;

        // REAL cx = this->cx;
        // REAL cy = this->cy;
        // REAL cz = this->cz;

        auto i = X_copy.begin();
        auto e = X_copy.end();
        for(; i<e; i++, output_ctr++){
          REAL u = (*i)[0];
          REAL v = (*i)[1];
          REAL w = (*i)[2];

          (*f_output)[output_ctr] = -(
            std::pow(
              p2(u) + TF1 * p2(v) + p2(w) - 1., TP3
            )
            - p2(u) * p3(w)
            - TF2 * p2(v) * p3(w)
          );

        }
    }
    virtual void eval_gradient(const vectorized_vect& X, vectorized_vect* output) const {

        vectorized_vect X_copy = X;
        matrix_vector_product(this->inv_transf_matrix, X_copy);

        // const REAL r = this->a*this->a;

        // REAL cx = this->cx;
        // REAL cy = this->cy;
        // REAL cz = this->cz;

        int output_ctr=0;
        auto i = X_copy.begin();
        auto e = X_copy.end();
        for(; i<e; i++, output_ctr++){

            REAL u = (*i)[0];
            REAL v = (*i)[1];
            REAL w = (*i)[2];

            REAL A = pow(p2(u) + TF1 * p2(v) + p2(w) - 1, 2);
            REAL g0 = -6. * u * A + 2. * u * p3(w);
            REAL g1 = -TF1b * v * A + TF2b * v * p3(w);
            REAL g2 = -6. * w * A + 3. * p2(u) * p2(w) + TF1c * p2(v) * p2(w);

            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;
        }
    }
    bool integrity_invariant() const {
      /*
      if(this->a < MIN_PRINTABLE_LENGTH || this->b < MIN_PRINTABLE_LENGTH || this->c < MIN_PRINTABLE_LENGTH)
        return false;
      else
      */
        return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        // Just to replicate the old logic/function faithfully as part of refactoring (knowing the values are incorrect)
        REAL a = 6.;
        REAL b = 2.5;
        REAL c = 1.;
        REAL max_size = norm_squared(a,b,c);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}  // namespace implicit_functions
}  // namespace mp5_implicit
