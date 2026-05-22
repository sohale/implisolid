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
    REAL TF1 = 9./4.; // alpha
    REAL TF2 = 9./200.; // beta
    /*
    REAL TF2b = 9./100.;
    REAL TF1b = (27./2);
    REAL TF1c = 27./200.;
    */
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
    /*
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
    */

    heart(REAL matrix[12], REAL alpha, REAL beta, REAL P3) {
        // this->a = 6.;
        // this->b = 2.5;
        // this->c = 1.;

        // this->cx = 0.;
        // this->cy = 0.;
        // this->cz = 0.;

        this->TF1 = alpha;
        this->TF2 = beta;
        this->TP3 = P3;

        this->init_the_transform_matrix_from_given(matrix);
        my_assert(this->integrity_invariant(), "");
    }

    /*
    // Was never used
    heart(
      // REAL radius_x, REAL radius_y, REAL radius_z,*/ /*REAL center_x, REAL center_y, REAL center_z
      ){
        // this->a = radius_x;
        // this->b = radius_y;
        // this->c = radius_z;
        // this->cx = center_x;
        // this->cy = center_y;
        // this->cz = center_z;

        init_the_transform_matrix_as_identity();
        my_assert(this->integrity_invariant(), "");
      }
      */


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

          REAL A = p2(u) + TF1 * p2(v) + p2(w) - 1.;
          REAL f = (
            - std::pow(A, TP3)
            + ( p2(u)  + TF2 * p2(v) ) * p3(w)
          );

          (*f_output)[output_ctr] = f;

        }
    }
    // helper
    inline static void transforme_back_and_copy(   REAL & out0, REAL & out1, REAL & out2 , REAL invmat[12], REAL g0, REAL g1, REAL g2 ) {
        // copy_3(g0,g1,g2, this->inv_transf_matrix, (*output)[output_ctr][0], (*output)[output_ctr][1], (*output)[output_ctr][2], )
        // (*output)[output_ctr][0] = ...
        constexpr size_t // Math notation: row-then-column, starting with 1
          _11 = 0, _21 = 1, _31 = 2, _41 = 3,
          _12 = 4, _22 = 5, _32 = 6, _42 = 7,
          _13 = 8, _23 = 9, _33 = 10, _43 = 11;
          // This reveals an issue: The size (4x3 or 3x4?),
          // and also, the translation is not used!
          // todo: + invmat[_41]
        out0 = invmat[_11] * g0 + invmat[_12] * g1 + invmat[_13] * g2;
        out1 = invmat[_21] * g0 + invmat[_22] * g1 + invmat[_23] * g2;
        out2 = invmat[_31] * g0 + invmat[_32] * g1 + invmat[_33] * g2;
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

            /*
            REAL A = pow(p2(u) + TF1 * p2(v) + p2(w) - 1, 2);
            REAL g0 = -6. * u * A + 2. * u * p3(w);
            REAL g1 = -TF1b * v * A + TF2b * v * p3(w);
            REAL g2 = -6. * w * A + 3. * p2(u) * p2(w) + TF1c * p2(v) * p2(w);
            */

            // todo: swap nameing of a, b
            REAL a = u * u + TF2 * v * v;
            REAL b = u * u + TF1 * v * v + w * w;

            REAL g0, g1, g2;
            // REAL b2 = p2(b);
            // REAL b2 = std::pow(b,2);
            // REAL b2 = std::pow(b, TP3 - 1) *3. / 3.;
            REAL b2 = std::pow(b, TP3 - 1) * TP3 / 3.;
            g0 = 2 * p3(w) * u - b2 * u;
            g1 = 2 * p3(w)*TF2*v - b2 * TF1 * v;
            g2 = p2(w) * a - b2 * w;

            g0 = g0 * 3 / 2.;
            g1 = g1 * 3 / 2.;
            g2 = g2 * 3 / 2.;

            transforme_back_and_copy(
               (*output)[output_ctr][0],
               (*output)[output_ctr][1],
               (*output)[output_ctr][2],
               this->inv_transf_matrix,
               g0, g1, g2
            );

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
