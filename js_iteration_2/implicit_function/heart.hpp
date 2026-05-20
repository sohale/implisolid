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
    size_t TP3 = 3;


public:
    heart(/*REAL radius_x, REAL radius_y, REAL radius_z*/){
        // this->a = radius_x;
        // this->b = radius_y;
        // this->c = radius_z;
        // this->cx = 0.;
        // this->cy = 0.;
        // this->cz = 0.;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];
        for (int i=0; i<12; i++){
          if(i==0 || i==5 || i==10){
            this->transf_matrix[i] = 1;
            this->inv_transf_matrix[i] = 1;
          }
          else{
            this->transf_matrix[i] = 0;
            this->inv_transf_matrix[i] = 0;
          }
        }
        my_assert(this->integrity_invariant(), "");
    }

    heart(REAL matrix[12]) {
        // this->a = 6.;
        // this->b = 2.5;
        // this->c = 1.;

        // this->cx = 0.;
        // this->cy = 0.;
        // this->cz = 0.;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");
    }

    heart(/*REAL radius_x, REAL radius_y, REAL radius_z,*/ /*REAL center_x, REAL center_y, REAL center_z*/){
        // this->a = radius_x;
        // this->b = radius_y;
        // this->c = radius_z;
        // this->cx = center_x;
        // this->cy = center_y;
        // this->cz = center_z;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];
        for (int i=0; i<12; i++){
          if(i==0 || i==5 || i==10){
            this->transf_matrix[i] = 1;
            this->inv_transf_matrix[i] = 1;
          }
          else{
            this->transf_matrix[i] = 0;
            this->inv_transf_matrix[i] = 0;
          }
        }
        my_assert(this->integrity_invariant(), "");
      }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {


        my_assert(this->integrity_invariant(), "");
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);
        // const REAL r = this->a*this->a;
        int output_ctr=0;

        // REAL cx = this->cx;
        // REAL cy = this->cy;
        // REAL cz = this->cz;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
          REAL i1 = (*i)[0];
          REAL i2 = (*i)[1];
          REAL i3 = (*i)[2];

          (*f_output)[output_ctr] = -(
            std::pow(
              i1*i1 + (this->TF1)*i2*i2 + i3*i3 - 1., this->TP3
            )
            - i1*i1*i3*i3*i3
            - (this->TF2) * i2*i2*i3*i3*i3
          );

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        // const REAL r = this->a*this->a;

        // REAL cx = this->cx;
        // REAL cy = this->cy;
        // REAL cz = this->cz;

        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){

            REAL g0;
            REAL g1;
            REAL g2;

            REAL i1 = (*i)[0];
            REAL i2 = (*i)[1];
            REAL i3 = (*i)[2];

            REAL a = pow(i1*i1 + (this->TF1)*i2*i2 + i3*i3 - 1, 2);
            g0 = -6.*i1*a + 2.*i1*i3*i3*i3;
            g1 = -(27./2)*i2*a + (this->TF2b)*i2*i3*i3*i3;
            g2 = -6.*i3*a + 3.*i1*i1*i3*i3 + (27./200.)*i2*i2*i3*i3;

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
        REAL max_size = norm_squared(a,b,c);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}  // namespace implicit_functions
}  // namespace mp5_implicit
