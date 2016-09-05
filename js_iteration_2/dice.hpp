#pragma once
#include "basic_data_structures.hpp"
#include "basic_functions.hpp"

namespace mp5_implicit {

class dice : public transformable_implicit_function {

protected:
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;


public:
    dice(REAL radius_x, REAL radius_y, REAL radius_z){
        this->a = radius_x;
        this->b = radius_y;
        this->c = radius_z;
        this->x = 0.;
        this->y = 0.;
        this->z = 0.;

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

    dice(REAL matrix[12]) {
        this->a = 6.;
        this->b = 2.5;
        this->c = 1.;

        this->x = 0.;
        this->y = 0.;
        this->z = 0.;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");
    }

    dice(REAL radius_x, REAL radius_y, REAL radius_z, REAL center_x, REAL center_y, REAL center_z){
        this->a = radius_x;
        this->b = radius_y;
        this->c = radius_z;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;

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

        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);
        const REAL r = this->a*this->a;
        int output_ctr=0;

        REAL cx = this->x;
        REAL cy = this->y;
        REAL cz = this->z;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
          REAL i1 = (*i)[0];
          REAL i2 = (*i)[1];
          REAL i3 = (*i)[2];

          // side of the cube
          REAL c1 = (i1 - cx - 0.5)*0.5*(-2.);
          REAL c2 = (i1 - cx + 0.5)*(-0.5)*(-2.);
          REAL c3 = (i2 - cy - 0.5)*0.5*(-2.);
          REAL c4 = (i2 - cy + 0.5)*(-0.5)*(-2.);
          REAL c5 = (i3 - cz - 0.5)*0.5*(-2.);
          REAL c6 = (i3 - cz + 0.5)*(-0.5)*(-2.);
          REAL cube = min(c1, min(c2, min(c3, min(c4, min(c5,c6)))));

          // face one
          REAL s1_1 = 0.01 - norm_squared(i1-0.5, i2-0., i3-0.);

          // face two
          REAL s2_1 = 0.01 - norm_squared(i1-0.2, i2-0.5, i3-0.2);
          REAL s2_2 = 0.01 - norm_squared(i1+0.2, i2-0.5, i3+0.2);


          //face four
          REAL s3_1 = 0.01 - norm_squared(i1+0.25, i2+0.25, i3+0.5);
          REAL s3_2 = 0.01 - norm_squared(i1-0., i2-0., i3+0.5);
          REAL s3_3 = 0.01 - norm_squared(i1-0.25, i2-0.25, i3+0.5);


          //face four
          REAL s4_1 = 0.01 - norm_squared(i1-0.2, i2+0.2, i3-0.5);
          REAL s4_2 = 0.01 - norm_squared(i1-0.2, i2-0.2, i3-0.5);
          REAL s4_3 = 0.01 - norm_squared(i1+0.2, i2+0.2, i3-0.5);
          REAL s4_4 = 0.01 - norm_squared(i1+0.2, i2-0.2, i3-0.5);


          //face five
          REAL s5_1 = 0.01 - norm_squared(i1-0.2, i2+0.5, i3-0.2);
          REAL s5_2 = 0.01 - norm_squared(i1+0.2, i2+0.5, i3-0.2);
          REAL s5_3 = 0.01 - norm_squared(i1+0.2, i2+0.5, i3+0.2);
          REAL s5_4 = 0.01 - norm_squared(i1-0.2, i2+0.5, i3+0.2);
          REAL s5_5 = 0.01 - norm_squared(i1+0., i2+0.5, i3-0.);


          //face six
          REAL s6_1 = 0.01 - norm_squared(i1+0.5, i2+0.25, i3-0.);
          REAL s6_2 = 0.01 - norm_squared(i1+0.5, i2+0.25, i3-0.25);
          REAL s6_3 = 0.01 - norm_squared(i1+0.5, i2+0.25, i3+0.25);
          REAL s6_4 = 0.01 - norm_squared(i1+0.5, i2-0.25, i3-0.);
          REAL s6_5 = 0.01 - norm_squared(i1+0.5, i2-0.25, i3-0.25);
          REAL s6_6 = 0.01 - norm_squared(i1+0.5, i2-0.25, i3+0.25);


          REAL spheres = max(s1_1,max(s2_1, max(s2_2,max(s6_1,max(s6_2, max(s6_3, max(s6_4, max(s6_5, max(s5_1, max(s5_2, max(s5_3, max(s5_4, max(s5_5,max(s4_1, max(s4_2, max(s4_3, max(s4_4, max(s3_1, max(s3_2, max(s3_3, s6_6))))))))))))))))))));

          if (cube < -spheres){
              (*f_output)[output_ctr] = cube;
          }
          else{
              (*f_output)[output_ctr] = -spheres;
          }
        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL r = this->a*this->a;

        REAL cx = this->x;
        REAL cy = this->y;
        REAL cz = this->z;

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

            // cube
            REAL c1 = (i1 - cx - 0.5)*0.5*(-2.);
            REAL c2 = (i1 - cx + 0.5)*(-0.5)*(-2.);
            REAL c3 = (i2 - cy - 0.5)*0.5*(-2.);
            REAL c4 = (i2 - cy + 0.5)*(-0.5)*(-2.);
            REAL c5 = (i3 - cz - 0.5)*0.5*(-2.);
            REAL c6 = (i3 - cz + 0.5)*(-0.5)*(-2.);

            REAL cube = min(c1, min(c2, min(c3, min(c4, min(c5,c6)))));

            // face one
            REAL s1_1 = 0.01 - norm_squared(i1-0.5, i2-0., i3-0.);

            // face two
            REAL s2_1 = 0.01 - norm_squared(i1-0.2, i2-0.5, i3-0.2);
            REAL s2_2 = 0.01 - norm_squared(i1+0.2, i2-0.5, i3+0.2);


            //face four
            REAL s3_1 = 0.01 - norm_squared(i1+0.25, i2+0.25, i3+0.5);
            REAL s3_2 = 0.01 - norm_squared(i1-0., i2-0., i3+0.5);
            REAL s3_3 = 0.01 - norm_squared(i1-0.25, i2-0.25, i3+0.5);


            //face four
            REAL s4_1 = 0.01 - norm_squared(i1-0.2, i2+0.2, i3-0.5);
            REAL s4_2 = 0.01 - norm_squared(i1-0.2, i2-0.2, i3-0.5);
            REAL s4_3 = 0.01 - norm_squared(i1+0.2, i2+0.2, i3-0.5);
            REAL s4_4 = 0.01 - norm_squared(i1+0.2, i2-0.2, i3-0.5);


            //face five
            REAL s5_1 = 0.01 - norm_squared(i1-0.2, i2+0.5, i3-0.2);
            REAL s5_2 = 0.01 - norm_squared(i1+0.2, i2+0.5, i3-0.2);
            REAL s5_3 = 0.01 - norm_squared(i1+0.2, i2+0.5, i3+0.2);
            REAL s5_4 = 0.01 - norm_squared(i1-0.2, i2+0.5, i3+0.2);
            REAL s5_5 = 0.01 - norm_squared(i1+0., i2+0.5, i3-0.);


            //face six
            REAL s6_1 = 0.01 - norm_squared(i1+0.5, i2+0.25, i3-0.);
            REAL s6_2 = 0.01 - norm_squared(i1+0.5, i2+0.25, i3-0.25);
            REAL s6_3 = 0.01 - norm_squared(i1+0.5, i2+0.25, i3+0.25);
            REAL s6_4 = 0.01 - norm_squared(i1+0.5, i2-0.25, i3-0.);
            REAL s6_5 = 0.01 - norm_squared(i1+0.5, i2-0.25, i3-0.25);
            REAL s6_6 = 0.01 - norm_squared(i1+0.5, i2-0.25, i3+0.25);


            REAL spheres = max(s1_1,max(s2_1, max(s2_2,max(s6_1,max(s6_2, max(s6_3, max(s6_4, max(s6_5, max(s5_1, max(s5_2, max(s5_3, max(s5_4, max(s5_5,max(s4_1, max(s4_2, max(s4_3, max(s4_4, max(s3_1, max(s3_2, max(s3_3, s6_6))))))))))))))))))));

            //substraction
            if (cube < - spheres){ // cube
                  int index = 0;
                  if (cube == c1){
                    g0 = -0.5;
                    g1 = 0.;
                    g2 = 0.;
                  }
                  else if(cube == c2){
                    g0 = +0.5;
                    g1 = 0.;
                    g2 = 0.;
                  }
                  else if (cube == c3){
                    g0 = 0.;
                    g1 = -0.5;
                    g2 = 0.;
                  }
                  else if (cube == c4){
                    g0 = 0.;
                    g1 = 0.5;
                    g2 = 0.;
                  }
                  else if (cube == c5){
                    g0 = 0.;
                    g1 = 0.;
                    g2 = -0.5;
                  }
                  else{
                    g0 = 0.;
                    g1 = 0.;
                    g2 = 0.5;
                  }

              }
            else{
              if (spheres == s1_1){
                  g0 = -2.*(i1 - 0.5);
                  g1 = -2.*(i2 - 0.);
                  g2 = -2.*(i3 - 0.);
              }
              else if(spheres == s2_1){
                  g0 = -2.*(i1 - 0.2);
                  g1 = -2.*(i2 - 0.5);
                  g2 = -2.*(i3 - 0.2);
              }
              else if(spheres == s2_2){
                  g0 = -2.*(i1 + 0.2);
                  g1 = -2.*(i2 - 0.5);
                  g2 = -2.*(i3 + 0.2);
              }
              else if(spheres == s3_1){
                  g0 = -2.*(i1 +0.25);
                  g1 = -2.*(i2 +0.25);
                  g2 = -2.*(i3 + 0.5);
              }
              else if(spheres == s3_2){
                  g0 = -2.*(i1 - 0.);
                  g1 = -2.*(i2 - 0.);
                  g2 = -2.*(i3 + 0.5);
              }
              else if(spheres == s3_3){
                  g0 = -2.*(i1 - 0.25);
                  g1 = -2.*(i2 - 0.25);
                  g2 = -2.*(i3 + 0.5);
              }
              else if(spheres == s4_1){
                  g0 = -2.*(i1 - 0.2);
                  g1 = -2.*(i2 + 0.2);
                  g2 = -2.*(i3 - 0.5);
              }
              else if(spheres == s4_2){
                  g0 = -2.*(i1 - 0.2);
                  g1 = -2.*(i2 - 0.2);
                  g2 = -2.*(i3 - 0.5);
              }
              else if(spheres == s4_3){
                  g0 = -2.*(i1 + 0.2);
                  g1 = -2.*(i2 + 0.2);
                  g2 = -2.*(i3 - 0.5);
              }
              else if(spheres == s4_4){
                  g0 = -2.*(i1 + 0.2);
                  g1 = -2.*(i2 - 0.2);
                  g2 = -2.*(i3 - 0.5);
              }
              else if(spheres == s5_1){
                  g0 = -2.*(i1 - 0.2);
                  g1 = -2.*(i2 + 0.5);
                  g2 = -2.*(i3 - 0.2);
              }
              else if(spheres == s5_2){
                  g0 = -2.*(i1 + 0.2);
                  g1 = -2.*(i2 + 0.5);
                  g2 = -2.*(i3 - 0.2);
              }
              else if(spheres == s5_3){
                  g0 = -2.*(i1 + 0.2);
                  g1 = -2.*(i2 + 0.5);
                  g2 = -2.*(i3 + 0.2);
              }
              else if(spheres == s5_4){
                  g0 = -2.*(i1 - 0.2);
                  g1 = -2.*(i2 + 0.5);
                  g2 = -2.*(i3 + 0.2);
              }
              else if(spheres == s5_5){
                  g0 = -2.*(i1 - 0.);
                  g1 = -2.*(i2 + 0.5);
                  g2 = -2.*(i3 - 0.);
              }
              else if(spheres == s6_1){
                  g0 = -2.*(i1 + 0.5);
                  g1 = -2.*(i2 + 0.25);
                  g2 = -2.*(i3 - 0.);
              }
              else if(spheres == s6_2){
                  g0 = -2.*(i1 + 0.5);
                  g1 = -2.*(i2 + 0.25);
                  g2 = -2.*(i3 - 0.25);
              }
              else if(spheres == s6_3){
                  g0 = -2.*(i1 + 0.5);
                  g1 = -2.*(i2 + 0.25);
                  g2 = -2.*(i3 + 0.25);
              }
              else if(spheres == s6_4){
                  g0 = -2.*(i1 + 0.5);
                  g1 = -2.*(i2 - 0.25);
                  g2 = -2.*(i3 - 0.);
              }
              else if(spheres == s6_5){
                  g0 = -2.*(i1 + 0.5);
                  g1 = -2.*(i2 - 0.25);
                  g2 = -2.*(i3 - 0.25);
              }
              else{
                  g0 = -2.*(i1 + 0.5);
                  g1 = -2.*(i2 - 0.25);
                  g2 = -2.*(i3 + 0.25);
              }
            }


            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;
        }
    }
    bool integrity_invariant() const {
      if(this->a < MIN_PRINTABLE_LENGTH || this->b < MIN_PRINTABLE_LENGTH || this->c < MIN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(a,b,c);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
