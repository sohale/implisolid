#pragma once
#include "basic_data_structures.hpp"
#include "basic_functions.hpp"

namespace mp5_implicit {

class legoland : public transformable_implicit_function {

protected:
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;


public:
    legoland(REAL radius_x, REAL radius_y, REAL radius_z){
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

    legoland(REAL matrix[12]) {
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

    legoland(REAL radius_x, REAL radius_y, REAL radius_z, REAL center_x, REAL center_y, REAL center_z){
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
          REAL c3 = (i2 - cy - 0.3)*0.3*(-2.);
          REAL c4 = (i2 - cy + 0.3)*(-0.3)*(-2.);
          REAL c5 = (i3 - cz - 0.2)*0.2*(-2.);
          REAL c6 = (i3 - cz + 0.2)*(-0.2)*(-2.);

          //first and second cylinder to be united
          REAL t0_1 = (i3);
          REAL t1_1 = 0.4 - t0_1;
          REAL r_1 = 0.15 - sqrt((i1  - 0.25)*(i1 - 0.25)
            + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

          REAL r_2 = 0.15 - sqrt((i1  + 0.25)*(i1 + 0.25)
            + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

          //third and fourth cylinder to be substracted
          REAL t0_3 = -(i3);
          REAL t1_3 =  0.5 - t0_1;
          REAL r_3 = 0.155 - sqrt((i1  - 0.25)*(i1 - 0.25)
            + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

          REAL r_4 = 0.155 - sqrt((i1  + 0.25)*(i1 + 0.25)
            + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

          // object implicit values
          REAL cube = min(c1, min(c2, min(c3, min(c4, min(c5,c6)))));
          REAL cyl_1 = min(t0_1, min(t1_1, r_1));
          REAL cyl_2 = min(t0_1, min(t1_1, r_2));
          REAL cyl_3 = min(min(t0_3, t1_3), r_3);
          REAL cyl_4 = min(min(t0_3, t1_3), r_4);

          // union and substraction
          if (max(cube, max(cyl_1, cyl_2)) < -max(cyl_3, cyl_4)){
            (*f_output)[output_ctr] = max(cube, max(cyl_1, cyl_2));
          }
          else{
              (*f_output)[output_ctr] = -max(cyl_3, cyl_4);
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
            REAL c3 = (i2 - cy - 0.3)*0.3*(-2.);
            REAL c4 = (i2 - cy + 0.3)*(-0.3)*(-2.);
            REAL c5 = (i3 - cz - 0.2)*0.2*(-2.);
            REAL c6 = (i3 - cz + 0.2)*(-0.2)*(-2.);

            // first and second cylinder to be united
            REAL t0_1 = (i3);
            REAL t1_1 = 0.4 - t0_1;
            REAL r_1 = 0.15 - sqrt((i1  - 0.25)*(i1 - 0.25)
              + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

            REAL r_2 = 0.15 - sqrt((i1  + 0.25)*(i1 + 0.25)
              + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

            // third and fourht cylinder to be substracted
            REAL t0_3 = -(i3);
            REAL t1_3 =  0.5 - t0_1;
            REAL r_3 = 0.155 - sqrt((i1  - 0.25)*(i1 - 0.25)
              + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

            REAL r_4 = 0.155 - sqrt((i1  + 0.25)*(i1 + 0.25)
              + (i2)*(i2) +(i3 - t0_1)*(i3 - t0_1));

            // object implicit values
            REAL cube = min(c1, min(c2, min(c3, min(c4, min(c5,c6)))));
            REAL cyl_1 = min(t0_1, min(t1_1, r_1));
            REAL cyl_2 = min(t0_1, min(t1_1, r_2));
            REAL cyl_3 = min(min(t0_3, t1_3), r_3);
            REAL cyl_4 = min(min(t0_3, t1_3), r_4);

            //substraction
            if (max(cube, max(cyl_1, cyl_2)) < -max(cyl_3, cyl_4)){
              //union
              if(cube > cyl_1 && cube > cyl_2){ // cube
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
                    g1 = -0.3;
                    g2 = 0.;
                  }
                  else if (cube == c4){
                    g0 = 0.;
                    g1 = 0.3;
                    g2 = 0.;
                  }
                  else if (cube == c5){
                    g0 = 0.;
                    g1 = 0.;
                    g2 = -0.2;
                  }
                  else{
                    g0 = 0.;
                    g1 = 0.;
                    g2 = 0.2;
                  }

              }
              else if (cyl_1 > cube && cyl_1 > cyl_2){ // first cylinder
                bool c_t0 = 0;
                bool c_t1 = 0;
                bool c_r = 0;

                if (t0_1 <= t1_1 && t0_1 <= r_1){
                  c_t0 = 1;
                }
                if (t1_1 <= t0_1 && t1_1 <= r_1){
                  c_t1 = 1;
                }
                if (r_1 <= t0_1 && r_1 <= t1_1){
                  c_r = 1;
                }

                g0 = + c_r*(0.25 - i1);
                g1 =  c_r*(- i2);
                g2 = c_t0*1. + c_t1*(-1.) + c_r*(1.*t0_1 - i3);

              }
              else{ // second cylinder
                bool c_t0 = 0;
                bool c_t1 = 0;
                bool c_r = 0;

                if (t0_1 <= t1_1 && t0_1 <= r_2){
                  c_t0 = 1;
                }
                if (t1_1 <= t0_1 && t1_1 <= r_2){
                  c_t1 = 1;
                }
                if (r_2 <= t0_1 && r_2 <= t1_1){
                  c_r = 1;
                }

                g0 = + c_r*(-0.25 - i1);
                g1 =  c_r*(- i2);
                g2 = c_t0*1. + c_t1*(-1.) + c_r*(1.*t0_1 - i3);

              }
            }
            else{ //substraction
              if (cyl_3 > cyl_4){ // - third cylinder
                bool c_t0 = 0;
                bool c_t1 = 0;
                bool c_r = 0;

                if (t0_3 <= t1_3 && t0_3 <= r_3){
                  c_t0 = 1;
                }
                if (t1_3 <= t0_3 && t1_3 <= r_3){
                  c_t1 = 1;
                }
                if (r_3 <= t0_3 && r_3 <= t1_3){
                  c_r = 1;
                }

                g0 = - c_r*(0.25 - i1);
                g1 =  -c_r*(- i2);
                g2 = -c_t0*1. - c_t1*(-1.) - c_r*(1.*t0_1 - i3);

              }
              else{ // -fourth cylinder
                bool c_t0 = 0;
                bool c_t1 = 0;
                bool c_r = 0;

                if (t0_3 <= t1_3 && t0_3 <= r_4){
                  c_t0 = 1;
                }
                if (t1_3 <= t0_3 && t1_3 <= r_4){
                  c_t1 = 1;
                }
                if (r_4 <= t0_3 && r_4 <= t1_3){
                  c_r = 1;
                }

                g0 = - c_r*(-0.25 - i1);
                g1 =  -c_r*(- i2);
                g2 = -c_t0*1. - c_t1*(-1.) - c_r*(1.*t0_1 - i3);

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
