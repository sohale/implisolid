#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class scone : public transformable_implicit_function {

protected:
    REAL h;
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;
    REAL* transf_matrix;
    REAL* inv_transf_matrix;

public:

  scone(REAL matrix12[12]){
      this->h = 1;
      this->a = 1;
      this->b = 1;
      this->c = 1;
      this->x = 0;
      this->y = 0;
      this->z = -0.5;

      this->transf_matrix = new REAL [12];
      this->inv_transf_matrix = new REAL [12];

      for (int i=0; i<12; i++){
          transf_matrix[i] = matrix12[i];
      }

      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
      my_assert(this->integrity_invariant(), "");
  }
    scone(REAL height, REAL radius_x, REAL radius_y, REAL radius_increase_speed ){
        this->h = height;
        this->a = radius_x;
        this->b = radius_y;
        this->c = 1/radius_increase_speed;
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
    }

    scone(REAL height, REAL radius_x, REAL radius_y, REAL radius_increase_speed, REAL center_x, REAL center_y, REAL center_z ){
        this->h = height;
        this->a = radius_x;
        this->b = radius_y;
        this->c = 1/radius_increase_speed;
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
    }

    virtual void rotate(const REAL angle, const vectorized_vect axis) const {
      REAL ca = cos(angle);
      REAL sa = sin(angle);
      REAL norm = sqrt(axis[0][0]*axis[0][0] + axis[0][1]*axis[0][1] + axis[0][2]*axis[0][2]);
      REAL a1 = axis[0][0]/norm;
      REAL a2 = axis[0][1]/norm;
      REAL a3 = axis[0][2]/norm;

      REAL rotation[12];
      rotation[0] = ca + a1*a1*(1.-ca);
      rotation[1] = a1*a2*(1.-ca) - a3*sa;
      rotation[2] = a1*a3*(1.-ca) + a2*sa;
      rotation[3] = 0.;
      rotation[4] = a1*a2*(1.-ca) + a3*sa;
      rotation[5] = ca + a2*a2*(1.-ca);
      rotation[6] = a2*a3*(1.-ca) - a1*sa;
      rotation[7] = 0.;
      rotation[8] = a1*a3*(1.-ca) - a2*sa;
      rotation[9] = a2*a3*(1.-ca) + a1*sa;
      rotation[10] = ca + a3*a3*(1.-ca);
      rotation[11] = 0.;

      matrix_matrix_product(this->transf_matrix, rotation);

      invert_matrix(this->transf_matrix, this->inv_transf_matrix);

    }

    virtual void move(const vectorized_vect direction) const{
      this->transf_matrix[3] += direction[0][0];
      this->transf_matrix[7] += direction[0][1];
      this->transf_matrix[11] += direction[0][2];
      invert_matrix(this->transf_matrix, this->inv_transf_matrix);

    }
    virtual void resize(const REAL ratio) const{
      for (int i=0; i<12; i++){
        if(i==3 || i==7 || i==11){
        }
        else{
        this->transf_matrix[i] *= ratio;
        }
      }
      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        vectorized_vect x_copy = x;

        Matrix_Vector_Product(this->inv_transf_matrix, x_copy);

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);

        const REAL x0 = this->x;
        const REAL y0 = this->y;
        const REAL z0 = this->z;

        int output_ctr=0;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL z = (*i)[2];
            REAL r = (1 - (z-z0)/this->c) / 2.;
            REAL f;
            if (((z-z0) >= 0.0)  && ((z-z0) < this->h) ) {  //  // Note that z is negated, because the cone will be reversed.
                f = - (x-x0)*(x-x0)/a2 - (y-y0)*(y-y0)/b2 + r*r;
            } else {
                f = -1.;
            }
            (*f_output)[output_ctr] = f;

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;

        Matrix_Vector_Product(this->inv_transf_matrix, x_copy);

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);
        const REAL x0 = this->x;
        const REAL y0 = this->y;
        const REAL z0 = this->z;

        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        const REAL r2 = squared(this->h);
        for(; i!=e; i++, output_ctr++){

            if((z-z0) <= this->h){  // wrong: need to check both 0 and h

                REAL x = (*i)[0];
                REAL y = (*i)[1];
                REAL z = (*i)[2];

                (*output)[output_ctr][0] = -2. * (x-x0)/a2;
                (*output)[output_ctr][1] = -2. * (y-y0)/b2;
                (*output)[output_ctr][2] = 2. * (z-z0)/c2;

            }
            else {

                (*output)[output_ctr][0] = 0.;
                (*output)[output_ctr][1] = 0.;
                (*output)[output_ctr][2] = -1.;
            }

            REAL g0 = (*output)[output_ctr][0];
            REAL g1 = (*output)[output_ctr][1];
            REAL g2 = (*output)[output_ctr][2];

            (*output)[output_ctr][0] = this->transf_matrix[0]*g0 + this->transf_matrix[4]*g1 + this->transf_matrix[8]*g2;
            (*output)[output_ctr][1] = this->transf_matrix[1]*g0 + this->transf_matrix[5]*g1 + this->transf_matrix[9]*g2;
            (*output)[output_ctr][2] = this->transf_matrix[2]*g0 + this->transf_matrix[6]*g1 + this->transf_matrix[10]*g2;

        }
    }
    bool integrity_invariant() const {
      if(this->h < MIN_PRINTABLE_LENGTH || this->a < MIN_PRINTABLE_LENGTH || this->b < MIN_PRINTABLE_LENGTH || this->c < MIN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(h,a+c*h,b+c*h);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
