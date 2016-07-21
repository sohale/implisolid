#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class scylinder : public transformable_implicit_function {

protected:
    REAL x; REAL y; REAL z;
    REAL* w; REAL c_len;
    REAL radius_u; REAL radius_v;
    REAL* transf_matrix;
    REAL* inv_transf_matrix;

public:
    scylinder(REAL matrix12[12]){
        this->radius_u = 0.5;
        this->radius_v = 0.5;
        this->c_len = 1.0;

        this->x = 0.;
        this->y = 0.;
        this->z = -0.5;

        this->w = new REAL[3];
        w[0] = 0;
        w[1] = 0.;
        w[2] = 1;
        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            transf_matrix[i] = matrix12[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");
    }

    scylinder(REAL W[3], REAL r_u, REAL r_v, REAL clen){
      this->radius_u = r_u;
      this->radius_v = r_v;
      this->c_len = clen;

      this->x = 0;
      this->y = 0;
      this->z = 0;

      this->w = new REAL[3];


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

    scylinder(REAL W[3], REAL r_u, REAL r_v, REAL clen, REAL center_x, REAL center_y, REAL center_z){
        this->radius_u = r_u;
        this->radius_v = r_v;
        this->c_len = clen;

        this->x = center_x;
        this->y = center_y;
        this->z = center_z;

        this->w = new REAL[3];

        for (int i=0; i<3; i++){
          this->w[i] = W[i];
        }

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

        int output_ctr=0;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
          REAL i0 = (*i)[0];
          REAL i1 = (*i)[1];
          REAL i2 = (*i)[2];
          REAL w0 = w[0];
          REAL w1 = w[1];
          REAL w2 = w[2];

          REAL t0 = (i0-this->x)*w0 + (i1-this->y)*w1 + (i2-this->z)*w2;
          REAL t1 = c_len - t0;
          REAL r_ = radius_u - sqrt((i0 - w0*t0 - this->x)*(i0 - w0*t0 - this->x)
            + (i1 - w1*t0- this->y)*(i1 - w1*t0 - this->y) +(i2 - w2*t0 - this->z)*(i2 - w2*t0 - this->z));

          (*f_output)[output_ctr] = min(t0,min(t1,r_));

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;

        Matrix_Vector_Product(this->inv_transf_matrix, x_copy);

        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();

        for(; i!=e; i++, output_ctr++){
          REAL i0 = (*i)[0];
          REAL i1 = (*i)[1];
          REAL i2 = (*i)[2];
          REAL w0 = w[0];
          REAL w1 = w[1];
          REAL w2 = w[2];

          REAL t0 = (i0-this->x)*w0 + (i1-this->y)*w1 + (i2-this->z)*w2;
          REAL t1 = c_len - t0;
          REAL r_ = radius_u - sqrt((i0 - w0*t0 - this->x)*(i0 - w0*t0 - this->x)
            + (i1 - w1*t0- this->y)*(i1 - w1*t0 - this->y) +(i2 - w2*t0 - this->z)*(i2 - w2*t0 - this->z));

          bool c_t0 = 0;
          bool c_t1 = 0;
          bool c_r = 0;

          if (t0 <= t1 && t0 <= r_){
            c_t0 = 1;
          }
          if (t1 <= t0 && t1 <= r_){
            c_t1 = 1;
          }
          if (r_ <= t0 && r_ <= t1){
            c_r = 1;
          }

          (*output)[output_ctr][0] = c_t0*w0 + c_t1*(-w0) + c_r*(w0*t0 + this->x - i0);
          (*output)[output_ctr][1] = c_t0*w1 + c_t1*(-w1) + c_r*(w1*t0 + this->y - i1);
          (*output)[output_ctr][2] = c_t0*w2 + c_t1*(-w2) + c_r*(w2*t0 + this->z - i2);

          REAL g0 = (*output)[output_ctr][0];
          REAL g1 = (*output)[output_ctr][1];
          REAL g2 = (*output)[output_ctr][2];

          (*output)[output_ctr][0] = this->transf_matrix[0]*g0 + this->transf_matrix[4]*g1 + this->transf_matrix[8]*g2;
          (*output)[output_ctr][1] = this->transf_matrix[1]*g0 + this->transf_matrix[5]*g1 + this->transf_matrix[9]*g2;
          (*output)[output_ctr][2] = this->transf_matrix[2]*g0 + this->transf_matrix[6]*g1 + this->transf_matrix[10]*g2;

        }
    }
    bool integrity_invariant() const {
      if(this->radius_u < MIN_PRINTABLE_LENGTH || this->c_len < MIN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(this->radius_u, this->c_len, 0.0);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
