#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class cube : public transformable_implicit_function {

protected:
    REAL p[18];
    REAL x; REAL y; REAL z;
    REAL* transf_matrix;
    REAL* inv_transf_matrix;

public:
    cube(REAL size_x, REAL size_y, REAL size_z){
        for (int i=0; i<18; i++){
          this->p[i] = 0;
        }
        this->p[0] = size_x;
        this->p[1] = size_x/2.;
        this->p[3] = -size_x;
        this->p[4] = -size_x/2.;
        this->p[8] = size_y/2.;
        this->p[7] = size_y;
        this->p[10] = -size_y;
        this->p[11] = -size_y/2.;
        this->p[14] = size_z;
        this->p[17] = -size_z;

        // REAL m = 50.;
        // this->p[0] = cos(1.*M_PI/5.)/m;
        // this->p[1] = sin(1.*M_PI/5)/m;
        // this->p[2] = 1./m;
        // this->p[3] = cos(3.*M_PI/5.)/m;
        // this->p[4] = sin(3*M_PI/5)/m;
        // this->p[5] = 1./m;
        // this->p[6] = cos(5.*M_PI/5.)/m;
        // this->p[7] = sin(5.*M_PI/5.)/m;
        // this->p[8] = 1./m;
        // this->p[9] = cos(7.*M_PI/5.)/m;
        // this->p[10] = sin(7.*M_PI/5.)/m;
        // this->p[11] = 1./m;
        // this->p[13] = cos(9.*M_PI/5.)/m;
        // this->p[12] = sin(9.*M_PI/5.)/m;
        // this->p[14] = 1./m;
        // this->p[17] = -1./4.;


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

    cube(REAL size_x, REAL size_y, REAL size_z, REAL center_x, REAL center_y, REAL center_z){
        for (int i=0; i<18; i++){
          this->p[i] = 0;
        }
        this->p[0] = size_x;
        this->p[3] = -size_x;
        this->p[7] = size_y;
        this->p[10] = -size_y;
        this->p[14] = size_z;
        this->p[17] = -size_z;

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

      Matrix_Matrix_Product(this->transf_matrix, rotation);

      InvertMatrix(this->transf_matrix, this->inv_transf_matrix);

    }

    virtual void move(const vectorized_vect direction) const{
      this->transf_matrix[3] += direction[0][0];
      this->transf_matrix[7] += direction[0][1];
      this->transf_matrix[11] += direction[0][2];
      InvertMatrix(this->transf_matrix, this->inv_transf_matrix);

    }
    virtual void resize(const REAL ratio) const{
      for (int i=0; i<12; i++){
        if(i==3 || i==7 || i==11){
        }
        else{
        this->transf_matrix[i] *= ratio;
        }
      }
      InvertMatrix(this->transf_matrix, this->inv_transf_matrix);
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output)const {

        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        vectorized_vect x_copy = x;

        Matrix_Vector_Product(this->inv_transf_matrix, x_copy);
        int output_ctr=0;

        REAL cx = this->x;
        REAL cy = this->y;
        REAL cz = this->z;

        REAL i1;
        REAL i2;
        REAL i3;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
          i1 = (*i)[0];
          i2 = (*i)[1];
          i3 = (*i)[2];
          (*f_output)[output_ctr] = min(
            (i1 - cx - p[0])*p[0]*(-2.) + (i2 - cy - p[1])*p[1]*(-2.) + (i3 - cz - p[2])*p[2]*(-2.),
            min((i1 - cx - p[0+1*3])*p[0+1*3]*(-2.) + (i2 - cy - p[1+1*3])*p[1+1*3]*(-2.) + (i3 - cz - p[2+1*3])*p[2+1*3]*(-2.),
            min((i1 - cx - p[0+2*3])*p[0+2*3]*(-2.) + (i2 - cy - p[1+2*3])*p[1+2*3]*(-2.) + (i3 - cz - p[2+2*3])*p[2+2*3]*(-2.),
            min((i1 - cx - p[0+3*3])*p[0+3*3]*(-2.) + (i2 - cy - p[1+3*3])*p[1+3*3]*(-2.) + (i3 - cz - p[2+3*3])*p[2+3*3]*(-2.),
            min((i1 - cx - p[0+4*3])*p[0+4*3]*(-2.) + (i2 - cy - p[1+4*3])*p[1+4*3]*(-2.) + (i3 - cz - p[2+4*3])*p[2+4*3]*(-2.),
            (i1 - cx - p[0+5*3])*p[0+5*3]*(-2.) + (i2 - cy - p[1+5*3])*p[1+5*3]*(-2.) + (i3 - cz - p[2+5*3])*p[2+5*3]*(-2.)
          )))));
        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {


      vectorized_vect x_copy = x;
      Matrix_Vector_Product(this->inv_transf_matrix, x_copy);

        int output_ctr=0;

      REAL cx = this->x;
      REAL cy = this->y;
      REAL cz = this->z;

      REAL i1;
      REAL i2;
      REAL i3;
      auto i = x_copy.begin();
      auto e = x_copy.end();
      for(; i<e; i++, output_ctr++){
        i1 = (*i)[0];
        i2 = (*i)[1];
        i3 = (*i)[2];
        int index = 0;
        REAL min = (i1 - cx - p[0])*p[0]*(-2.) + (i2 - cy - p[1])*p[1]*(-2.) + (i3 - cz - p[2])*p[2]*(-2.);
        for (int i=0; i<6; i++){
          if ((i1 - cx - p[0+i*3])*p[0+i*3]*(-2.) + (i2 - cy - p[1+i*3])*p[1+i*3]*(-2.) + (i3 - cz - p[2+i*3])*p[2+i*3]*(-2.) < min){
            index = i;
            min = (i1 - cx - p[0+i*3])*p[0+i*3]*(-2.) + (i2 - cy - p[1+i*3])*p[1+i*3]*(-2.) + (i3 - cz - p[2+i*3])*p[2+i*3]*(-2.);
          }
        }
        (*output)[output_ctr][0] = p[index*3+0];
        (*output)[output_ctr][1] = p[index*3+1];
        (*output)[output_ctr][2] = p[index*3+2];
      }

    }
    bool integrity_invariant() const {
      // if(this->p[0] < MEAN_PRINTABLE_LENGTH || this->p[7] < MEAN_PRINTABLE_LENGTH) //this->p[14] < MEAN_PRINTABLE_LENGTH)
      //   return false;
      // else
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(p[0],p[7],p[14]);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
