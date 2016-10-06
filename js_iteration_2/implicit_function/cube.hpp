#pragma once
//#include "../basic_data_structures.hpp"
//#include "../basic_functions.hpp"

namespace mp5_implicit {

class cube : public transformable_implicit_function {

protected:
    REAL p[18];
    REAL x; REAL y; REAL z;

public:
    cube(REAL matrix12[12]){
       for (int i=0; i<18; i++){
          this->p[i] = 0;
        }
        this->p[0] = 0.5;
        this->p[3] = -0.5;
        this->p[7] = 0.5;
        this->p[10] = -0.5;
        this->p[14] = 0.5;
        this->p[17] = -0.5;

        this->x = 0;
        this->y = 0;
        this->z = 0;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];
        for (int i=0; i<12; i++){

            this->transf_matrix[i] = matrix12[i];


        }
          invert_matrix(this->transf_matrix, this->inv_transf_matrix);
    }
    cube(REAL size_x, REAL size_y, REAL size_z){
        for (int i=0; i<18; i++){
          this->p[i] = 0;
        }
        this->p[0] = 0.5;
        this->p[3] = -0.5;
        this->p[7] = 0.5;
        this->p[10] = -0.5;
        this->p[14] = 0.5;
        this->p[17] = -0.5;

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


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output)const {

        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);
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
            // (X-C) . P
                (i1 - cx - p[0+0*3])*p[0+0*3]*(-2.) + (i2 - cy - p[1])*p[1]*(-2.) + (i3 - cz - p[2])*p[2]*(-2.),
            min((i1 - cx - p[0+1*3])*p[0+1*3]*(-2.) + (i2 - cy - p[1+1*3])*p[1+1*3]*(-2.) + (i3 - cz - p[2+1*3])*p[2+1*3]*(-2.),
            min((i1 - cx - p[0+2*3])*p[0+2*3]*(-2.) + (i2 - cy - p[1+2*3])*p[1+2*3]*(-2.) + (i3 - cz - p[2+2*3])*p[2+2*3]*(-2.),
            min((i1 - cx - p[0+3*3])*p[0+3*3]*(-2.) + (i2 - cy - p[1+3*3])*p[1+3*3]*(-2.) + (i3 - cz - p[2+3*3])*p[2+3*3]*(-2.),
            min((i1 - cx - p[0+4*3])*p[0+4*3]*(-2.) + (i2 - cy - p[1+4*3])*p[1+4*3]*(-2.) + (i3 - cz - p[2+4*3])*p[2+4*3]*(-2.),
            (i1 - cx - p[0+5*3])*p[0+5*3]*(-2.) + (i2 - cy - p[1+5*3])*p[1+5*3]*(-2.) + (i3 - cz - p[2+5*3])*p[2+5*3]*(-2.)
          )))));
        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {


      vectorized_vect x_copy = vectorized_vect(x);
      matrix_vector_product(this->inv_transf_matrix, x_copy);

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

        (*output)[output_ctr][0] = -p[index*3+0];
        (*output)[output_ctr][1] = -p[index*3+1];
        (*output)[output_ctr][2] = -p[index*3+2];

        REAL g0 = (*output)[output_ctr][0];
        REAL g1 = (*output)[output_ctr][1];
        REAL g2 = (*output)[output_ctr][2];

        (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
        (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
        (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;

      }

    }
    bool integrity_invariant() const {
      // if(this->p[0] < MIN_PRINTABLE_LENGTH || this->p[7] < MIN_PRINTABLE_LENGTH) || this->p[14] < MIN_PRINTABLE_LENGTH)
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
