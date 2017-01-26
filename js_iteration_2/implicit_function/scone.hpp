#pragma once
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"

namespace mp5_implicit {
namespace implicit_functions {

class scone : public transformable_implicit_function {

protected:
    REAL h;
    REAL r1; REAL r2;
    REAL x0; REAL y0; REAL z0;
    static const bool allow_zero_r1 = true;

public:

  scone(REAL matrix12[12]){
      this->h = 1.;
      this->r1 = 0.0;
      this->r2 = 0.5;
      this->x0 = 0.;
      this->y0 = 0.;
      this->z0 = +0.5;

      this->transf_matrix = new REAL [12];
      this->inv_transf_matrix = new REAL [12];

      for (int i=0; i<12; i++){
          transf_matrix[i] = matrix12[i];
      }

      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
      my_assert(this->integrity_invariant(), "");
  }
    scone(REAL height, REAL radius_min, REAL radius_max){
        this->h = height;
        this->r1 = radius_min;
        this->r2 = radius_max;
        this->x0 = 0.;
        this->y0 = 0.;
        this->z0 = 0.;

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

    scone(REAL height, REAL radius_min, REAL radius_max, REAL center_x, REAL center_y, REAL center_z ){
        this->h = height;
        this->r1 = radius_min;
        this->r2 = radius_max;
        this->x0 = center_x;
        this->y0 = center_y;
        this->z0 = center_z;

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


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        //vectorized_vect x_copy = x;
        //matrix_vector_product(this->inv_transf_matrix, x_copy);
        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix, x);

        const REAL a2 =this->r2/this->h;
        const REAL x0 = this->x0;
        const REAL y0 = this->y0;
        const REAL z0 = this->z0;

        int output_ctr=0;


        auto e = local_x.end();
        for(auto i = local_x.begin(); i<e; i++, output_ctr++){
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL z = (*i)[2];

            /*
            //REAL f = -(x-x0)*(x-x0) - (y-y0)*(y-y0) + (z-z0)*(z-z0)*a2;   // old, caused problems
            REAL f = - std::sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) 
                + std::sqrt((z-z0)*(z-z0)*a2);
            REAL uperside = -(z-z0)-r1;
            REAL lowerside = (z-z0)+h;

            (*f_output)[output_ctr] = min(f,min(uperside,lowerside));*/
            //REAL c1 = 1/sqrt(1+a2);
            //REAL c2 = sqrt(1-a2);
            //REAL e1 = std::min(euclidean_norm(,y)*c1 - z*c2, euclidean_norm(x,y)*c1 + z*c2);
            //REAL e2 = std::min(-(z-z0)+this->h, z-z0);
            //REAL e2 = std::max((z-z0)+this->h,-(z-z0));
            REAL c = 0.9;
            /*REAL a = -(c*euclidean_norm(x-x0,y-y0)+sqrt(1-c*c)*(z-z0));
            REAL b = z+0.5;*/
            (*f_output)[output_ctr] = -(c*euclidean_norm(x-x0,y-y0)+sqrt(1-c*c)*(z-z0));

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        //vectorized_vect local_x = x;
        //matrix_vector_product(this->inv_transf_matrix, local_x);
        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix, x);

        REAL a2 = squared(this->r2/this->h);

        const REAL x0 = this->x0;
        const REAL y0 = this->y0;
        const REAL z0 = this->z0;

        int output_ctr=0;

        auto e = local_x.end();
        //const REAL r2 = squared(this->h);
        for(auto i = local_x.begin(); i!=e; i++, output_ctr++){

            //shane
            /*
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL z = (*i)[2];

            REAL f = -(x-x0)*(x-x0)/a2 - (y-y0)*(y-y0)/a2 + (z-z0)*(z-z0);
            REAL uperside = -(z-z0)-r1;
            REAL lowerside = (z-z0)+h;

            REAL gx;
            REAL gy;
            REAL gz;

            if(uperside < f && uperside < lowerside){
                gx = 0.;
                gy = 0.;
                gz = -1.;
              }
            else if(lowerside < f && lowerside < uperside){
                gx = 0.;
                gy = 0.;
                gz = +1.;
              }
            else{
              gx = -2*(x-x0)/a2;
              gy = -2*(y-y0)/a2;
              gz = 2*(z-z0);
            }*/
            
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL z = (*i)[2];

            REAL gx;
            REAL gy;
            REAL gz;

            //REAL c1 = 1/sqrt(1+a2);
            //REAL c2 = sqrt(1-a2);
            REAL c = 0.9;

            gx = -((x-x0)*c)/euclidean_norm(x-x0,y-y0);
            gy = -((y-y0)*c)/euclidean_norm(x-x0,y-y0);
            gz = -sqrt(1-c*c);


            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*gx + this->inv_transf_matrix[4]*gy + this->inv_transf_matrix[8]*gz;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*gx + this->inv_transf_matrix[5]*gy + this->inv_transf_matrix[9]*gz;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*gx + this->inv_transf_matrix[6]*gy + this->inv_transf_matrix[10]*gz;

        }
    }
    bool integrity_invariant() const {
      if(
          this->h < MIN_PRINTABLE_LENGTH ||
          ((this->r1 < MIN_PRINTABLE_LENGTH) && (!scone::allow_zero_r1)) ||
          this->r2 < MIN_PRINTABLE_LENGTH
        )
          return false;
      else
          return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(h,r2,r2);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}  // namespace implicit_functions
}  // namespace mp5_implicit
