#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class double_mushroom : public transformable_implicit_function {

protected:
    REAL r;
    REAL a; REAL b; REAL c;// TODO convert to REAL_mm
    REAL x; REAL y; REAL z;

public:
    double_mushroom(REAL height_mm, REAL center_radius_x_mm, REAL center_radius_y_mm, REAL curvature){
        this->r = height_mm/2;
        this->a = center_radius_x_mm;
        this->b = center_radius_y_mm;
        this->c = 1/curvature;
        this->x = 0;
        this->y = 0;
        this->z = 0;
        //works with r>3********

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

    double_mushroom(REAL height_mm, REAL center_radius_x_mm, REAL center_radius_y_mm, REAL curvature, REAL center_x, REAL center_y, REAL center_z){
        this->r = height_mm/2;
        this->a = center_radius_x_mm;
        this->b = center_radius_y_mm;
        this->c = 1/curvature;
        this->x = center_x;
        this->y = center_y;
        this->z = center_z;
        //works with r>3********

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
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);

        int output_ctr=0;
        bool test = this->integrity_invariant();
        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){
            // REAL f = -(pow((*i)[0]-this->x,2)/a2+pow((*i)[1]-this->y,2)/b2-pow((*i)[2]-this->z,2)/c2-1);
            // REAL uperside = r-((*i)[2]-this->z);
            // REAL lowerside = r+((*i)[2]-this->z);
            //
            // (*f_output)[output_ctr] = min(f,min(uperside,lowerside));


            if ((*i)[2]>r){
              (*f_output)[output_ctr] = r - (*i)[2];
            }
            else if ((*i)[2]<-r){
              (*f_output)[output_ctr] = r + (*i)[2];
            }
            else{
              (*f_output)[output_ctr] = -(pow((*i)[0]-this->x,2)/a2+pow((*i)[1]-this->y,2)/b2-pow((*i)[2]-this->z,2)/c2-1);
            }

        }
    }
      virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL a2 = squared(this->a);
        const REAL b2 = squared(this->b);
        const REAL c2 = squared(this->c);

        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){

            REAL g0;
            REAL g1;
            REAL g2;

            REAL f = -(pow((*i)[0]-this->x,2)/a2+pow((*i)[1]-this->y,2)/b2-pow((*i)[2]-this->z,2)/c2-1);
            REAL uperside = r-((*i)[2]-this->z);
            REAL lowerside = r+((*i)[2]-this->z);

            if((*i)[2]<-r){
              g2 = 1;
              g1 = 0;
              g0 = 0;
            }
            else if ((*i)[2]>r){
              g2 = -1;
              g1 = 0;
              g0 = 0;
            }
            else{

              g0 = -2*((*i)[0]-this->x)/a2;
              g1 = -2*((*i)[1]-this->y)/b2;
              g2 =  2*((*i)[2]-this->z)/c2;
            }

            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*g0 + this->inv_transf_matrix[4]*g1 + this->inv_transf_matrix[8]*g2;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*g0 + this->inv_transf_matrix[5]*g1 + this->inv_transf_matrix[9]*g2;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*g0 + this->inv_transf_matrix[6]*g1 + this->inv_transf_matrix[10]*g2;
        }
    }
    bool integrity_invariant() const {
      if(MIN_PRINTABLE_LENGTH > this->r || MIN_PRINTABLE_LENGTH > this->a || MIN_PRINTABLE_LENGTH  > this->b)
        return false;
      else
        return true;
    }
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL max_size = norm_squared(r,a+c*r,b+c*r);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }
};

}
