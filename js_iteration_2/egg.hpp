#pragma once
#include "basic_data_structures.hpp"
namespace mp5_implicit {

class egg : public transformable_implicit_function {

protected:
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;
    REAL* transf_matrix;
    REAL* inv_transf_matrix;

public:
    egg(REAL radius_x, REAL radius_y, REAL radius_z){
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
    }

    egg(REAL radius_x, REAL radius_y, REAL radius_z, REAL center_x, REAL center_y, REAL center_z){
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

      }

    virtual void rotate(const REAL angle, const vectorized_vect axis) const {

    }
    virtual void move(const vectorized_vect direction) const{
      transf_matrix[3] += direction[0][0];
      transf_matrix[7] += direction[0][1];
      transf_matrix[11] += direction[0][2];
      InvertMatrix(transf_matrix, inv_transf_matrix);

    }
    virtual void resize(const REAL ratio) const{

    }
    virtual void eval_implicit(vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        Matrix_Vector_Product(this->inv_transf_matrix, x);

        int output_ctr=0;

        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*f_output)[output_ctr] = 1 - norm_squared(((*i)[0]-this->x)/this->a, ((*i)[1]-this->y)/this->b, ((*i)[2]-this->z)/this->c);
        }
    }
    virtual void eval_gradient(vectorized_vect& x, vectorized_vect* output) const {

        Matrix_Vector_Product(this->inv_transf_matrix, x);
        const REAL a2 = pow(this->a,2);
        const REAL b2 = pow(this->b,2);
        const REAL c2 = pow(this->c,2);
        int output_ctr=0;
        auto i = x.begin();
        auto e = x.end();
        for(; i<e; i++, output_ctr++){
            (*output)[output_ctr][0] = -2. * ((*i)[0]-this->x)/a2;
            (*output)[output_ctr][1] = -2. * ((*i)[1]-this->y)/b2;
            (*output)[output_ctr][2] = -2. * ((*i)[2]-this->z)/c2;
        }
    }
    bool integrity_invariant() const {
      if(this->a < MEAN_PRINTABLE_LENGTH || this->b < MEAN_PRINTABLE_LENGTH || this->c < MEAN_PRINTABLE_LENGTH)
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
