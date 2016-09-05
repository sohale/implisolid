#pragma once
#include "basic_data_structures.hpp"
#include "basic_functions.hpp"
namespace mp5_implicit {

class honey_comb : public transformable_implicit_function {

protected:
    REAL a; REAL b; REAL c;
    REAL x; REAL y; REAL z;


public:
    honey_comb(REAL radius_x, REAL radius_y, REAL radius_z){
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

    honey_comb(REAL matrix[12]) {
        this->a = 0.5;
        this->b = 30.;
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

    honey_comb(REAL radius_x, REAL radius_y, REAL radius_z, REAL center_x, REAL center_y, REAL center_z){
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

        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){

          //*******HERE ARE ALL THE SHAPES THAT WROK AND MAY PROVE USEFULL ONE DAY************//

                                        //Honey Comb//

            REAL f = sin(this->b*(*i)[0])+sin(this->b*(*i)[1])-sin(this->b*(*i)[2]);
            REAL bouding_sphere = -(*i)[0]*(*i)[0] - (*i)[1]*(*i)[1] -(*i)[2]*(*i)[2] +r;
            (*f_output)[output_ctr] = min(f, bouding_sphere);

                                        //Cube of spheres//

            // REAL f1 = c*(1-4*(pow((*i)[0],2)+pow((*i)[1],2)+pow((*i)[2],2))/(9*a*a*2.4)+ 17*(pow((*i)[0],2)+pow((*i)[1],2)+pow((*i)[2],2))/(9*a*a*2.4)-22*(pow((*i)[0],2)+pow((*i)[1],2)+pow((*i)[2],2))/(9*a*a*2.4));
            // REAL f2 = c*(1-4*(pow((*i)[0]-3,2)+pow((*i)[1]-3,2)+pow((*i)[2]-3,2))/(9*a*a)+ 17*(pow((*i)[0]-3,2)+pow((*i)[1]-3,2)+pow((*i)[2]-3,2))/(9*a*a)-22*(pow((*i)[0]-3,2)+pow((*i)[1]-3,2)+pow((*i)[2]-3,2))/(9*a*a));
            // REAL f3 = c*(1-4*(pow((*i)[0]+3,2)+pow((*i)[1]+3,2)+pow((*i)[2]+3,2))/(9*a*a)+ 17*(pow((*i)[0]+3,2)+pow((*i)[1]+3,2)+pow((*i)[2]+3,2))/(9*a*a)-22*(pow((*i)[0]+3,2)+pow((*i)[1]+3,2)+pow((*i)[2]+3,2))/(9*a*a));
            // REAL f4 = c*(1-4*(pow((*i)[0]-3,2)+pow((*i)[1]+3,2)+pow((*i)[2]+3,2))/(9*a*a)+ 17*(pow((*i)[0]-3,2)+pow((*i)[1]+3,2)+pow((*i)[2]+3,2))/(9*a*a)-22*(pow((*i)[0]-3,2)+pow((*i)[1]+3,2)+pow((*i)[2]+3,2))/(9*a*a));
            // REAL f5 = c*(1-4*(pow((*i)[0]+3,2)+pow((*i)[1]+3,2)+pow((*i)[2]-3,2))/(9*a*a)+ 17*(pow((*i)[0]+3,2)+pow((*i)[1]+3,2)+pow((*i)[2]-3,2))/(9*a*a)-22*(pow((*i)[0]+3,2)+pow((*i)[1]+3,2)+pow((*i)[2]-3,2))/(9*a*a));
            // REAL f6 = c*(1-4*(pow((*i)[0]+3,2)+pow((*i)[1]-3,2)+pow((*i)[2]+3,2))/(9*a*a)+ 17*(pow((*i)[0]+3,2)+pow((*i)[1]-3,2)+pow((*i)[2]+3,2))/(9*a*a)-22*(pow((*i)[0]+3,2)+pow((*i)[1]-3,2)+pow((*i)[2]+3,2))/(9*a*a));
            // REAL f7 = c*(1-4*(pow((*i)[0]-3,2)+pow((*i)[1]-3,2)+pow((*i)[2]+3,2))/(9*a*a)+ 17*(pow((*i)[0]-3,2)+pow((*i)[1]-3,2)+pow((*i)[2]+3,2))/(9*a*a)-22*(pow((*i)[0]-3,2)+pow((*i)[1]-3,2)+pow((*i)[2]+3,2))/(9*a*a));
            // REAL f8 = c*(1-4*(pow((*i)[0]+3,2)+pow((*i)[1]-3,2)+pow((*i)[2]-3,2))/(9*a*a)+ 17*(pow((*i)[0]+3,2)+pow((*i)[1]-3,2)+pow((*i)[2]-3,2))/(9*a*a)-22*(pow((*i)[0]+3,2)+pow((*i)[1]-3,2)+pow((*i)[2]-3,2))/(9*a*a));
            // REAL f9 = c*(1-4*(pow((*i)[0]-3,2)+pow((*i)[1]+3,2)+pow((*i)[2]-3,2))/(9*a*a)+ 17*(pow((*i)[0]-3,2)+pow((*i)[1]+3,2)+pow((*i)[2]-3,2))/(9*a*a)-22*(pow((*i)[0]-3,2)+pow((*i)[1]+3,2)+pow((*i)[2]-3,2))/(9*a*a));
            //
            // (*f_output)[output_ctr] = max(f1,max(f2,max(f3,max(f4,max(f5,max(f6,max(f7,max(f8,f9))))))));

                                        //MAETBALLL//

            //
            // REAL f1 = exp(-(pow((*i)[0]+3.5,2)+pow((*i)[1]+3.5,2)+pow((*i)[2]+3.5,2))/9)-0.05;
            // REAL f2 = exp(-(pow((*i)[0]-3.5,2)+pow((*i)[1]-3.5,2)+pow((*i)[2]-3.5,2))/9);
            // REAL f3 = exp(-(pow((*i)[0]-3.5,2)+pow((*i)[1]+3.5,2)+pow((*i)[2]-3.5,2))/9);
            // REAL f4 = exp(-(pow((*i)[0]+3.5,2)+pow((*i)[1]-3.5,2)+pow((*i)[2]+3.5,2))/9);
            // (*f_output)[output_ctr] = f1+f2+f3+f4;

                                        //?????//

        }
    }
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        const REAL r = this->a*this->a;

        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        for(; i<e; i++, output_ctr++){

                                  // Gradient for the honey comb //

          REAL f = sin(this->b*(*i)[0])+sin(this->b*(*i)[1])-sin(this->b*(*i)[2]);
          REAL bouding_sphere = -(*i)[0]*(*i)[0] - (*i)[1]*(*i)[1] -(*i)[2]*(*i)[2] +r;

          if(f<bouding_sphere){
            (*output)[output_ctr][0] = this->b*cos(this->b*(*i)[0]);
            (*output)[output_ctr][1] = this->b*cos(this->b*(*i)[1]);
            (*output)[output_ctr][2] = -this->b*cos(this->b*(*i)[2]);
          }
          else{
            (*output)[output_ctr][0] = -2*(*i)[0];
            (*output)[output_ctr][1] = -2*(*i)[0];
            (*output)[output_ctr][2] = -2*(*i)[0];
          }

            REAL g0 = (*output)[output_ctr][0];
            REAL g1 = (*output)[output_ctr][1];
            REAL g2 = (*output)[output_ctr][2];

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
