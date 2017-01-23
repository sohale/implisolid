
/**
 * File:
 * 		torus.hpp
 * Description:
 * 		Defines class torus, an implicit object with its implicit function
 *   	and implicit gradient.
 */

#pragma once
//#include "basic_data_structures.hpp"
//#include "basic_functions.hpp"

/* part of the namespace mp5_implicit */
namespace mp5_implicit {
namespace implicit_functions {


class torus : public transformable_implicit_function {

protected:
    REAL r, rx ,ry, rz;
    // what else should be on protected ?

public:

    torus(REAL r, REAL rx, REAL ry, REAL rz){
      this->r = r;
      this->rx = rx;
      this->ry = ry;
      this->rz = rz;
      this->transf_matrix = new REAL [12];
      this->inv_transf_matrix = new REAL [12];

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

    torus(REAL matrix12[12]){
        this->r = 1.;
        this->rx = 0.2;
        this->ry = 0.3;
        this->rz = 0.2;

        this->transf_matrix = new REAL [12];
        this->inv_transf_matrix = new REAL [12];

        for (int i=0; i<12; i++){
            this->transf_matrix[i] = matrix12[i];
        }

        invert_matrix(this->transf_matrix, this->inv_transf_matrix);
        my_assert(this->integrity_invariant(), "");

    }


    virtual void eval_implicit(const vectorized_vect & x, vectorized_scalar * f_output) const {
        my_assert(assert_implicit_function_io(x, * f_output), "");
        my_assert(this->integrity_invariant(), "");


        // declarations
        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix ,x);
        int output_ctr = 0;
        const REAL r = this->r;
        const REAL rx = this->rx;
        const REAL ry = this->ry;
        const REAL rz = this->rz;

        auto e = local_x.end();

        /* main loop for function */
        for (auto i = local_x.begin(); i != e; i++, output_ctr++){

            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL z = (*i)[2];

            /* the equation is: 1 - (r - ((x / rx)**2 + (y / ry)**2)**0.5)**2 - (z / rz)**2 */

            //REAL f = (1 - std::pow( r - std::pow((std::pow(x / rx, 2) + std::pow(y / ry,2)),0.5),2) - std::pow(z / rz, 2));
            //SDF
            (*f_output)[output_ctr] = 0.2 - euclidean_norm((euclidean_norm(x,y) -0.5), z);
        }
    }

    virtual void eval_gradient(const vectorized_vect & x, vectorized_vect * output) const {

        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix,x);
        int output_ctr = 0;
        const REAL r = this->r;
        const REAL rx = this->rx;
        const REAL ry = this->ry;
        const REAL rz = this->rz;

        auto e = local_x.end();

        /* main for-loop for gradient */

        for (auto i = local_x.begin(); i != e; i++, output_ctr++){

            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL z = (*i)[2];
            /*
            REAL a = std::pow(x*x/(rx*rx) + y*y/(ry*ry), 0.5);

            REAL gx = (2 * x /(rx*rx*a)) * (r - a);
            REAL gy = (2 * y / (ry*ry*a)) * (r - a);
            REAL gz = -2*z/(rz*rz);*/

            REAL A = euclidean_norm((euclidean_norm(x,y) -0.5), z);
            REAL B = euclidean_norm(x,y) -0.5;
            REAL C = euclidean_norm(x,y);
            REAL gx = - (x*B)/(C*A);
            REAL gy = -(y*B)/(C*A);
            REAL gz = - z/A;

            (*output)[output_ctr][0] = this->inv_transf_matrix[0]*gx + this->inv_transf_matrix[4]*gy + this->inv_transf_matrix[8]*gz;
            (*output)[output_ctr][1] = this->inv_transf_matrix[1]*gx + this->inv_transf_matrix[5]*gy + this->inv_transf_matrix[9]*gz;
            (*output)[output_ctr][2] = this->inv_transf_matrix[2]*gx + this->inv_transf_matrix[6]*gy + this->inv_transf_matrix[10]*gz;
        }

    }

    /* this acts as a consistency test, for example we know that r > 0
        Note: currently not implemented for torus
    */

    bool integrity_incariant() const {
        return true;
    }

    virtual mp5_implicit::bounding_box get_boundingbox() const{
        return mp5_implicit::bounding_box{-2, 2, -2, 2, -1,1};
    }

};

}  // namespace implicit_functions
}  // namespace mp5_implicit
