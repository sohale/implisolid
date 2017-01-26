#pragma once
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
#include <exception>
namespace mp5_implicit {
namespace implicit_functions {


class half_plane_z_SDF : public transformable_implicit_function {

protected:
    REAL z0;
    int direction; //0 = z0 - z >=0, 1= z - z0 >= 0

public:

    half_plane_z_SDF(REAL matrix12[12], REAL z0,int direction)
    {

      this->z0 = z0;
      this->direction = direction;
      try{
        if(this->direction != 0 && this->direction != 1){
            throw "Wrong direction argument";
        }
      }catch(std::string& e){
        std::cout<< e <<std::endl;
      }
      this->transf_matrix = new REAL [12];
      this->inv_transf_matrix = new REAL [12];

      for (int i=0; i<12; i++){
          transf_matrix[i] = matrix12[i];
      }

      invert_matrix(this->transf_matrix, this->inv_transf_matrix);
      my_assert(this->integrity_invariant(), "");
    }


    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {

        // todo: apply matrix_vector_product()
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");

        vectorized_vect x_copy = x;

        matrix_vector_product(this->inv_transf_matrix, x_copy);

        int output_ctr=0;

        auto i = x_copy.begin();
        auto e = x_copy.end();
        if(this->direction == 1){
          for(; i<e; i++, output_ctr++){
              (*f_output)[output_ctr] = (*i)[0]-z0;
          }
        }else{
          for(; i<e; i++, output_ctr++){
              (*f_output)[output_ctr] = -((*i)[0]-z0);
          }
        }
    }


    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {

        // todo: apply matrix_vector_product()

        // todo: apply inv_transf_matrix

        vectorized_vect x_copy = x;
        matrix_vector_product(this->inv_transf_matrix, x_copy);

        int output_ctr=0;
        auto i = x_copy.begin();
        auto e = x_copy.end();
        if(this->direction == 1){
          for(; i<e; i++, output_ctr++){

              const REAL gx = 0.f;
              const REAL gy = 0.f;
              const REAL gz = 1.f;

              (*output)[output_ctr][0] = this->inv_transf_matrix[0]*gx + this->inv_transf_matrix[4]*gy + this->inv_transf_matrix[8]*gz;
              (*output)[output_ctr][1] = this->inv_transf_matrix[1]*gx + this->inv_transf_matrix[5]*gy + this->inv_transf_matrix[9]*gz;
              (*output)[output_ctr][2] = this->inv_transf_matrix[2]*gx + this->inv_transf_matrix[6]*gy + this->inv_transf_matrix[10]*gz;
          }
        }else{
          for(; i<e; i++, output_ctr++){

              const REAL gx = 0.f;
              const REAL gy = 0.f;
              const REAL gz = -1.f;

              (*output)[output_ctr][0] = this->inv_transf_matrix[0]*gx + this->inv_transf_matrix[4]*gy + this->inv_transf_matrix[8]*gz;
              (*output)[output_ctr][1] = this->inv_transf_matrix[1]*gx + this->inv_transf_matrix[5]*gy + this->inv_transf_matrix[9]*gz;
              (*output)[output_ctr][2] = this->inv_transf_matrix[2]*gx + this->inv_transf_matrix[6]*gy + this->inv_transf_matrix[10]*gz;
          }
        }
    }

    static void getHalfZParameters( REAL& z0,
                    const pt::ptree& shapeparams_dict){

        std::cout << "big thing: " << shapeparams_dict.get<float>("z0") << std::endl;
        z0 = shapeparams_dict.get<float>("z0");

    }

};

}  // namespace implicit_functions
}  // namespace mp5_implicit
