#pragma once

#include "basic_data_structures.hpp"
#include "basic_functions.hpp"

#include "transformed.hpp"

namespace mp5_implicit {

class linearly_transformed : public transformed {

private:
    const implicit_function *child;

public:
    linearly_transformed (implicit_function* object, REAL matrix[12])
        : transformed(matrix), child(object)
    {
        //my_assert(
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        const vectorized_scalar::index nsize = x.shape()[0];

        const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix, x);
        auto sf = make_shape_1d(nsize);
        vectorized_scalar f1 = vectorized_scalar(sf);  // first function
        child->eval_implicit(local_x, f_output);
    };

    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
      const vectorized_scalar::index nsize = x.shape()[0];

      auto sf = make_shape_1d(nsize);

      /*auto shape = boost::array<vectorized_vect::index, 2>{nsize, 3};
      vectorized_vect grad1 = vectorized_vect(shape);
      */

      const vectorized_vect local_x = prepare_inner_vectors(this->inv_transf_matrix, x);

      child->eval_gradient(local_x, output);


      vectorized_scalar::index output_ctr = 0;

      for (auto i = local_x.begin(), e = local_x.end(); i < e; i++, output_ctr++){
          REAL gx = (*output)[output_ctr][0];
          REAL gy = (*output)[output_ctr][1];
          REAL gz = (*output)[output_ctr][2];
          (*output)[output_ctr][0] = this->inv_transf_matrix[0]*gx + this->inv_transf_matrix[4]*gy + this->inv_transf_matrix[8]*gz;
          (*output)[output_ctr][1] = this->inv_transf_matrix[1]*gx + this->inv_transf_matrix[5]*gy + this->inv_transf_matrix[9]*gz;
          (*output)[output_ctr][2] = this->inv_transf_matrix[2]*gx + this->inv_transf_matrix[6]*gy + this->inv_transf_matrix[10]*gz;
      }
    };


    bool integrity_invariant() const {
        REAL svd_s1 = 1;
        REAL svd_s2 = 1;
        REAL svd_s3 = 1;
        if(
          svd_s1 < MIN_PRINTABLE_LENGTH ||
          svd_s2 < MIN_PRINTABLE_LENGTH ||
          svd_s3 < MIN_PRINTABLE_LENGTH
        )
            return false;
        else
            return true;
    }

    /*
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        //REAL max_size = norm_squared(h,r2,r2);
        //-max_size, max_size, -max_size, max_size, -max_size, max_size};
        REAL incorrect_dummy = 10.0;
        return mp5_implicit::bounding_box{-incorrect_dummy,incorrect_dummy,-incorrect_dummy,incorrect_dummy,-incorrect_dummy,incorrect_dummy};
    }
    */
    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        REAL incorrect_dummy = 1;  // NaN cannot be used in JS
        return mp5_implicit::bounding_box{-incorrect_dummy,incorrect_dummy,-incorrect_dummy,incorrect_dummy,-incorrect_dummy,incorrect_dummy};
    }

};

}
