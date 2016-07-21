#pragma once
#include <vector>
#include "basic_data_structures.hpp"
#include "transformed.hpp"

namespace mp5_implicit {

class transformed_union : public transformed {

    std::vector<implicit_function*> children;

public:
    transformed_union (std::vector<implicit_function*> children, REAL matrix[12])
        : transformed(matrix), children(children)
    {
        my_assert(children.size() == 2, "for now only works on two objects.");
    }

    virtual void eval_implicit(const vectorized_vect& x, vectorized_scalar* f_output) const {
        my_assert(assert_implicit_function_io(x, *f_output), "");
        my_assert(this->integrity_invariant(), "");
        const vectorized_scalar::index nsize = x.shape()[0];

        const vectorized_vect local_x = prepare_inner_vectors(x);
        auto sf = make_shape_1d(nsize);
        vectorized_scalar f1 = vectorized_scalar(sf);  // first function
        children[0]->eval_implicit(local_x, &f1);

        vectorized_scalar f2 = vectorized_scalar(sf);  // second function
        children[1]->eval_implicit(local_x, &f2);

        vectorized_scalar::index output_ctr = 0;
        auto e = x.end();
        for (auto i = x.begin(); i < e; i++, output_ctr++){
            (*f_output)[output_ctr] = (f1[output_ctr] > f2[output_ctr]) ? (f1[output_ctr]): f2[output_ctr];
        }
    };
    virtual void eval_gradient(const vectorized_vect& x, vectorized_vect* output) const {
      const vectorized_scalar::index nsize = x.shape()[0];

      auto sf = make_shape_1d(nsize);

      vectorized_scalar f1 = vectorized_scalar(sf);  // first function
      vectorized_scalar f2 = vectorized_scalar(sf);  // second function

      auto shape = boost::array<vectorized_vect::index, 2>{nsize, 3};
      vectorized_vect grad1 = vectorized_vect(shape);
      vectorized_vect grad2 = vectorized_vect(shape);

      children[0]->eval_implicit(x, &f1);
      children[1]->eval_implicit(x, &f2);

      children[0]->eval_gradient(x, &grad1);
      children[1]->eval_gradient(x, &grad2);

      vectorized_scalar::index output_ctr = 0;

      auto e = x.end();

      for (auto i = x.begin(); i < e; i++, output_ctr++){
          (*output)[output_ctr] = (f1[output_ctr] > f2[output_ctr]) ? (grad1[output_ctr]): grad2[output_ctr];
      }
    };
    bool integrity_invariant() const {

        //transformed(this).integrity_invariant();
        //???????????

        //TODO(sohail): Check the minimum size based on the SVD of the matrix.
        return true;
    }

    virtual mp5_implicit::bounding_box  get_boundingbox() const {
        for ( auto it = children.begin(); it < children.end(); it++ ){

        }
        REAL max_size = norm_squared(transf_matrix[0], transf_matrix[4], transf_matrix[8]);
        return mp5_implicit::bounding_box{-max_size, max_size, -max_size, max_size, -max_size, max_size};
    }

};

};
