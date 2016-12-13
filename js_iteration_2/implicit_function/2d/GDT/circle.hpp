#pragma once
#include "../../../configs.hpp"
#include "../basic_data_structures_2d.hpp"
//#include "../basic_functions_2d.hpp"
namespace mp5_implicit {

class circle : public implicit_function_2d {

protected:
    REAL center_x, center_y, radius;

public:
    circle(REAL center_x, REAL center_y, REAL radius) {
        this->center_x = center_x;
        this->center_y = center_y;
        this->radius = radius;

        my_assert(this->integrity_invariant(), "");
    }

    ~circle() {
    }


    virtual void eval_implicit(const vectorized_vect_2d& x, vectorized_scalar* f_output) const {

        my_assert(assert_implicit_function_io_2d(x, f_output), "");
        my_assert(this->integrity_invariant(), "");

        int output_ctr=0;
        for(auto i = x.begin(), e = x.end(); i < e; i++, output_ctr++){
            REAL x2 = std::pow((*i)[0] - this->center_x, 2);
            REAL y2 = std::pow((*i)[1] - this->center_y, 2);
            (*f_output)[output_ctr] = radius - std::sqrt(x2 + y2);
        }
    }

    virtual void eval_gradient(const vectorized_vect_2d& x, vectorized_vect_2d* output) const {

        int output_ctr=0;
        for(auto i = x.begin(), e = x.end(); i < e; i++, output_ctr++) {
            REAL x = (*i)[0] - this->center_x;
            REAL y = (*i)[1] - this->center_y;

            (*output)[output_ctr][0] = 2 * x;
            (*output)[output_ctr][1] = 2 * y;
        }
    }

    bool integrity_invariant() const {
      if(this->radius < CONFIG.MIN_PRINTABLE_LENGTH/2.0)
        return false;
      else
        return true;
    }

    virtual mp5_implicit::bounding_box_2d  get_boundingbox() const {
        return bounding_box_2d{center_x - radius, center_x + radius, center_y - radius, center_y + radius };
    }
};

}
