#pragma once
// #include "../../../basic_data_structures.hpp"  // REAL
#include "../basic_data_structures_2d.hpp"
//#include "../basic_functions_2d.hpp"
namespace mp5_implicit {

class upright_rectangle : public implicit_function_2d {

protected:
    REAL xmin, ymin, xmax, ymax;

public:
    upright_rectangle (REAL xmin, REAL ymin, REAL xmax, REAL ymax) {
        this->xmin = xmin;
        this->ymin = ymin;
        this->xmax = xmax;
        this->ymax = ymax;

        my_assert(this->integrity_invariant(), "");
    }

    ~upright_rectangle () {
    }


    virtual void eval_implicit(const vectorized_vect_2d& x, vectorized_scalar* f_output) const {

        my_assert(assert_implicit_function_io_2d(x, f_output), "");
        my_assert(this->integrity_invariant(), "");

        int output_ctr=0;
        for(auto i = x.begin(), e = x.end(); i < e; i++, output_ctr++) {
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL xm = std::min(x - this->xmin, this->xmax - x);
            REAL ym = std::min(y - ymin, ymax - y);
            (*f_output)[output_ctr] = std::min(xm, ym);
        }
    }

    virtual void eval_gradient(const vectorized_vect_2d& x, vectorized_vect_2d* output) const {
        // rotate
        int output_ctr=0;
        for(auto i = x.begin(), e = x.end(); i < e; i++, output_ctr++) {
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            REAL x0 = x - xmin;
            REAL x1 = xmax - x;
            REAL y0 = y - ymin;
            REAL y1 = ymax - y;
            REAL xm = std::min(x0, x1);
            REAL ym = std::min(y0, y1);
            REAL vx, vy;
            if (xm <= ym) {
                // x
                if (x0 <= x1) {
                    // use x0 == x - xmin
                    vx = +1.0; vy = 0.0;
                } else {
                    // use x1 == xmax - x
                    vx = -1.0; vy = 0.0;
                }
            } else {
                // y
                if (y0 <= y1) {
                    // use y0  == y - ymin
                    vx = 0.0; vy = +1.0;
                } else {
                    // use y1 == ymax - y
                    vx = 0.0; vy = -1.0;
                }
            }

            (*output)[output_ctr][0] = vx;
            (*output)[output_ctr][1] = vy;
        }
    }

    bool integrity_invariant() const {
      if(this->xmax - this->xmin < MIN_PRINTABLE_LENGTH)
        return false;
      if(this->ymax - this->ymin < MIN_PRINTABLE_LENGTH)
        return false;
      else
        return true;
    }

    virtual mp5_implicit::bounding_box_2d  get_boundingbox() const {
        return bounding_box_2d{xmin, xmax, ymin, ymax };
    }
};

}
