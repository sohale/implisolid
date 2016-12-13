#pragma once
// #include "../../../basic_data_structures.hpp"  // REAL
#include "../basic_data_structures_2d.hpp"
//#include "../basic_functions_2d.hpp"
namespace mp5_implicit {

// todo: For convex version, see HyperFun's review by Pasko & Adzhiev (incl. Implicit Curved Polygons - HyperFun)

class concave_polygon : public implicit_function_2d {
public:
    typedef struct { REAL x, y; } xy_t;
    std::vector<REAL> nx;
    std::vector<REAL> ny;
    std::vector<REAL> n0;

protected:
    std::vector<xy_t> corners;

public:
    /* Corners have to be arranged counter-clockwise, and form a convex polygon. */
    concave_polygon (std::vector<xy_t> corners)
        :
        nx(corners.size()),
        ny(corners.size()),
        n0(corners.size())
    {
        this->corners = corners;
        int nc = corners.size();

        std::vector<REAL> dx(nc);
        std::vector<REAL> dy(nc);

        std::vector<REAL> a_x(nc);  // y = a_y * x + b_y
        std::vector<REAL> b_x(nc);

        std::vector<REAL> a_y(nc);  // x = a_x * y + b_x
        std::vector<REAL> b_y(nc);

        std::vector<bool> whichx(nc);  // if |dy| > |dx|

        /*
        std::vector<REAL> nx(nc);
        std::vector<REAL> ny(nc);
        std::vector<REAL> n0(nc);
        */

        for (int i = 0; i < nc; ++i) {
            int next_i = (i < nc)? i+1 : 0;
            dx[i] = corners[next_i].x - corners[i].x;
            dy[i] = corners[next_i].y - corners[i].y;

            a_x[i] = dx[i] / dy[i];
            b_x[i] = corners[i].y;

            a_y[i] = dy[i] / dx[i];
            b_y[i] = corners[i].x;

            whichx[i] = std::abs(dy[i]) > std::abs(dx[i]);
            // if whichy => divide by y => use:   x = a_x * y + b_x

            REAL d = std::sqrt( dx[i]*dx[i] + dy[i] * dy[i] );
            REAL dinv = (d > 0.00000001)? 1.0 / d : 0.0;
            nx[i] = + dy[i] * dinv;
            ny[i] = - dx[i] * dinv;  // outward if counter-clockwise

            n0[i] = corners[i].x * nx[i] + corners[i].y * ny[i];
            // (nx,ny) points outwards =>
            // x * nx + y * ny - n0 < 0 ==> inside
        }

        my_assert(this->integrity_invariant(), "");
    }

    ~concave_polygon () {
    }

protected:
    inline auto calc_minv(REAL x, REAL y) const {
        REAL minv;
        int which_min = -1;
        assert(n0.size() > 0);
        for (int j = 0; j < n0.size(); ++j) {
            // (nx,ny) points outwards =>
            // x * nx + y * ny - n0 < 0 ==> inside
            REAL v = - ( x * nx[j] + y * ny[j] - n0[j] );
            if (v < minv || which_min < 0) {
                minv = v;
                which_min = j;
            }
        }
        assert(which_min >= 0);
        return std::make_pair(minv, which_min);
    }
public:
    virtual void eval_implicit(const vectorized_vect_2d& x, vectorized_scalar* f_output) const {

        my_assert(assert_implicit_function_io_2d(x, f_output), "");
        my_assert(this->integrity_invariant(), "");

        int output_ctr=0;
        for (auto i = x.begin(), e = x.end(); i < e; i++, output_ctr++) {
            REAL x = (*i)[0];
            REAL y = (*i)[1];
            auto min_val_which_pair = this->calc_minv(x, y);
            (*f_output)[output_ctr] = min_val_which_pair.first;
        }
    }

    virtual void eval_gradient(const vectorized_vect_2d& x, vectorized_vect_2d* output) const {
        // rotate
        int output_ctr=0;
        for (auto j = x.begin(), e = x.end(); j < e; j++, output_ctr++) {
            REAL x = (*j)[0];
            REAL y = (*j)[1];
            auto min_val_which_pair = this->calc_minv(x, y);
            (*output)[output_ctr][0] = nx[min_val_which_pair.second];
            (*output)[output_ctr][1] = ny[min_val_which_pair.second];
        }
    }

    bool integrity_invariant() const {
        if (n0.size() < 3) {
            return false;
        }
        for (int i = 0; i < n0.size(); ++i) {
            int next_i = (i < corners.size())? i+1 : 0;
            REAL dx = corners[next_i].x - corners[i].x;
            REAL dy = corners[next_i].y - corners[i].y;

            if(std::abs(std::sqrt(dx*dx + dy*dy)) < MIN_PRINTABLE_LENGTH/ 10.0)
              return false;
           //todo: check convexity
           //todo: check counter-clockwise
        }
        return true;
    }

    virtual mp5_implicit::bounding_box_2d  get_boundingbox() const {
        // todo: fixme
        return bounding_box_2d{-100, 100, -100, 100 };
    }
};

}
