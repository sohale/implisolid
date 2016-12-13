#pragma once
// #include "../../../basic_data_structures.hpp"  // REAL
#include "../basic_data_structures_2d.hpp"
#include "../../../configs.hpp"
//#include "../basic_functions_2d.hpp"

#include "../../../vectorised_algorithms/cross_product.hpp"
using mp5_implicit::vectorised_algorithms::cross_product;

namespace mp5_implicit {

// todo: For convex version, see HyperFun's review by Pasko & Adzhiev (incl. Implicit Curved Polygons - HyperFun)

class concave_polygon : public implicit_function_2d {
public:
    typedef struct { REAL x, y; } xy_t;
    std::vector<REAL> nx;
    std::vector<REAL> ny;
    std::vector<REAL> n0;
    bool invalidate = true;  // i.e. not consolidated

protected:
    std::vector<xy_t> corners;

public:
    static std::vector<xy_t> conv1(std::vector<REAL> corners_x, std::vector<REAL> corners_y) {
        int nc = corners_x.size();
        std::vector<xy_t> corners(nc);
        assert(corners_x.size() == corners_y.size());
        //this->corners = corners;
        for (int i = 0; i < nc; ++i) {
            corners[i].x = corners_x[i];
            corners[i].y = corners_y[i];
        }
        return corners;
    }
    concave_polygon (std::vector<REAL> corners_x, std::vector<REAL> corners_y)
        :
        //nx(corners_x.size()),
        //ny(corners_x.size()),
        //n0(corners_x.size()),
        //corners(conv1(corners_x, corners_y)),
        concave_polygon(conv1(corners_x, corners_y))
    {
    }

    /* Corners have to be arranged counter-clockwise, and form a convex polygon. */
    concave_polygon (std::vector<xy_t> _corners)
        :
        nx(_corners.size()),
        ny(_corners.size()),
        n0(_corners.size()),
        corners(_corners)
    {
        //cout << "aaaaaaaaaaaa" << std::flush << std::endl;
        // _corners = conv1(corners_x, corners_y);

        /*
        int nc = corners.size();

        std::vector<REAL> dx(nc);
        std::vector<REAL> dy(nc);

        std::vector<REAL> a_x(nc);  // y = a_y * x + b_y
        std::vector<REAL> b_x(nc);

        std::vector<REAL> a_y(nc);  // x = a_x * y + b_x
        std::vector<REAL> b_y(nc);

        std::vector<bool> whichx(nc);  // if |dy| > |dx|
        */

        /*
        std::vector<REAL> nx(nc);
        std::vector<REAL> ny(nc);
        std::vector<REAL> n0(nc);
        */
        //cout << "llllllll nc = " << nc << std::flush << std::endl;

        this->invalidate = false;
        this->update_inner_data();

        my_assert(this->integrity_invariant(), "");
    }

    ~concave_polygon () {
    }

public:
    /* no change in size*/
    void update_inner_data() {
        int nc = corners.size();
        for (int i = 0; i < nc; ++i) {

            // cout << ">>>>" << i << std::endl << std::flush ;

            // int next_i = (i < nc)? i+1 : 0;
            int next_i = (i < nc-1)? i+1 : 0;
            REAL dx = corners[next_i].x - corners[i].x;
            REAL dy = corners[next_i].y - corners[i].y;

            //a_x[i] = dx[i] / dy[i];
            //b_x[i] = corners[i].y;

            //a_y[i] = dy[i] / dx[i];
            //b_y[i] = corners[i].x;

            //whichx[i] = std::abs(dy[i]) > std::abs(dx[i]);
            // if whichy => divide by y => use:   x = a_x * y + b_x

            REAL d = std::sqrt( dx*dx + dy * dy );
            REAL dinv = (d > 0.00000001)? 1.0 / d : 0.0;
            nx[i] = + dy * dinv;
            ny[i] = - dx * dinv;  // outward if counter-clockwise

            n0[i] = corners[i].x * nx[i] + corners[i].y * ny[i];
            // (nx,ny) points outwards =>
            // x * nx + y * ny - n0 < 0 ==> inside
        }
        this->invalidate = false;
        my_assert(this->integrity_invariant(), "");
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
            //cout << "Evaluating: " << x << "," << y << "  ";
            auto min_val_which_pair = this->calc_minv(x, y);
            (*f_output)[output_ctr] = min_val_which_pair.first;
            //cout << " => " << min_val_which_pair.first << " (" << (*f_output)[output_ctr]   << ")"<< std::endl;
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
    // getter/setter/ref
    REAL& getX(int corner_index) {
        this->invalidate = true;
        return corners[corner_index].x;
    }
    REAL& getY(int corner_index) {
        this->invalidate = true;
        return corners[corner_index].y;
    }
    // todo: add point, delete point
    // todo: move this to a different class that keeps a set of points. It will kow whether it is convex or not.
    // It wil have a separate interface

    bool is_counter_clockwise() const {
        //REAL last_dx = std::nan("");
        //REAL last_dy = std::nan("");
        bool outcome_counter_clockwise = true;
        // cout << std::endl << "CC: " << corners.size() << "   ";

        int nc = corners.size();
        for (int i = 0; i < nc + 1 + 5; ++i) {

            // cout << ">>>>" << i ; // << std::endl << std::flush ;

            int ci = i % nc;  // circular i
            int prev_i = (ci > 0)? (ci - 1) : (nc - 1);
            int next_i = (ci + 1) % nc;
            REAL dx = corners[ci].x - corners[prev_i].x;
            REAL dy = corners[ci].y - corners[prev_i].y;
            REAL last_dx = corners[next_i].x - corners[ci].x;
            REAL last_dy = corners[next_i].y - corners[ci].y;

            //if (i > 0)
            {
                // cross_product(const vectorized_vect& A, const vectorized_vect& B, vectorized_vect &C)
                vectorized_vect A{boost::extents[1][3]};
                vectorized_vect B{boost::extents[1][3]};
                vectorized_vect C{boost::extents[1][3]};
                A[0][0] = dx;
                A[0][1] = dy;
                A[0][2] = 0;
                B[0][0] = last_dx;
                B[0][1] = last_dy;
                B[0][2] = 0;
                cross_product(A,B,C);
                //cross = cross(dx,dy,0, last_dx, last_dy,0);
                //cout << "cross: " << cross[0] << " "<< cross[1] << " "<< cross[2] << std::endl;
                /*
                cout << "ci:" << ci << " prev_i:" << prev_i << ":  ";
                cross = A[0]; //C[0];
                cout << "-------" << ci << "----cross: "
                  << cross[0] << " "<< cross[1] << " "<< cross[2]
                  ; //<< std::endl;
              cross = B[0]; //C[0];
                cout << "-------" << ci << "----cross: "
                  << cross[0] << " "<< cross[1] << " "<< cross[2]
                  ; //<< std::endl;
              cross = C[0]; //C[0];
                cout << "-------" << ci << "----cross: "
                  << cross[0] << " "<< cross[1] << " "<< cross[2]
                  << std::endl;
                */

                auto cross = C[0];
                const int _Z = 2;
                // cout << "  cross[_Z] = " << cross[_Z] << "  ";
                if (cross[_Z] < -0.0000001 ) {
                    // return false;
                    outcome_counter_clockwise = false;
                }
            }
            // last_dx = dx;
            // last_dy = dy;
        }
        //return true;
        // cout << std::endl;
        return outcome_counter_clockwise;
    }
    bool integrity_invariant() const {
        if (n0.size() < 3) {
            return false;
        }
        for (int i = 0; i < n0.size(); ++i) {
            int next_i = (i < corners.size()-1)? i+1 : 0;
            REAL dx = corners[next_i].x - corners[i].x;
            REAL dy = corners[next_i].y - corners[i].y;

            if(std::abs(std::sqrt(dx*dx + dy*dy)) < CONFIG.MIN_PRINTABLE_LENGTH/ 10.0)
              return false;
           //todo: check convexity
           //todo: check counter-clockwise
        }
        if (!this->is_counter_clockwise()) {
            cout << "integrity_invariant() : The polygon is not counter-clockwise" << std::endl;
            return false;
        } else {
            cout << "integrity_invariant() : It's counter-clockwise" << std::endl;
            // return true;
        }
        return true;
    }

    virtual mp5_implicit::bounding_box_2d  get_boundingbox() const {
        // todo: fixme
        return bounding_box_2d{-100, 100, -100, 100 };
    }
};

}
