#pragma once
// #include "../../../basic_data_structures.hpp"  // REAL
#include "../basic_data_structures_2d.hpp"
#include "../../../configs.hpp"
//#include "../basic_functions_2d.hpp"

//#include "../../../vectorised_algorithms/cross_product.hpp"
//using mp5_implicit::vectorised_algorithms::cross_product;

#include "polygon_handler.hpp"

namespace mp5_implicit {

// todo: For convex version, see HyperFun's review by Pasko & Adzhiev (incl. Implicit Curved Polygons - HyperFun)

class convex_polygon : public implicit_function_2d {
public:
    // typedef struct { REAL x, y; } xy_t;
    std::vector<REAL> nx;
    std::vector<REAL> ny;
    std::vector<REAL> n0;
    bool invalidate = true;  // i.e. not consolidated
    polygon_handler polygon;  // ordered set of points (point sequence; closed)

//protected:
//    std::vector<xy_t> corners;

public:
    /*
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
    */

    convex_polygon (std::vector<REAL> corners_x, std::vector<REAL> corners_y)
        :
        //nx(corners_x.size()),
        //ny(corners_x.size()),
        //n0(corners_x.size()),
        //corners(conv1(corners_x, corners_y)),
        //convex_polygon(conv1(corners_x, corners_y))

        //polygon(corners_x, corners_y),
        //convex_polygon(polygon)
        // polygon(corners_x, corners_y),
        convex_polygon(polygon_handler{corners_x, corners_y})
    {
    }


    /* Corners have to be arranged counter-clockwise, and form a convex polygon. */
    //convex_polygon (std::vector<xy_t> _corners)
    convex_polygon (polygon_handler _polygon)
        :
        polygon(_polygon),
        nx(_polygon.corners.size()),
        ny(_polygon.corners.size()),
        n0(_polygon.corners.size())
        //corners(_corners)
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
        //polygon.consolidate

        my_assert(this->integrity_invariant(), "");
    }

    ~convex_polygon () {
    }

public:
    /* no change in size*/
    void update_inner_data() {
        int nc = polygon.corners.size();
        for (int i = 0; i < nc; ++i) {

            // cout << ">>>>" << i << std::endl << std::flush ;

            // int next_i = (i < nc)? i+1 : 0;
            int next_i = (i < nc-1)? i+1 : 0;
            REAL dx = polygon.corners[next_i].x - polygon.corners[i].x;
            REAL dy = polygon.corners[next_i].y - polygon.corners[i].y;

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

            n0[i] = polygon.corners[i].x * nx[i] + polygon.corners[i].y * ny[i];
            // (nx,ny) points outwards =>
            // x * nx + y * ny - n0 < 0 ==> inside
        }
        this->invalidate = false;
        my_assert(this->integrity_invariant(), "");
    }

    bool is_counter_clockwise() const {
        cout << "Going to test CCW" << std::endl << std::flush;
        return polygon.is_counter_clockwise();
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
            (*output)[output_ctr][0] = -nx[min_val_which_pair.second];
            (*output)[output_ctr][1] = -ny[min_val_which_pair.second];
        }
    }

    bool integrity_invariant() const {
        if (n0.size() < 3) {
            return false;
        }
        for (int i = 0; i < n0.size(); ++i) {
            int next_i = (i < polygon.corners.size()-1)? i+1 : 0;
            REAL dx = polygon.corners[next_i].x - polygon.corners[i].x;
            REAL dy = polygon.corners[next_i].y - polygon.corners[i].y;

            if(std::abs(std::sqrt(dx*dx + dy*dy)) < CONFIG.MIN_PRINTABLE_LENGTH/ 10.0)
              return false;
           //todo: check convexity
           //todo: check counter-clockwise
        }
        if (!this->polygon.is_counter_clockwise()) {
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
