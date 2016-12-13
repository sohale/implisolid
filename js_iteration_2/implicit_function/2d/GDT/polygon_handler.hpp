#pragma once
// #include "../../../basic_data_structures.hpp"  // REAL
#include "../basic_data_structures_2d.hpp"
#include "../../../configs.hpp"
//#include "../basic_functions_2d.hpp"

#include "../../../vectorised_algorithms/cross_product.hpp"
using mp5_implicit::vectorised_algorithms::cross_product;

namespace mp5_implicit {

// todo: For convex version, see HyperFun's review by Pasko & Adzhiev (incl. Implicit Curved Polygons - HyperFun)

// polygon ADT (Abstract Data Type)  was: polygon_adt.hpp
// Handles a set of (free) points. No validation test. Contains algorithms such as being counter-clockwise, or being convex (to be done).
// Contains datastructures, accessors and constructors to store points.
class polygon_handler {
public:
    typedef struct { REAL x, y; } xy_t;

public:
//protected:
    std::vector<xy_t> corners;
    bool consolidated = false;

public:
    static std::vector<xy_t> conv1(std::vector<REAL> corners_x, std::vector<REAL> corners_y) {
        int nc = corners_x.size();
        std::vector<xy_t> corners(nc);
        assert(corners_x.size() == corners_y.size());
        for (int i = 0; i < nc; ++i) {
            corners[i].x = corners_x[i];
            corners[i].y = corners_y[i];
        }
        return corners;
    }
    polygon_handler (std::vector<REAL> corners_x, std::vector<REAL> corners_y)
        :
        polygon_handler(conv1(corners_x, corners_y))
    {
    }

    /* Corners have to be arranged counter-clockwise, and form a convex polygon. */
    polygon_handler (std::vector<xy_t> _corners)
        :
        corners(_corners)
    {
    }


protected:
public:
    // getter/setter/ref
    REAL& getX(int corner_index) {
        this->consolidated = false;
        return corners[corner_index].x;
    }
    REAL& getY(int corner_index) {
        this->consolidated = false;
        return corners[corner_index].y;
    }
    void consolidate() {
        this->consolidated = true;
    }
    bool is_changed() {
        return !this->consolidated;
    }
    void addPoint(REAL x, REAL y) {
        corners.push_back(xy_t{x, y});
    }
    void removePoint(int i) {
        auto to_remove_iter = corners.begin() + i;
        if (to_remove_iter < corners.begin() || to_remove_iter >= corners.end()) {cerr << "Error: removing a point with index out of bounds.";}
        corners.erase(to_remove_iter);
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
                if (!(cross[_Z] > +0.0000001 )) {
                    // return false;
                    // cout << " ???" << cross[_Z] << " ??????";
                    outcome_counter_clockwise = false;
                }
                // else if (cross[_Z] < -0.0000001 ) ...
            }
            // last_dx = dx;
            // last_dy = dy;
        }
        //return true;
        //cout << std::endl;
        return outcome_counter_clockwise;
    }

    mp5_implicit::bounding_box_2d  get_boundingbox() const {
        REAL xmin, xmax, ymin, ymax;
        for (int i = 0; i < corners.size(); ++i) {
            if (i == 0 || xmin > corners[i].x) xmin = corners[i].x;
            if (i == 0 || xmax < corners[i].x) xmax = corners[i].x;
            if (i == 0 || ymin > corners[i].y) ymin = corners[i].y;
            if (i == 0 || ymax < corners[i].y) ymax = corners[i].y;
        }
        return bounding_box_2d{xmin, xmax, ymin, ymax};
    }
};

}
