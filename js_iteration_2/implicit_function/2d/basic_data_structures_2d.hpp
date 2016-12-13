#pragma once

#include "../../basic_data_structures.hpp"

namespace mp5_implicit {

    typedef boost::multi_array<REAL, 2>  vectorized_vect_2d;

    //struct bounding_box_2d {
    //    struct {REAL x, REAL y} min;
    //    struct {REAL x, REAL y} max;
    //};
    struct bounding_box_2d {
        REAL xmin, xmax, ymin, ymax;
    };

}
