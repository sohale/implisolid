#pragma once

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

bool assert_implicit_function_io_2d(const vectorized_vect_2d& x, vectorized_scalar* f_output) {
    my_assert(x.shape()[1] == 2, "Size needs to be N x 2");
    my_assert(x.shape()[0] == f_output->shape()[0], "Sizes of the vectorised input and output arrays don't match");
    return true;
}


