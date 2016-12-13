#pragma once

namespace mp5_implicit {

bool assert_implicit_function_io_2d(const vectorized_vect_2d& x, vectorized_scalar* f_output) {
    my_assert(x.shape()[1] == 2, "Size needs to be N x 2");
    my_assert(x.shape()[0] == f_output->shape()[0], "Sizes of the vectorised input and output arrays don't match");
    return true;
}
/* Also see vectorized_vect  make_empty_x(const int); */
vectorized_vect_2d  make_empty_x_2d(const int nsize) {
    boost::array<int, 2> values_shape = {{ nsize, 2 }};
    vectorized_vect_2d  values (values_shape);
    return values;
}
}
