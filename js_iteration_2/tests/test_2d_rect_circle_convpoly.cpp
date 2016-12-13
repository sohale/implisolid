/**
*/
#include "gtest/gtest.h"

#include "../implicit_function/2d/basic_data_structures_2d.hpp"
#include "../implicit_function/2d/basic_functions_2d.hpp"
#include "../implicit_function/2d/primitives_2d.hpp"

using mp5_implicit::make_empty_x_2d;

TEST(CircleTests, CenterPoint) {
/**
 */
    const int number_of_points = 3;    // size of input vector
    auto shape_tuple = make_shape_1d(number_of_points);
    vectorized_scalar f = vectorized_scalar(shape_tuple);
    vectorized_vect x = make_empty_x_2d(number_of_points);

    // mp5_implicit::twodim::circle s1(2.0);
    mp5_implicit::circle cc(0,0, 1.3);

    x[0][0] = 0.0;
    x[0][1] = 0.0;

    x[1][0] = 1.0;
    x[1][1] = 1.0;

    x[2][0] = -1.0;
    x[2][1] = -1.0;


    cc.eval_implicit(x, &f);
    cout << "f0:" << f[0];
    cout << "f1:" << f[1];
    cout << "f2:" << f[2] << std::endl;

    EXPECT_GT( f[0], +ROOT_TOLERANCE );
    EXPECT_LT( f[1], -ROOT_TOLERANCE );
    EXPECT_LT( f[2], -ROOT_TOLERANCE );
}

TEST(ConvexPolygonTests, Triangle) {
/**
 */
    const int number_of_points = 4;    // size of input vector
    auto shape_tuple = make_shape_1d(number_of_points);
    vectorized_scalar f = vectorized_scalar(shape_tuple);
    vectorized_vect x = make_empty_x_2d(number_of_points);

    int corners = 3;
    std::vector<REAL> xarray {0, 1, 2};
    std::vector<REAL> yarray {0, 1, 0};

    cout << "constructor" << std::flush << std::endl;


    // mp5_implicit::twodim::circle s1(2.0);
    mp5_implicit::concave_polygon cc (xarray, yarray);

    cout << "x:" << std::flush << std::endl;
    x[0][0] = 0.0;
    x[0][1] = 0.0;

    x[1][0] = 1.0;
    x[1][1] = 1.0;

    x[2][0] = -1.0;
    x[2][1] = -1.0;

    x[3][0] = +1.0;
    x[3][1] = -1.0;

    cout << "going to evaluate:" << std::flush << std::endl;

    cc.eval_implicit(x, &f);
    cout << "f0:" << f[0];
    cout << "f1:" << f[1];
    cout << "f2:" << f[2] << std::endl;

    EXPECT_GT( f[0], +ROOT_TOLERANCE );
    EXPECT_LT( f[1], -ROOT_TOLERANCE );
    EXPECT_LT( f[2], -ROOT_TOLERANCE );
}


/*
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
