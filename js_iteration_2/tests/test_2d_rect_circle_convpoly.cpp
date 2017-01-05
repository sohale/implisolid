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
    const int number_of_test_points = 3;    // size of input vector
    auto shape_tuple = make_shape_1d(number_of_test_points);
    vectorized_scalar f = vectorized_scalar(shape_tuple);
    vectorized_vect x = make_empty_x_2d(number_of_test_points);

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


TEST(ConvexPolygonTests, ClockwiseNess) {
    {
    std::vector<REAL> xarray {0, 1, 2};
    std::vector<REAL> yarray {0, 1, 0};
    mp5_implicit::polygon_handler clockwise_poly (xarray, yarray);
    EXPECT_TRUE( !clockwise_poly.is_counter_clockwise() );
    }

    {
    std::vector<REAL> xarray {0, 2, 1};
    std::vector<REAL> yarray {0, 0, 1};
    mp5_implicit::polygon_handler counter_clockwise_poly (xarray, yarray);
    EXPECT_TRUE( counter_clockwise_poly.is_counter_clockwise() );
    }
}

TEST(ConvexPolygonTests, Triangle) {
    //std::vector<REAL> xarray {0, 1, 2};
    //std::vector<REAL> yarray {0, 1, 0};
    std::vector<REAL> xarray {0, 2, 1};
    std::vector<REAL> yarray {0, 0, 1};
    mp5_implicit::convex_polygon cc (xarray, yarray);


    const int number_of_test_points = 6;    // size of input vector
    auto shape_tuple = make_shape_1d(number_of_test_points);
    vectorized_scalar f = vectorized_scalar(shape_tuple);
    vectorized_vect x = make_empty_x_2d(number_of_test_points);

    std::vector<bool> should(number_of_test_points);

    cout << "x:" << std::flush << std::endl;
    x[0][0] = 0.0 + 0.1;
    x[0][1] = 0.0 + 0.1/2;
    should[0] = true;

    x[1][0] = 1.0;
    x[1][1] = 1.0 - 0.1;
    should[1] = true;

    x[2][0] = 1.0;
    x[2][1] = 0.0 + 0.1;
    should[2] = true;


    x[3][0] = 2.0 - 0.1;
    x[3][1] = 0.0 + 0.1/2;
    should[3] = true;


    x[4][0] = 1.0;
    x[4][1] = -1.0;
    should[4] = false;

    x[5][0] = 1.0;
    x[5][1] = 2.0;
    should[5] = false;

    cout << "going to evaluate:" << std::flush << std::endl;

    cc.eval_implicit(x, &f);
    for (int i = 0; i < should.size(); ++i) {
        cout << "f["<< i<<"] = " << f[i] << " (" << (should[i]?"inner":"outside") <<") " << std::endl;
    }
    cout << std::endl <<  std::flush;

    EXPECT_TRUE( cc.is_counter_clockwise() );

    for (int i = 0; i < should.size(); ++i) {
        cout << "f["<< i<<"] = " << f[i] << " (" << (should[i]?"inner":"outside") <<") " << std::endl;

        if (should[i]) {
            EXPECT_GT( f[i], +ROOT_TOLERANCE );
        } else {
            EXPECT_LT( f[i], -ROOT_TOLERANCE );
        }
    }
    cout << std::endl <<  std::flush;
}


/*
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
