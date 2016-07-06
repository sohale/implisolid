/**
* File: test_crisp_subtract.cpp
*------------------------------
*
* In this file there are tests for the crisp subtraction operation.
*/

#include "../primitives.cpp"
#include "../basic_data_structures.hpp"
#include "../crisp_subtract.hpp"
#include "../unit_sphere.hpp"
#include "gtest/gtest.h"

// #include "../timer.hpp"

TEST(CrispSubtractTests, TwoSpheres) {

    const int nsize = 3;    // size of input vector
    auto shape_tuple = make_shape_1d(nsize);
    vectorized_scalar f = vectorized_scalar(shape_tuple);

    vectorized_vect x = make_empty_x(nsize);

    // 1st point, the origin.
    x[0][0] = 0.0;
    x[0][1] = 0.0;
    x[0][2] = 0.0;

    // // 2nd point,
    // x[0][0] = 0.0;
    // x[0][1] = 0.0;
    // x[0][2] = 0.0;

    mp5_implicit::unit_sphere s1(2.0);
    mp5_implicit::unit_sphere s2(1.3);

    CrispSubtract crs = CrispSubtract(s1,s2);
    crs.eval_implicit(x, &f);
    // see the output for the first point
    EXPECT_LT( f[0], -ROOT_TOLERANCE );
}

void linspace(REAL x_start, REAL x_end, int num_steps){

}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    std::cout << "Good bye." << std::endl;
    return RUN_ALL_TESTS();
}
