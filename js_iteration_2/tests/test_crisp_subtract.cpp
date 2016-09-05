/**
* File: test_crisp_subtract.cpp
*------------------------------
*
* In this file there are tests for the crisp subtraction operation.
*/

#include "../primitives.cpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
#include "../crisp_subtract.hpp"
#include "../unit_sphere.hpp"
#include "gtest/gtest.h"

// #include "../timer.hpp"

TEST(CrispSubtractTests, TwoSpheresOneInsideOther) {
/**
 * Description:
 * Creates 2 spheres, takes the subtraction and tests the sign of the implicit
 * function correctness on 3 points.
 */
    const int nsize = 3;    // size of input vector
    auto shape_tuple = make_shape_1d(nsize);
    vectorized_scalar f = vectorized_scalar(shape_tuple);

    vectorized_vect x = make_empty_x(nsize);

    mp5_implicit::unit_sphere s1(2.0);
    mp5_implicit::unit_sphere s2(1.3);

    mp5_implicit::CrispSubtract crs = mp5_implicit::CrispSubtract(s1,s2);

    // 1st point x1**2 + y1**2 + z1**2 = 0,inside both spheres, should output -
    x[0][0] = 0.0;
    x[0][1] = 0.0;
    x[0][2] = 0.0;

    // 2nd point x2**2 + y2**2 + z2**2 = 6,outside both spheres,should output -
    x[1][0] = 1.0;
    x[1][1] = 1.0;
    x[1][2] = 2.0;

    // 3rd point , x3**2 + y3**2 + z3**2 = 3 , inside big, outside of small,
    // should output +
    x[2][0] = 1.0;
    x[2][1] = 1.0;
    x[2][2] = 1.0;


    crs.eval_implicit(x, &f);
    // see the output for the first point
    EXPECT_LT( f[0], -ROOT_TOLERANCE );
    EXPECT_LT( f[1], -ROOT_TOLERANCE );
    EXPECT_GT( f[2], +ROOT_TOLERANCE );
}

TEST(CrispSubtractTests, TwoSpheresOverlapping){
    /**
     * Description:
     * Creates 2 spheres, takes the subtraction and tests the sign of the
     * implicit function correctnes on 3 points. The difference with the above
     * test is that the two spheres now overlap instead of one being inside the
     * other.
     */

     const int nsize = 3;    // size of input vector
     auto shape_tuple = make_shape_1d(nsize);
     vectorized_scalar f = vectorized_scalar(shape_tuple);

     vectorized_vect x = make_empty_x(nsize);

     mp5_implicit::unit_sphere s1(1.5, -1, 0, 0);
     mp5_implicit::unit_sphere s2(1.5, 1, 0, 0);

     mp5_implicit::CrispSubtract crs = mp5_implicit::CrispSubtract(s1,s2);

     // 1st point outside both spheres, should output -
     x[0][0] = -3.0;
     x[0][1] = 0.0;
     x[0][2] = 0.0;

     // 2nd point should output +
     x[1][0] = -1.75;
     x[1][1] = 0.2;
     x[1][2] = 0.2;

     // 3rd point should output -
     x[2][0] = 0.0;
     x[2][1] = 0.0;
     x[2][2] = 0.0;


     crs.eval_implicit(x, &f);
     // see the output for the first point
     EXPECT_LT( f[0], -ROOT_TOLERANCE );
     EXPECT_GT( f[1], +ROOT_TOLERANCE );
     EXPECT_LT( f[2], -ROOT_TOLERANCE );
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
