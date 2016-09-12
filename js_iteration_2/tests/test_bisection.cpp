#include "gtest/gtest.h"

// #include "../implicit_function.hpp"
#include "../primitives.cpp"
#include "../vectorised_algorithms/make_random_pm1.hpp"
#include "../vectorised_algorithms/normalise_inplace.hpp"
#include "../vectorised_algorithms/assert_are_normalised.hpp"
#include "../implicit_vectorised_algorithms.hpp"

#include "../polygoniser/bisection.hpp"


using mp5_implicit::vectorised_algorithms::assert_are_normalised;

using mp5_implicit::vectorised_algorithms::make_random_pm1;
using mp5_implicit::vectorised_algorithms::normalise_forced_inplace;

//using mp5_implicit::vectorised_algorithms::
// using mp5_implicit;

// mp5_implicit::implicit_function
mp5_implicit::CrispSubtract  two_spheres(REAL r1, REAL r2) {


    mp5_implicit::unit_sphere s1(r1);
    mp5_implicit::unit_sphere s2(r2);

    mp5_implicit::CrispSubtract crs = mp5_implicit::CrispSubtract(s2, s1);
    return crs;
}

void smult(vectorized_vect & v, REAL ampl) {
    for( auto i = v.begin(), e = v.end(); i != e; ++i ) {
        (*i)[0] *= ampl;
        (*i)[1] *= ampl;
        (*i)[2] *= ampl;
    }
}

TEST(BisectionTests1, on_sphere1) {

    //mp5_implicit::CrispSubtract object = two_spheres(1.3, 2.0);
    REAL r1 = 1.3;
    REAL r2 = 2.0;
    REAL ROOT_TOLERANCE = 0.0001;

    mp5_implicit::unit_sphere s1(r1);
    mp5_implicit::unit_sphere s2(r2);

    mp5_implicit::CrispSubtract object = mp5_implicit::CrispSubtract(s2, s1);


    int nsize = 100; //00;

    // make_empty_x
    vectorized_vect x1 = make_random_pm1(nsize, 3, 1.0);   // uniform(-1,1)
    // todo: use gaussian
    normalise_forced_inplace(x1, 0.00001);

    assert(assert_are_normalised(x1));


    vectorized_vect x2 = x1;

    REAL a1 = 2.5;
    REAL a2 = 1.6;

    smult(x1, a1);
    smult(x2, a2);



    ASSERT_TRUE(test_if_points_are_outside(x1, object, ROOT_TOLERANCE, true));
    ASSERT_TRUE(test_if_points_are_inside(x2, object, ROOT_TOLERANCE, true));

    vectorized_vect result_x {x1};

    std::cout << "Starting the vectorised bisection:" << std::endl << std::flush;

    bisection(&object, result_x, x1, x2, ROOT_TOLERANCE );

    std::cout << "Finished the vectorised bisection." << std::endl << std::flush;

    auto shape_tuple = make_shape_1d(nsize);
    vectorized_scalar f = vectorized_scalar(shape_tuple);
    object.eval_implicit(result_x, &f);

    // see the output for the first point
    EXPECT_LT( std::abs(f[0]), ROOT_TOLERANCE );
    EXPECT_LT( std::abs(f[1]), ROOT_TOLERANCE );


    EXPECT_TRUE( ASSERT_USED );

}

/*
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
