#include "../primitives.cpp"
#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
#include "../crisp_subtract.hpp"
#include "../unit_sphere.hpp"

#include "gtest/gtest.h"

// #include "../../js_iteration_1/timer.hpp"

TEST(CrispSubtractTests, TwoSpheresOneInsideOther) {

    EXPECT_LT( 1, 2 );
    EXPECT_GT( 2, 1 );
}

/*
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/