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

TEST(MarchingCubes, BoundingBox) {
/**
 * Description:
 */

    //json
    //do marching cubes
    // call make_geometry
    //get min and max of verts
    //check min is almost == xmin
    //repeat for a range of xmin and max
    //EXPECT_LT( f[0], -ROOT_TOLERANCE );

}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    std::cout << "Good bye." << std::endl;
    return RUN_ALL_TESTS();
}
