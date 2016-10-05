#include <boost/random/uniform_01.hpp>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int_distribution.hpp"

#include <algorithm>
#include <iostream>


#include "../basic_data_structures.hpp"
#include "../basic_functions.hpp"
//#include "../basic_algorithms.hpp"

#include "../vectorised_algorithms/make_random_pm1.hpp"


//vectorized algorithms
// valgorithmss   vectorized_algorithms  vectalg  vectalgo  num vectnum  vectcpp  multiarray arrayalg numarray vectcalc arithm looper  vectorized_calc
// mkdir algovec algovect  iter  basic_algorithms



#include "../vectorised_algorithms/add_inplace.hpp"



#include "gtest/gtest.h"

TEST(VectorizsedAlgorithmsTests, make_random_pm1__mean) {

    vectorized_vect result = mp5_implicit::vectorised_algorithms::make_random_pm1(4, 3, 1);
    for (vindex_t i = 0; i < std::min<unsigned int>(result.shape()[0], 10); ++i) {
        std::clog
            << result[i][0] << " "
            << result[i][1] << " "
            << result[i][2] << " "
            << std::endl;
    }

    unsigned int n = 400000/100;

    vectorized_vect A = mp5_implicit::vectorised_algorithms::make_random_pm1(n, 3, 1);
    REAL x=0, y=0, z=0;
    REAL x2=0, y2=0, z2=0;
    for (vindex_t i = 0; i < std::min<unsigned int>(A.shape()[0], 10); ++i) {
        x += A[i][0];
        y += A[i][1];
        z += A[i][2];
        x2 += A[i][0] * A[i][0];
        y2 += A[i][1] * A[i][1];
        z2 += A[i][2] * A[i][2];
    }
    std::clog << "mean: " << x/n << " " << y/n << " " << z/n << " " << std::endl;
    std::clog << "var: " << x2/(n-1) << " " << y2/(n-1) << " " << z2/(n-1) << " " << std::endl;
    std::clog << "std: " << std::sqrt(x2/(n-1)) << " " << std::sqrt(y2/(n-1)) << " " << std::sqrt(z2/(n-1)) << " " << std::endl;


    mp5_implicit::vectorised_algorithms::add_inplace(A, A);
    //todo: calculate mean


    EXPECT_LT( 1, 2 );
    EXPECT_GT( 2, 1 );
}


/*
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/

/*
#include "gtest/gtest.h"
TEST(VectorizsedAlgorithmsTests, make_random_pm1__mean) {

    EXPECT_LT( 1, 2 );
    EXPECT_GT( 2, 1 );
}
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
