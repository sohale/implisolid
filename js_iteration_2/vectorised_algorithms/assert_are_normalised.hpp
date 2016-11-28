#pragma once

#include "../basic_functions.hpp"

namespace mp5_implicit {
namespace vectorised_algorithms {


inline bool vector_are_zero(const REAL& x, const REAL& y, const REAL& z) {

    return (
        std::abs(x) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS &&
        std::abs(y) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS &&
        std::abs(z) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS
    );
}

inline bool assert_are_normalised(const vectorized_vect& A) {

    for (vindex_t i = 0, e = A.shape()[0]; i < e; i++) {
        REAL norm2 = norm_squared(A[i][0], A[i][1], A[i][2]);
        if (vector_are_zero(A[i][0], A[i][1], A[i][2]) || 
            std::abs(norm2 - 1.0) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS) {

        } else {
            clog << " assersion failed at i=" << i << " " << A[i][0] << "," << A[i][1] << "," << A[i][2] << std::endl;
            return false;
        }
    }
    return true;
}

inline vindex_t first_not_normalised(const vectorized_vect& A) {
    /*
    For debugging only
    */
    for (vindex_t i = 0, e = A.shape()[0]; i < e; i++) {
        REAL norm2 = norm_squared(A[i][0], A[i][1], A[i][2]);
        if(std::abs(norm2 - 1.0) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS)
            {}
        else
            return i;
    }
    return -1;
}


/*
inline void assert_are_normalised(const vectorized_vect& A) {

    // static_assert(ASSERT_USED, "Use assert_are_normalised() in assertion mode only. Use #if ASSERT_USED  before calling assert_are_normalised()");

    #if ASSERT_USED
    for (vindex_t i = 0, e = A.shape()[0]; i < e; i++) {
        REAL norm2 = norm_squared(A[i][0], A[i][1], A[i][2]);
        assert(std::abs(norm2 - 1.0) < mp5_implicit::vectorised_algorithms::ALL_CLOSE_EPS);
    }
    #else
        // error
    #endif
}
*/


}
}
