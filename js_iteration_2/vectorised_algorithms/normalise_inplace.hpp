#pragma once

namespace mp5_implicit {
namespace vectorised_algorithms {



inline REAL norm_2(REAL x, REAL y, REAL z){
  REAL norm = sqrt(x*x + y*y + z*z);
  return norm;
}

//move to basic functions
inline REAL norm_2_squared(REAL x, REAL y, REAL z){
  REAL norms = x*x + y*y + z*z;
  return norms;
}

/*
    Normalisation method 1:
        Leave the vectors that are too small remain as they are.

*/
void normalise_inplace(verts_t& A, REAL minimum_acceptable_norm) {
    assert(A.shape()[1] == 3);

    for (int i = 0, e = A.shape()[0]; i < e; i++) {
        REAL norm = norm_2(A[i][0], A[i][1], A[i][2]);
        // constexpr REAL MIN_NORM = 1.0;  // norm of the gradients that are zero
        #if ASSERT_USED
            const REAL debug_norm0 = norm;
            auto debug_A = A[i];
        #endif
        norm = (norm < minimum_acceptable_norm)?
            1.0:
            norm;
        REAL factor = 1.0 / norm;
        A[i][0]=A[i][0] * factor;
        A[i][1]=A[i][1] * factor;
        A[i][2]=A[i][2] * factor;
        #if ASSERT_USED
            for (int j = 0; j < 3; j++) {
                if (!(A[i][j] <= 1.0)) {
                    std::cout << A[i][j] << " factor:" << factor << " norm:" << norm <<
                        " debug: pre-norm:" << debug_norm0 << " pre-A:" << debug_A[0] << "," << debug_A[1] << "," << debug_A[2] << std::endl;
                }
                assert(A[i][j] <= 1.0);
                assert(A[i][j] >= -1.0);
            }
        #endif
    }
}

/*
    Normalisation method 2:
        Simply divide by norm. Perhaps it assumes that the length is neither zero nor near zero.
*/
void normalize_1111(verts_t & A) {
    for (int i = 0; i < A.shape()[0]; i++) {
        REAL norm = norm_2(A[i][0], A[i][1], A[i][2]);
        assert(norm != 0.);
        for (int j = 0; j < 3; j++) {
            A[i][j] = A[i][j] / norm;
            assert(A[i][j] <= 1.00);
            assert(A[i][j] >= -1.00);
        }
    }
}

/*
    Normalisation method 3:
        See randomized ones.

*/
// todo: move code here


}
}
