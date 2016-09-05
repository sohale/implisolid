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


void normalise_inplace(verts_t& A) {
    assert(A.shape()[1] == 3);

    for (int i = 0, e = A.shape()[0]; i < e; i++) {
        REAL norm = norm_2(A[i][0], A[i][1], A[i][2]);
        // constexpr REAL MIN_NORM = 1.0;  // norm of the gradients that are zero
        REAL debug_norm0 = norm;
        auto debug_A = A[i];
        norm = (norm < mp5_implicit::CONFIG_C::center_projection::min_gradient_len)? norm: 1.0;
        REAL factor = 1.0 / norm;
        A[i][0]=A[i][0] * factor;
        A[i][1]=A[i][1] * factor;
        A[i][2]=A[i][2] * factor;
        #if ASSERT_USED
            for (int j = 0; j < 3; j++) {
                if (!(A[i][j] <= 1.))
                    std::cout << A[i][j] << " factor:" << factor << " norm:" << norm <<
                " debug: pre-norm:" << debug_norm0 << " pre-A:" << debug_A[0] << "," << debug_A[1] << "," << debug_A[2] << std::endl;
                assert(A[i][j] <= 1.0);
                assert(A[i][j] >= -1.0);
            }
        #endif
    }
}

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





}
}
