// matrix_operations.hpp
// linear_algebera.hpp

#include "Eigen/Core"
#include "Eigen/SVD"

using Eigen::Matrix;


    // A singular value will be considered nonzero if its value is strictly greater than
    // | singular_value |  <=  threshold x | max_singular_value |
    // accepted:
    // |singular_value| / |max_singular_value| <=  threshold

    /*
    The type of S is:
    const SingularValuesType' (aka 'const Eigen::Matrix<float, 3, 1, 0, 3, 1>')

    */
/*
inline int SVD(
        const Matrix<REAL, 3, 3> &A,
        Matrix<REAL, 3, 3> & U,
        Eigen::SingularValuesType & S,
        Matrix<REAL, 3, 3> & V,
        REAL threshold
    )
{

    Eigen::JacobiSVD<Matrix<REAL, 3, 3>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    svd.setThreshold( threshold );
    S = svd.singularValues();
    U = svd.matrixU();
    V = svd.matrixV();
    return svd.rank();
}
*/

