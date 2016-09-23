#include "gtest/gtest.h"

#include "Eigen/Core"
#include "Eigen/SVD"
using Eigen::Matrix;

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

TEST(GeneralSparse_Eigen, CSR) {
    // matlab's csr_matrix
}

TEST(GeneralSparse_Boost, CSR) {
    /*!
       \brief csr_matrix
       Example usage of Boost's CompressedMatrix, equivalent to matlab's csr_matrix

       \see http://www.boost.org/doc/libs/1_45_0/libs/numeric/ublas/doc/matrix_sparse.htm#2CompressedMatrix
       For Doxygen, \see http://stackoverflow.com/questions/51667/best-tips-for-documenting-code-using-doxygen

       Also see  matrix_container<>
    */


    // compressed_matrix<double, row_major> is the default one
    using namespace boost::numeric::ublas;
    compressed_matrix<double> m (3, 3, 3 * 3);
    for (unsigned i = 0; i < m.size1 (); ++ i)
        for (unsigned j = 0; j < m.size2 (); ++ j)
            m (i, j) = 3 * i + j;
    std::cout << m << std::endl;
    std::cout << "m mm mm mm mm mm mm mm mm mm m" << std::endl;
}
