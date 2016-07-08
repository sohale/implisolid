/** File: basic_data_structures.hpp
 *  -------------------------------
 * This file contains data structures and functions that are repeatedly used,
 * so they should be all in a common place, with their documentation and usage
 * messages.
 */


#ifndef IMPLICIT_BASIC_DATASTRUCTURES_HPP
#define IMPLICIT_BASIC_DATASTRUCTURES_HPP


#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <math.h>
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/io.hpp"



using namespace boost::numeric::ublas;
using namespace std;

//#include "boost/assert.hpp"

//namespace implicit {


/*
================================================================
=                     Useful Data Structures                   =
================================================================
*/




typedef unsigned short int dim_t;
typedef float REAL;

typedef boost::multi_array<REAL, 1>  array1d;
typedef boost::array<array1d::index, 1>  array_shape_t;
typedef array1d::index  index_t;


typedef array1d::index  vertex_t;

/** Implementation Note
 *  Data type: vectorized_scalar
 *  -----------------------------
 *  The type vectorized_scalar is a container for floats that has one dimension
 *  or else, a vector.
 *
 */

typedef boost::multi_array<REAL, 1>  vectorized_scalar;

/** Implementation Note
 * Data type: vectorized_vect
 * --------------------------
 *
 * The type vectorized_vect is a container for floats and has two dimensions.
 * It aims to be a container for 3d coordinates, so a variable of this type
 * will be a N x 3 array.
 */

typedef boost::multi_array<REAL, 2>  vectorized_vect;


#include ".//my_assert.hpp"


inline array_shape_t make_shape_1d(array1d::index size) {
    my_assert(size>=0, "");
    array_shape_t shape = { size, };
    return shape;
}


typedef boost::multi_array<REAL, 2>  array2d;


/*
================================================================
=                       Useful Functions                       =
================================================================
*/

/**
 * Function: make_empty_x
 * Usage: boost::multi_array<REAL, 2> x = make_empty_x(100)
 * --------------------------------------------------------
 * Creates an empty array with dimensions N x 3, whose elements are of floating
 * point numbers(REAL).
 *
 */

boost::multi_array<REAL, 2>  make_empty_x(const int nsize){
    auto sf = make_shape_1d(nsize);
    //vectorized_scalar  f = vectorized_scalar(sf);

    boost::array<int, 2> values_shape = {{ nsize, 3 }};
    boost::multi_array<REAL, 2> values (values_shape);
    return values;
}

// Create a InvertMatrix function

 /* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool InvertMatrix(const matrix<T>& input, matrix<T>& inverse)
{
	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<T> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

// void Matrix_Vector_Product(const matrix<REAL>& matou, REAL *vectou){
//   vectou[0] = matou[0][0]*vectou[0] + matou[0][1]*vectou[1] + matou[0][2]*vectou[2] + matou[0][3]*vectou[3];
//   vectou[1] = matou[1][0]*vectou[0] + matou[1][1]*vectou[1] + matou[1][2]*vectou[2] + matou[1][3]*vectou[3];
//   vectou[2] = matou[2][0]*vectou[0] + matou[2][1]*vectou[1] + matou[2][2]*vectou[2] + matou[2][3]*vectou[3];
// // we do not compute vectou[3] because it is not relevant.
//
// }

namespace mp5_implicit {

    struct bounding_box {
        REAL xmin, xmax, ymin, ymax, zmin, zmax;
    };

}

/*
================================================================
=           Configuration Parameters                           =
================================================================
*/


/**
 * Configuration parameters
 * --------------------------
 * Contains global configuration options for tolerances, epsilon and other
 * parameters
 *
 * See: implicit_config.py
 */

const REAL ROOT_TOLERANCE = 0.001;  //
const REAL MEAN_PRINTABLE_LENGTH = 0.01;  //


#endif // IMPLICIT_BASIC_DATASTRUCTURES_HPP
