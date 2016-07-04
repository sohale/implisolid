#ifndef IMPLICIT_BASIC_DATASTRUCTURES_HPP
#define IMPLICIT_BASIC_DATASTRUCTURES_HPP


#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <math.h>
//#include "boost/assert.hpp"

//namespace implicit {


/*
================================================================
=                       Useful Data Structures                 =
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
 *  so it can be thought of as a vector.This data structure aims to be a container
 *  for the evaluated implicit function values.
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
 * Creates an array with dimensions N x 3. The argument is of type real.
 *
 */

boost::multi_array<REAL, 2>  make_empty_x(const int nsize){
    auto sf = make_shape_1d(nsize);
    //vectorized_scalar  f = vectorized_scalar(sf);

    boost::array<int, 2> values_shape = {{ nsize, 3 }};
    boost::multi_array<REAL, 2> values (values_shape);
    return values;
}

//}
#endif // IMPLICIT_BASIC_DATASTRUCTURES_HPP
