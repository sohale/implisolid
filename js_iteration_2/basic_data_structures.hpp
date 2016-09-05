/** File: basic_data_structures.hpp
 *  -------------------------------
 * This file contains data structures and functions that are repeatedly used,
 * so they should be all in a common place, with their documentation and usage
 * messages.
 */


#pragma once


#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <math.h>

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/io.hpp"

#include "svd.cpp"



using namespace boost::numeric::ublas;
using namespace std;


//#include "boost/assert.hpp"

#ifdef NDEBUG
    #define NO_ASSERT true
#else
    #define ASSERT_USED true
#endif



//namespace implicit {


/*
================================================================
=                     Useful Data Structures                   =
================================================================
*/




typedef unsigned short int dim_t;

typedef float REAL;     // heavily used (can be changed to double from here at any time)

const REAL NaN = std::numeric_limits<REAL>::quiet_NaN();
//const REAL NaN = std::numeric_limits<REAL>::signaling_NaN();

inline bool isNaN(REAL x){return isnan(x);};

typedef boost::multi_array<REAL, 1>  array1d;

typedef boost::array<array1d::index, 1>  array_shape_t;

typedef array1d::index  index_t;

/* define types for vertices, faces and indexes of them */
typedef boost::multi_array<REAL, 2> verts_t;
typedef boost::multi_array<int, 2> faces_t;

typedef verts_t::index vindex_t;   // used for arrays of verts or centroids
//typedef verts_t::size_type vindex_t;
typedef verts_t::index eindex_t; // may be long-er than vindex_t, becasue each edgepair has twice number of vertices. (not in MC, but in O&B)


typedef short int bool_t;

const bool_t b_true = 1;
const bool_t b_false = 0;

//typedef boost::array<vectorized_vect::index, 2>  shape_t;


typedef array1d::index  vertex_t;

/** Implementation Note
 *  Data type: vectorized_scalar
 *  -----------------------------
 *  The type vectorized_scalar is a container for floats that has one dimension
 *  or else, a vector.
 *
 */

typedef boost::multi_array<REAL, 1>  vectorized_scalar;
typedef boost::array<int, 1>  vectorized_scalar_shape;

/** Implementation Note
 * Data type: vectorized_vect
 * --------------------------
 *
 * The type vectorized_vect is a container for floats and has two dimensions.
 * It aims to be a container for 3d coordinates, so a variable of this type
 * will be a N x 3 array.
 */

typedef boost::multi_array<REAL, 2>  vectorized_vect;

typedef boost::array<vectorized_vect::index, 2>  vectorized_vect_shape;
//typedef boost::array<int, 2>  vectorized_vect_shape;
// boost::array<unsigned int, 2>


// boost::multi_array<bool, 1>
typedef boost::multi_array<bool_t, 1>  vectorized_bool;
//typedef boost::array<unsigned int, 1>  vectorized_bool_shape;
typedef boost::array<vectorized_bool::index, 1>  vectorized_bool_shape;

// auto& loger = std::cerr;

#include ".//my_assert.hpp"


inline array_shape_t make_shape_1d(array1d::index size) {
    my_assert(size>=0, "");
    array_shape_t shape = { size, };
    return shape;
}


typedef boost::multi_array<REAL, 2>  array2d;


/*
    BOUNDING BOXES
*/

namespace mp5_implicit {

    struct bounding_box {
        REAL xmin, xmax, ymin, ymax, zmin, zmax;
    };

}

#include "configs.hpp"
