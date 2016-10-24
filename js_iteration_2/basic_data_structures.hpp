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



//using namespace boost::numeric::ublas;
using namespace std;


//#include "boost/assert.hpp"

#ifdef NDEBUG
    #define NO_ASSERT true
#else
    #define ASSERT_USED true
#endif


#ifndef ASSERT_USED
    #define ASSERT_USED false
#endif

// Enables long print out of variables and arrays contents for debugging
#define DEBUG_VERBOSE false


//namespace implicit {


/*
================================================================
=                     Useful Data Structures                   =
================================================================
*/




typedef unsigned short int dim_t;

typedef float REAL;     // heavily used (can be changed to double from here at any time)

typedef boost::multi_array<REAL, 1>  array1d;

typedef boost::array<array1d::index, 1>  array_shape_t;

typedef array1d::index  index_t;

/* define types for vertices, faces and indexes of them */
// typedef boost::multi_array<REAL, 2> verts_t;
//typedef boost::multi_array<int, 2> faces_t;

//typedef verts_t::index vindex_t;   // used for arrays of verts or centroids
//typedef verts_t::size_type vindex_t;
//typedef verts_t::index eindex_t; // may be long-er than vindex_t, becasue each edgepair has twice number of vertices. (not in MC, but in O&B)


/*
C::element  is not C::value_type   when C is multi-dimensional
*/
//typedef  faces_t::value_type  vertexindex_type;
//typedef  faces_t::element  vertexindex_type;
// Also see vertexindex_type_


typedef short int bool_t;

const bool_t  b_true = 1;
const bool_t  b_false = 0;

//typedef boost::array<vectorized_vect::index, 2>  shape_t;


typedef array1d::index  vertex_t;

/** Implementation Note
 *  Data type: vectorized_scalar
 *  -----------------------------
 *  The type vectorized_scalar is a container for floats that has one dimension
 *  or else, a vector.
 *
 */

// Used for implicit values only
typedef boost::multi_array<REAL, 1>  vectorized_scalar;
typedef boost::array<int, 1>  vectorized_scalar_shape;

// Used for curvatures, etc.
typedef boost::multi_array<REAL, 1>  vectrorized_real;

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
typedef boost::multi_array<bool_t, 2>  vectorized_bool_2d;


// vector_of_indices   array_of_indices   indices_array
typedef boost::multi_array<int, 1>  array_of_indices;
typedef boost::array<array_of_indices::index, 1>  array_of_indices_shape;
// boost::multi_array<vindex_t,1>

class array_of_indices_struct {
    array_of_indices  array;
    array_of_indices::index  effective_length;

    array_of_indices_struct(array_of_indices_shape shape)
        :array(shape),
        effective_length(shape[0])
    {
    }
};


// typedef boost::multi_array<vectorized_vect::index, 2>  vectorized_faces;
// typedef boost::multi_array<short int, 2>  vectorized_faces;

/* Note that this type is also used in Marching Cubes. But MC should use a different versionof the type, internally.
 * short int did not work, and it may be because of this reason.
 */
typedef  int  vertexindex_type;
// history: int, short int, int

typedef boost::multi_array<vertexindex_type, 2>  vectorized_faces;

typedef boost::array<vectorized_faces::index, 2>  vectorized_faces_shape;

// typedef vectorized_vect::index  vertexindex_type_;
//typedef vectorized_faces::element  vertexindex_type;

typedef vectorized_faces::index  faceindex_type;

//assert_static(vertexindex_type_ == vertexindex_type);
//assert_static(vindex_t == vertexindex_type);


// The following are used in MarchingCubes
typedef vectorized_vect::index vindex_t;   // used for arrays of verts or centroids
// wrong:
typedef vectorized_vect::index eindex_t; // may be long-er than vindex_t, becasue each edgepair has twice number of vertices. (not in MC, but in O&B)

/*
Refactoring todo: replace:
    vertexindex_type:
        faces_t::element -> vectorized_vect::index
    faces_t -> vectorized_faces
    verts_t -> vectorized_vect
    vertexindex_type -> vertexindex_type_
*/


// *******************************************************
/*!
    Types related to mesh subdivision and mesh_algorithm s
*/
typedef long edge_pair_type;
typedef std::map<edge_pair_type, vertexindex_type>  eulookup_map_type;
// vertexindex_type


typedef short int short_edge_type;
typedef boost::multi_array<short_edge_type, 2>   edges_of_xxx_type;


typedef boost::multi_array<faceindex_type, 2>   faces_of_xxx_type;

/* Subdivision */
namespace mp5_implicit {
namespace subdivision {

    typedef std::map<edge_pair_type, vertexindex_type> midpointmap_type;


    typedef  boost::multi_array<faceindex_type, 1>  faces_subset_type;

} // namespace subdivision
} // namespace mp5_implicit

/*
namespace mp5_implicit {
    // edgepair_Base  EdgecodeBase
    constexpr edge_pair_type  edgepair_Base = 1000000L;  // A typical sculpture on MMF has 800 K faces => 400K vertices.
}
*/


//remove this:
//typedef boost::multi_array<vectorized_vector::element_type, 1>  array_of_________indices;
//typedef boost::array<array_of_indices::index, 1>  array_of_indices_shape;

// *******************************************************

// auto& loger = std::cerr;

#include "./my_assert.hpp"


inline array_shape_t make_shape_1d(array1d::index size) {
    my_assert(size>=0, "");
    array_shape_t shape = { size, };
    return shape;
}

// don't use. deprecate.
typedef boost::multi_array<REAL, 2>  array2d;


/*
    BOUNDING BOXES
*/

namespace mp5_implicit {

    struct bounding_box {
        REAL xmin, xmax, ymin, ymax, zmin, zmax;
    };

}

//#include "configs.hpp"
