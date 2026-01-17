#pragma once

#include <boost/array.hpp>

typedef float REAL;     // heavily used (can be changed to double from here at any time)

// moved from mcc2.cpp:

// typedef float REAL;

// typedef unsigned short int size_t;
// typedef unsigned short int dim_t; // From MS version/fork.
typedef uint16_t dim_t;  // small integers for example the size of one side of the grid

// typedef unsigned long int index_t;

// boost::array will not work becasue the size of a boost:array has to be known in compile-time (static).
typedef  boost::multi_array<REAL, 1>  array1d;
// typedef array1d::index  array_index_t;
typedef boost::array<array1d::index, 1>  array_shape_t;
// #define array1d  boost::multi_array<REAL, 1>
typedef array1d::index  index_t;


// MS version/fork: js_iteration_1/legacy/mcc2_MS.cpp
// type ideas: (may be removed)
typedef index_t index3_t; // Range of the element type has to be large enough, larger than (size^3)*3.
typedef boost::multi_array<index3_t, 1>   array1d_e3;
typedef std::map<index3_t,int>  e3map_t;


// from js_iteration_2/basic_data_structures.hpp

typedef unsigned short int dim_t;

// typedef float REAL;     // heavily used (can be changed to double from here at any time)

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
