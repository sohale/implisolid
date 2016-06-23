#include "boost/multi_array.hpp"
#include "boost/array.hpp"
//#include "boost/assert.hpp"


typedef unsigned short int dim_t;
typedef float REAL;

typedef boost::multi_array<REAL, 1>  array1d;
typedef boost::array<array1d::index, 1>  array_shape_t;
typedef array1d::index  index_t;


typedef array1d::index  vertex_t;
typedef boost::multi_array<REAL, 1>  vectorized_scalar;
//typedef boost::multi_array<REAL, 2>  vectorized_x;
//typedef boost::multi_array<REAL, 2>  vectorized_vect3d;
// vectorized_vect3d   vectorized_x  assert_3x
typedef boost::multi_array<REAL, 2>  vectorized_vect;


//vectorized_vect z = {{-1.,4,5}};

/*
typedef boost::multi_array<index3_t, 1>   array1d_e3;
typedef std::map<index3_t,int>  e3map_t;
*/

#include ".//my_assert.hpp"


inline array_shape_t make_shape_1d(array1d::index size) {
    my_assert(size>=0, "");
    array_shape_t shape = { size, };
    return shape;
}


typedef boost::multi_array<REAL, 2>  array2d;

/*
inline bool assert_3x_output(const array2d& x, array1d& output){
    assert(x.shape()[1] == output.shape()[0]);
    assert(x.shape()[1] == output.shape()[0]);
}
*/
/*
inline bool assert_3x(const array2d& x){
    int n = x.shape()[0];
    assert(x.shape()[1] == output.shape()[0]);
    assert(x.shape()[1] == output.shape()[0]);
}
*/


