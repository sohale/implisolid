
#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include "../js_iteration_2/primitives.cpp"
//#include "../js_iteration_2/basic_data_structures.hpp"
//#include "../js_iteration_2/unit_sphere.hpp"

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

using namespace std;

const bool VERBOSE = false;
const bool REPORT_STATS = false;


typedef unsigned short int dim_t; //small integers for example the size of one side of the grid
typedef float REAL;

typedef  boost::multi_array<REAL, 1>  array1d;
typedef boost::array<array1d::index, 1>  array_shape_t;
typedef array1d::index  index_t;

typedef index_t index3_t; //Range of the element type has to be large enough, larger than (size^3)*3.
typedef boost::multi_array<index3_t, 1>   array1d_e3;
typedef std::map<index3_t,int>  e3map_t;


struct callback_t { void call (void*) const { } callback_t(){} };

// array_shape_t make_shape_1d(int size) {
//     array_shape_t shape = {{ size, }};
//     return shape;
// }


REAL lerp(REAL a, REAL b, REAL t ) {
    return a + ( b - a ) * t;
}


class MarchingCubes{
    bool enableUvs, enableColors;
    dim_t resolution;

    index_t  yd, zd;
    index_t  yd_global, zd_global;
    REAL halfsize;
    REAL delta;

    static const dim_t queueSize = 4096;

 protected:
    //Buffers:
    index_t  temp_buffer_size = 12;
    // temp buffers used in polygonize
    array1d vlist_buffer;
    array1d nlist_buffer;
    array1d_e3 e3list_buffer;


    //Queues:
    int queue_counter = 0;

    bool hasPositions = false;
    bool hasColors = false;
    bool hasUvs = false;

    array1d positionQueue;
    array1d_e3 e3Queue;

    array1d *colorQueue = 0;
    array1d *uvQueue = 0;

    void kill();

    static const int mc_edge_lookup_table[256];
    static const int mc_triangles_table[256*16];

 protected:
    void init(dim_t resolution);

inline void VIntX(index_t q, array1d &pout, array1d &nout,    int offset, REAL isol, REAL x, REAL y, REAL z, REAL valp1, REAL valp2,  index_t ijk, array1d_e3& e3out);
inline void VIntY(index_t q, array1d& pout, array1d& nout,    int offset, REAL isol, REAL x, REAL y, REAL z, REAL valp1, REAL valp2,  index_t ijk, array1d_e3& e3out);
inline void VIntZ(index_t q, array1d& pout, array1d& nout,    int offset, REAL isol, REAL x, REAL y, REAL z, REAL valp1, REAL valp2,  index_t ijk, array1d_e3& e3out);

void compNorm( index_t q );
void posnormtriv( array1d& pos__vlist, array1d& norm__nlist, array1d_e3& e3__e3list, int o1, int o2, int o3, const callback_t& renderCallback );

void begin_queue();
void finish_queue( const callback_t& renderCallback );


public:
    index_t size, size2, size3;
    array1d field;

    MarchingCubes( dim_t resolution, bool enableUvs, bool enableColors );
    ~MarchingCubes();

    REAL isolation;

    void flush_geometry_queue(std::ostream& cout, int& marching_cube_start, std::vector<REAL> &verts3, std::vector<int> &faces3, e3map_t &e3map, int& next_unique_vect_counter);
    void reset_result();

    int polygonize_cube( REAL fx, REAL fy, REAL fz, index_t q, REAL isol, const callback_t& callback );

//shape:
    void addBall( REAL ballx, REAL bally, REAL ballz, REAL strength, REAL subtract );
    void addPlaneX( REAL strength, REAL subtract );
    void addPlaneZ( REAL strength, REAL subtract );
    void addPlaneY( REAL strength, REAL subtract );
    void seal_exterior(const REAL exterior_value = -1.);

//field
    void reset();

//geometry/threejs interface side.
    void render_geometry(const callback_t& renderCallback );
    void sow();

// output. filled using sow()
    int resultqueue_faces_start = 0;

    std::vector<REAL> result_verts;
    std::vector<int> result_faces;
    e3map_t result_e3map;
    int next_unique_vect_counter = 0;
};


int EXCESS = 0;
MarchingCubes::MarchingCubes( dim_t resolution, bool enableUvs=false, bool enableColors=false )
    :   //constructor's initialisation list: pre-constructor code
        //All memory allocation code is here. Because the size of arrays is determined in run-time.
        field(array1d( array_shape_t ({{ resolution*resolution*resolution }}) )),
        vlist_buffer(array1d( array_shape_t( {temp_buffer_size * 3} ) )),
        e3list_buffer(array1d_e3(  make_shape_1d(temp_buffer_size)   )),

        positionQueue(array1d(make_shape_1d(MarchingCubes::queueSize * 3 + EXCESS))),
        e3Queue(array1d_e3(make_shape_1d(MarchingCubes::queueSize )))

{

    this->enableUvs = enableUvs;
    this->enableColors = enableColors;

    if(VERBOSE)
        std::cout << resolution << " init"<< std::endl;

    this->init( resolution );

}


void MarchingCubes::init( dim_t resolution ) {
        // May throw  std::bad_alloc. See #include <new>
        // init() is only called by the constructor

        this->resolution = resolution;

        this->isolation = 80.0;

        this->size = resolution;
        this->size2 = this->size * this->size;
        this->size3 = this->size2 * this->size;
        this->halfsize = ((REAL)this->size) / 2.0;

        // deltas
        this->delta = 2.0 / (REAL)this->size;
        this->yd = this->size;
        this->zd = this->size2;
        this->yd_global = this->size;
        this->zd_global = this->size2;

        array_shape_t size = {(int)this->size3};
        this->field = array1d(size);

        assert(this->size3 < 10000000);
        assert(this->size3 > 0);


        if(false){
            this->vlist_buffer = array1d( make_shape_1d( temp_buffer_size * 3 ) );
        }

        this->queue_counter = 0;

        this->hasPositions = false;
        this->hasColors = false;
        this->hasUvs = false;


        auto shape_maxCount_x_3 = make_shape_1d(MarchingCubes::queueSize * 3);
        this->positionQueue = array1d(shape_maxCount_x_3);
        auto shape_maxCount_x_1 = make_shape_1d(MarchingCubes::queueSize * 1);
        this->e3Queue = array1d_e3(shape_maxCount_x_1);

        auto shape_maxCount_x_2 = make_shape_1d(MarchingCubes::queueSize * 2);

        if ( this->enableUvs ) {

            this->uvQueue = 0;
            this->uvQueue = new array1d(shape_maxCount_x_2);

        }


        if ( this->enableColors ) {
            this->colorQueue = 0;
            this->colorQueue = new array1d(shape_maxCount_x_3);

        }

}

MarchingCubes::~MarchingCubes()
{
    if(VERBOSE)
        std::cout << "Destructor: ~MarchingCubes" << std::endl;

    if ( this->enableUvs )
    {
        if(this->uvQueue){
            delete this->uvQueue;
            this->uvQueue = 0;
        }
    }
    if ( this->enableColors )
    {
        if(this->colorQueue) // if is not necessary for colorQueue,but keep this if.
        {
            delete this->colorQueue;
            this->colorQueue = 0;

        }
    }
}

void MarchingCubes::kill()

{
    ;
}



inline void MarchingCubes:: VIntX(
    index_t q, array1d &pout, array1d &nout,
    int offset,
    REAL isol,
    REAL x, REAL y, REAL z,
    REAL valp1,
    REAL valp2,
    index_t ijk, array1d_e3& e3out )
{

    REAL mu = ( isol - valp1 ) / ( valp2 - valp1 );

    pout[ offset ]     = x + mu * this->delta;
    pout[ offset + 1 ] = y;
    pout[ offset + 2 ] = z;


    index3_t e3x = ijk*3;

    e3out[offset/3] = e3x;

}

inline void fp(){
    std::cout << "it";
}

void (*fpp)() = fp;

inline void MarchingCubes:: VIntY (index_t q, array1d& pout, array1d& nout, int offset, REAL isol, REAL x, REAL y, REAL z, REAL valp1, REAL valp2,
    index_t ijk, array1d_e3& e3out )
{

    REAL mu = ( isol - valp1 ) / ( valp2 - valp1 );

    pout[ offset ]     = x;
    pout[ offset + 1 ] = y + mu * this->delta;
    pout[ offset + 2 ] = z;

    index3_t e3x = ijk*3+1;

    e3out[offset/3] = e3x;

}

inline void MarchingCubes:: VIntZ(index_t q, array1d& pout, array1d& nout, int offset, REAL isol, REAL x, REAL y, REAL z, REAL valp1, REAL valp2,
    index_t ijk, array1d_e3& e3out )
{


    REAL mu = ( isol - valp1 ) / ( valp2 - valp1 );


    pout[ offset ]     = x;
    pout[ offset + 1 ] = y;
    pout[ offset + 2 ] = z + mu * this->delta;

    index3_t e3x = ijk*3+2;
    e3out[offset/3] = e3x;
}


// Returns total number of triangles. Fills triangles.

inline int MarchingCubes::polygonize_cube( REAL fx, REAL fy, REAL fz, index_t q, REAL isol, const callback_t& renderCallback ) {
    /** Polygonise a single cube in the grid. */

    // cache indices
    index_t qx = q + 1,
        qy = q + this->yd,
        qz = q + this->zd,
        qxy = qx + this->yd,
        qxz = qx + this->zd,
        qyz = q + this->yd + this->zd,
        qxyz = qx + this->yd + this->zd;

    index3_t ijk_0 = q;
    index3_t
             ijk_x = ijk_0 + 1 ,
             ijk_y = ijk_0 + this->yd_global ,
             ijk_z = ijk_0 + this->zd_global ,
             ijk_xy = ijk_0 + 1 + this->yd_global ,
             ijk_yz = ijk_0 + this->yd_global + this->zd_global ,
             ijk_xz = ijk_0 + 1 + this->zd_global ;

    unsigned int cubeindex = 0;

    REAL
        field0 = this->field[ q ],
        field1 = this->field[ qx ],
        field2 = this->field[ qy ],
        field3 = this->field[ qxy ],
        field4 = this->field[ qz ],
        field5 = this->field[ qxz ],
        field6 = this->field[ qyz ],
        field7 = this->field[ qxyz ];

    if ( field0 < isol ) cubeindex |= 1;
    if ( field1 < isol ) cubeindex |= 2;
    if ( field2 < isol ) cubeindex |= 8;
    if ( field3 < isol ) cubeindex |= 4;
    if ( field4 < isol ) cubeindex |= 16;
    if ( field5 < isol ) cubeindex |= 32;
    if ( field6 < isol ) cubeindex |= 128;
    if ( field7 < isol ) cubeindex |= 64;


    // if cube is entirely in/out of the surface - bail, nothing to draw

    int bits = mc_edge_lookup_table[ cubeindex ];
    if ( bits == 0x00 ) return 0;

    REAL d = this->delta,
        fx2 = fx + d,
        fy2 = fy + d,
        fz2 = fz + d;

    index_t ijk = q;

    // top of the cube

    if ( bits & 1 ) {
        this->VIntX( q * 3, this->vlist_buffer, this->nlist_buffer, 0, isol, fx, fy, fz, field0, field1,  ijk_0, this->e3list_buffer);

    }

    if ( bits & 2 ) {
        this->VIntY( qx * 3, this->vlist_buffer, this->nlist_buffer, 3, isol, fx2, fy, fz, field1, field3,  ijk_x, this->e3list_buffer);

    }

    if ( bits & 4 ) {
        this->VIntX( qy * 3, this->vlist_buffer, this->nlist_buffer, 6, isol, fx, fy2, fz, field2, field3,  ijk_y , this->e3list_buffer);

    }

    if ( bits & 8 ) {
        this->VIntY( q * 3, this->vlist_buffer, this->nlist_buffer, 9, isol, fx, fy, fz, field0, field2,  ijk_0, this->e3list_buffer);

    }

    // bottom of the cube

    if ( bits & 16 ) {
        this->VIntX( qz * 3, this->vlist_buffer, this->nlist_buffer, 12, isol, fx, fy, fz2, field4, field5,  ijk_z, this->e3list_buffer);

    }

    if ( bits & 32 ) {
        this->VIntY( qxz * 3,  this->vlist_buffer, this->nlist_buffer, 15, isol, fx2, fy, fz2, field5, field7,  ijk_xz, this->e3list_buffer );

    }

    if ( bits & 64 ) {
        this->VIntX( qyz * 3, this->vlist_buffer, this->nlist_buffer, 18, isol, fx, fy2, fz2, field6, field7,  ijk_yz, this->e3list_buffer );

    }

    if ( bits & 128 ) {
        this->VIntY( qz * 3,  this->vlist_buffer, this->nlist_buffer, 21, isol, fx, fy, fz2, field4, field6,  ijk_z, this->e3list_buffer );

    }

    // vertical lines of the cube

    if ( bits & 256 ) {
        this->VIntZ( q * 3, this->vlist_buffer, this->nlist_buffer, 24, isol, fx, fy, fz, field0, field4,  ijk_0, this->e3list_buffer);

    }

    if ( bits & 512 ) {
        this->VIntZ( qx * 3,  this->vlist_buffer, this->nlist_buffer, 27, isol, fx2, fy,  fz, field1, field5,  ijk_x, this->e3list_buffer );

    }

    if ( bits & 1024 ) {
        this->VIntZ( qxy * 3, this->vlist_buffer, this->nlist_buffer, 30, isol, fx2, fy2, fz, field3, field7,  ijk_xy, this->e3list_buffer);

    }

    if ( bits & 2048 ) {
        this->VIntZ( qy * 3, this->vlist_buffer, this->nlist_buffer, 33, isol, fx,  fy2, fz, field2, field6,  ijk_y, this->e3list_buffer );

    }

    cubeindex <<= 4;

    int o1, o2, o3, numtris = 0, i = 0;

    // here is where triangles are created

    while ( mc_triangles_table[ cubeindex + i ] != - 1 ) {
        o1 = cubeindex + i;
        o2 = o1 + 1;
        o3 = o1 + 2;

        //stores the triangles into the buffers
        this->posnormtriv(
            this->vlist_buffer, this->nlist_buffer, this->e3list_buffer,
            3 * MarchingCubes::mc_triangles_table[ o1 ],
            3 * MarchingCubes::mc_triangles_table[ o2 ],
            3 * MarchingCubes::mc_triangles_table[ o3 ],
            renderCallback );

        i += 3;
        numtris++;
    }
    return numtris;
}

#define DEBUG_PA001(positionQueue , c)   {std::cout << " >" << positionQueue[ (c) ] << positionQueue[ (c) + 1 ] <<    positionQueue[ (c) + 2 ] << "< ";}

/////////////////////////////////////
// Immediate-render mode simulator
/////////////////////////////////////

void MarchingCubes::posnormtriv(
    array1d& pos__vlist, array1d& norm__nlist, array1d_e3& e3__e3list,
    int o1, int o2, int o3,
    const callback_t& renderCallback ) {

    int c = this->queue_counter * 3;

    // positions
    this->positionQueue[ c ]     = pos__vlist[ o1 ];
    this->positionQueue[ c + 1 ] = pos__vlist[ o1 + 1 ];
    this->positionQueue[ c + 2 ] = pos__vlist[ o1 + 2 ];

    this->positionQueue[ c + 3 ] = pos__vlist[ o2 ];
    this->positionQueue[ c + 4 ] = pos__vlist[ o2 + 1 ];
    this->positionQueue[ c + 5 ] = pos__vlist[ o2 + 2 ];

    this->positionQueue[ c + 6 ] = pos__vlist[ o3 ];
    this->positionQueue[ c + 7 ] = pos__vlist[ o3 + 1 ];
    this->positionQueue[ c + 8 ] = pos__vlist[ o3 + 2 ];

    int c_div_3 = this->queue_counter;
    this->e3Queue[ c_div_3 ]     = e3__e3list[ o1/3 ];
    this->e3Queue[ c_div_3 + 1 ] = e3__e3list[ o2/3 ];
    this->e3Queue[ c_div_3 + 2 ] = e3__e3list[ o3/3 ];

    // uvs

    if ( this->enableUvs ) {
        int d = this->queue_counter * 2;

        (*this->uvQueue)[ d ]     = pos__vlist[ o1 ];
        (*this->uvQueue)[ d + 1 ] = pos__vlist[ o1 + 2 ];

        (*this->uvQueue)[ d + 2 ] = pos__vlist[ o2 ];
        (*this->uvQueue)[ d + 3 ] = pos__vlist[ o2 + 2 ];

        (*this->uvQueue)[ d + 4 ] = pos__vlist[ o3 ];
        (*this->uvQueue)[ d + 5 ] = pos__vlist[ o3 + 2 ];
    }

    // colors

    if ( this->enableColors ) {
        (*this->colorQueue)[ c ]     = pos__vlist[ o1 ];
        (*this->colorQueue)[ c + 1 ] = pos__vlist[ o1 + 1 ];
        (*this->colorQueue)[ c + 2 ] = pos__vlist[ o1 + 2 ];

        (*this->colorQueue)[ c + 3 ] = pos__vlist[ o2 ];
        (*this->colorQueue)[ c + 4 ] = pos__vlist[ o2 + 1 ];
        (*this->colorQueue)[ c + 5 ] = pos__vlist[ o2 + 2 ];

        (*this->colorQueue)[ c + 6 ] = pos__vlist[ o3 ];
        (*this->colorQueue)[ c + 7 ] = pos__vlist[ o3 + 1 ];
        (*this->colorQueue)[ c + 8 ] = pos__vlist[ o3 + 2 ];
    }

    this->queue_counter += 3;

    if ( this->queue_counter >= this->queueSize - 3 ) {
        this->hasPositions = true;
        if ( this->enableUvs ) {
            this->hasUvs = true;
        }
        if ( this->enableColors ) {
            this->hasColors = true;
        }
        renderCallback.call( (void*)this );
        this->sow();
    }
}


void MarchingCubes::sow() {

    this->flush_geometry_queue(std::cout, this->resultqueue_faces_start, this->result_verts, this->result_faces,  this->result_e3map, this->next_unique_vect_counter);
}

void MarchingCubes::begin_queue() {
    /** resets the queue. */
    this->queue_counter = 0;

    this->hasPositions = false;
    this->hasUvs = false;
    this->hasColors = false;
}

//
void MarchingCubes::finish_queue( const callback_t& renderCallback ) {

    if ( this->queue_counter == 0 ) return;

    std::fill(this->positionQueue.begin() + (this->queue_counter * 3), this->positionQueue.end(), 0.0 );

    this->hasPositions = true;

    if ( this->enableUvs ) {
        this->hasUvs = true;
    }

    if ( this->enableColors ) {
        this->hasColors = true;
    }

    renderCallback.call(this);
    sow();
}



/////////////////////////////////////
// Metaballs
/////////////////////////////////////


void MarchingCubes::addBall(
        REAL ballx, REAL bally, REAL ballz,
        REAL strength, REAL subtract) {
    REAL radius = this->size * sqrt(strength / subtract);

    REAL
        zs = ballz * this->size,
        ys = bally * this->size,
        xs = ballx * this->size;

    int min_z = floor( zs - radius ); if ( min_z < 1 ) min_z = 1;
    int max_z = floor( zs + radius ); if ( max_z > this->size - 1 ) max_z = this->size - 1;
    int min_y = floor( ys - radius ); if ( min_y < 1 ) min_y = 1;
    int max_y = floor( ys + radius ); if ( max_y > this->size - 1 ) max_y = this->size - 1;
    int min_x = floor( xs - radius ); if ( min_x < 1  ) min_x = 1;
    int max_x = floor( xs + radius ); if ( max_x > this->size - 1 ) max_x = this->size - 1;


    int x, y, z;
    REAL fx, fy, fz, fz2, fy2, val;
    int y_offset, z_offset;

    for ( z = min_z; z < max_z; z++ ) {

        z_offset = this->size2 * z,
        fz = z / (REAL)this->size - ballz,
        fz2 = fz * fz;

        for ( y = min_y; y < max_y; y++ ) {

            y_offset = z_offset + this->size * y;
            fy = y / (REAL)this->size - bally;
            fy2 = fy * fy;

            for ( x = min_x; x < max_x; x++ ) {

                fx = x / (REAL)this->size - ballx;
                val = strength / ( (REAL)0.000001 + fx * fx + fy2 + fz2 ) - subtract;
                if ( val > 0.0 ) this->field[ y_offset + x ] += val;
            }
        }
    }
}

void MarchingCubes::addPlaneX(REAL strength, REAL subtract ) {
    int x, y, z;
    REAL val;
    REAL xx, xdiv;
    int cxy;

    // cache attribute lookups
    int yd = this->yd;
    int size = this->size;
    int zd = this->zd;
    array1d& field = this->field;
    REAL dist = size * sqrt(strength / (REAL)subtract);

    if ( dist > size ) dist = size;
    for ( x = 0; x < dist; x++ ) {
        xdiv = x / (REAL)size;
        xx = xdiv * xdiv;
        val = strength / (REAL)( 0.0001 + xx ) - subtract;
        if ( val > 0.0 ) {
            for ( y = 0; y < size; y++ ) {
                cxy = x + y * yd;
                for ( z = 0; z < size; z++ ) {
                    field[ zd * z + cxy ] += val;
                }
            }
        }
    }
}

void MarchingCubes::seal_exterior(const REAL exterior_value) {

    int x, y, z;
    int cxy;

    // cache attribute lookups
    int yd = this->yd;
    int size = this->size;
    int zd = this->zd;
    array1d& field = this->field;

    for ( x = 0; x < size; x++ ) {
        for ( y = 0; y < size; y++ ) {
            cxy = x + y * yd;
            for ( z = 0; z < size; z++ ) {
                bool border = (x == 0) || (x == size-1) || (y == 0) || (y == size-1) || (z == 0) || (z == size-1);
                if(border){
                    field[ zd * z + cxy ] = exterior_value;
                }
                if (z == 4 && x == 2)
                    field[ zd * z + cxy ] = +2.;
            }
        }
    }
}

void MarchingCubes::addPlaneY(REAL strength, REAL subtract ) {
    int x, y, z;
    REAL yy;
    REAL val;
    REAL ydiv;
    int cy;
    int cxy;

    // cache attribute lookups
    int size = this->size;
    int yd = this->yd;
    int zd = this->zd;
    array1d& field = this->field;
    REAL dist = size * sqrt(strength / subtract);

    if ( dist > size ) dist = size;

    for ( y = 0; y < dist; y++ ) {
        ydiv = y / (REAL)size;
        yy = ydiv * ydiv;
        val = strength / (REAL)( 0.0001 + yy ) - subtract;
        if ( val > 0.0 ) {
            cy = y * yd;
            for ( x = 0; x < size; x++ ) {
                cxy = cy + x;
                for ( z = 0; z < size; z++ )
                    field[ zd * z + cxy ] += val;
            }
        }
    }
}

void MarchingCubes::addPlaneZ( REAL strength, REAL subtract )
{
    int x, y, z;
    REAL zz, val, zdiv;
    int cz, cyz;

    // cache attribute lookups
    int size = this->size;
    int yd = this->yd;
    int zd = this->zd;
    array1d& field = this->field;
    REAL dist = size * sqrt( strength / subtract );

    if ( dist > size ) dist = size;
    for ( z = 0; z < dist; z++ ) {
        zdiv = z / (REAL)size;
        zz = zdiv * zdiv;
        val = strength / (REAL)( 0.0001 + zz ) - subtract;
        if ( val > 0.0 ) {
            cz = zd * z;
            for ( y = 0; y < size; y++ ) {
                cyz = cz + y * yd;
                for ( x = 0; x < size; x++ )
                    field[ cyz + x ] += val;
            }
        }
    }
}



/////////////////////////////////////
// Updates
/////////////////////////////////////

void MarchingCubes::reset()
{

    for (int i = 0; i < this->size3; i++ ) {
        this->field[ i ] = 0.0;
    }
}


void MarchingCubes::reset_result() {

    this->next_unique_vect_counter = 0;

    //preallocate
    int expected_vertices = 10;
    int expected_faces = 10;
    this->result_verts.reserve(expected_vertices*3);
    this->result_faces.reserve(expected_faces*3);


}

void MarchingCubes::render_geometry(const callback_t& renderCallback ) {
    this->reset_result();  //receiver of the queue
    this->begin_queue();

    int smin2 = this->size - 2;

    for ( int z = 1; z < smin2; z++ ) {

        index_t z_offset = this->size2 * z;
        REAL fz = ( z - this->halfsize ) / (REAL)this->halfsize; //+ 1

        for ( int y = 1; y < smin2; y++ ) {

            index_t y_offset = z_offset + this->size * y;
            REAL fy = ( y - this->halfsize ) / (REAL)this->halfsize; //+ 1

            for ( int x = 1; x < smin2; x++ ) {

                REAL fx = ( x - this->halfsize ) / (REAL)this->halfsize; //+ 1
                index_t q = y_offset + x;

                this->polygonize_cube( fx, fy, fz, q, this->isolation, renderCallback );

            }
        }
    }
    this->finish_queue(renderCallback);
}

/////////////////////////////////////
// Marching cubes lookup tables
/////////////////////////////////////

const int MarchingCubes::mc_edge_lookup_table[256] = {
    0x000, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x099, 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x033, 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0x0aa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x066, 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0x0ff, 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x055, 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0x0cc,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0x0cc, 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x055, 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0x0ff, 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x066, 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0x0aa, 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x033, 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x099, 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x000
};

const int MarchingCubes::mc_triangles_table[256*16] = {
                                          -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 1, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 8, 3,   9, 8, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,10,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 3,   1, 2,10,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 9, 2,10,   0, 2, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 2, 8, 3,   2,10, 8,  10, 9, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3,11, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0,11, 2,   8,11, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 9, 0,   2, 3,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1,11, 2,   1, 9,11,   9, 8,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3,10, 1,  11,10, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0,10, 1,   0, 8,10,   8,11,10,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3, 9, 0,   3,11, 9,  11,10, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 8,10,  10, 8,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 7, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 3, 0,   7, 3, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 1, 9,   8, 4, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 1, 9,   4, 7, 1,   7, 3, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,10,   8, 4, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 3, 4, 7,   3, 0, 4,   1, 2,10,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 2,10,   9, 0, 2,   8, 4, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
 2,10, 9,   2, 9, 7,   2, 7, 3,   7, 9, 4,                                            -1,-1,-1,  -1,
 8, 4, 7,   3,11, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
11, 4, 7,  11, 2, 4,   2, 0, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 0, 1,   8, 4, 7,   2, 3,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 4, 7,11,   9, 4,11,   9,11, 2,   9, 2, 1,                                            -1,-1,-1,  -1,
 3,10, 1,   3,11,10,   7, 8, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1,11,10,   1, 4,11,   1, 0, 4,   7,11, 4,                                            -1,-1,-1,  -1,
 4, 7, 8,   9, 0,11,   9,11,10,  11, 0, 3,                                            -1,-1,-1,  -1,
 4, 7,11,   4,11, 9,   9,11,10,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 5, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 9, 5, 4,   0, 8, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 5, 4,   1, 5, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 8, 5, 4,   8, 3, 5,   3, 1, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,10,   9, 5, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 3, 0, 8,   1, 2,10,   4, 9, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5, 2,10,   5, 4, 2,   4, 0, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,
 2,10, 5,   3, 2, 5,   3, 5, 4,   3, 4, 8,                                            -1,-1,-1,  -1,
 9, 5, 4,   2, 3,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0,11, 2,   0, 8,11,   4, 9, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 5, 4,   0, 1, 5,   2, 3,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 2, 1, 5,   2, 5, 8,   2, 8,11,   4, 8, 5,                                            -1,-1,-1,  -1,
10, 3,11,  10, 1, 3,   9, 5, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,
 4, 9, 5,   0, 8, 1,   8,10, 1,   8,11,10,                                            -1,-1,-1,  -1,
 5, 4, 0,   5, 0,11,   5,11,10,  11, 0, 3,                                            -1,-1,-1,  -1,
 5, 4, 8,   5, 8,10,  10, 8,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 7, 8,   5, 7, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 9, 3, 0,   9, 5, 3,   5, 7, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 7, 8,   0, 1, 7,   1, 5, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 5, 3,   3, 5, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 9, 7, 8,   9, 5, 7,  10, 1, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 1, 2,   9, 5, 0,   5, 3, 0,   5, 7, 3,                                            -1,-1,-1,  -1,
 8, 0, 2,   8, 2, 5,   8, 5, 7,  10, 5, 2,                                            -1,-1,-1,  -1,
 2,10, 5,   2, 5, 3,   3, 5, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
 7, 9, 5,   7, 8, 9,   3,11, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 5, 7,   9, 7, 2,   9, 2, 0,   2, 7,11,                                            -1,-1,-1,  -1,
 2, 3,11,   0, 1, 8,   1, 7, 8,   1, 5, 7,                                            -1,-1,-1,  -1,
11, 2, 1,  11, 1, 7,   7, 1, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 5, 8,   8, 5, 7,  10, 1, 3,  10, 3,11,                                            -1,-1,-1,  -1,
 5, 7, 0,   5, 0, 9,   7,11, 0,   1, 0,10,  11,10, 0,                                            -1,
11,10, 0,  11, 0, 3,  10, 5, 0,   8, 0, 7,   5, 7, 0,                                            -1,
11,10, 5,   7,11, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
10, 6, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 3,   5,10, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 9, 0, 1,   5,10, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 8, 3,   1, 9, 8,   5,10, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 6, 5,   2, 6, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 6, 5,   1, 2, 6,   3, 0, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 6, 5,   9, 0, 6,   0, 2, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5, 9, 8,   5, 8, 2,   5, 2, 6,   3, 2, 8,                                            -1,-1,-1,  -1,
 2, 3,11,  10, 6, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
11, 0, 8,  11, 2, 0,  10, 6, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 1, 9,   2, 3,11,   5,10, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5,10, 6,   1, 9, 2,   9,11, 2,   9, 8,11,                                            -1,-1,-1,  -1,
 6, 3,11,   6, 5, 3,   5, 1, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 8,11,   0,11, 5,   0, 5, 1,   5,11, 6,                                            -1,-1,-1,  -1,
 3,11, 6,   0, 3, 6,   0, 6, 5,   0, 5, 9,                                            -1,-1,-1,  -1,
 6, 5, 9,   6, 9,11,  11, 9, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5,10, 6,   4, 7, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 3, 0,   4, 7, 3,   6, 5,10,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 9, 0,   5,10, 6,   8, 4, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 6, 5,   1, 9, 7,   1, 7, 3,   7, 9, 4,                                            -1,-1,-1,  -1,
 6, 1, 2,   6, 5, 1,   4, 7, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 2, 5,   5, 2, 6,   3, 0, 4,   3, 4, 7,                                            -1,-1,-1,  -1,
 8, 4, 7,   9, 0, 5,   0, 6, 5,   0, 2, 6,                                            -1,-1,-1,  -1,
 7, 3, 9,   7, 9, 4,   3, 2, 9,   5, 9, 6,   2, 6, 9,                                            -1,
 3,11, 2,   7, 8, 4,  10, 6, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5,10, 6,   4, 7, 2,   4, 2, 0,   2, 7,11,                                            -1,-1,-1,  -1,
 0, 1, 9,   4, 7, 8,   2, 3,11,   5,10, 6,                                            -1,-1,-1,  -1,
 9, 2, 1,   9,11, 2,   9, 4,11,   7,11, 4,   5,10, 6,                                            -1,
 8, 4, 7,   3,11, 5,   3, 5, 1,   5,11, 6,                                            -1,-1,-1,  -1,
 5, 1,11,   5,11, 6,   1, 0,11,   7,11, 4,   0, 4,11,                                            -1,
 0, 5, 9,   0, 6, 5,   0, 3, 6,  11, 6, 3,   8, 4, 7,                                            -1,
 6, 5, 9,   6, 9,11,   4, 7, 9,   7,11, 9,                                            -1,-1,-1,  -1,
10, 4, 9,   6, 4,10,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4,10, 6,   4, 9,10,   0, 8, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 0, 1,  10, 6, 0,   6, 4, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,
 8, 3, 1,   8, 1, 6,   8, 6, 4,   6, 1,10,                                            -1,-1,-1,  -1,
 1, 4, 9,   1, 2, 4,   2, 6, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3, 0, 8,   1, 2, 9,   2, 4, 9,   2, 6, 4,                                            -1,-1,-1,  -1,
 0, 2, 4,   4, 2, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 8, 3, 2,   8, 2, 4,   4, 2, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 4, 9,  10, 6, 4,  11, 2, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 2,   2, 8,11,   4, 9,10,   4,10, 6,                                            -1,-1,-1,  -1,
 3,11, 2,   0, 1, 6,   0, 6, 4,   6, 1,10,                                            -1,-1,-1,  -1,
 6, 4, 1,   6, 1,10,   4, 8, 1,   2, 1,11,   8,11, 1,                                            -1,
 9, 6, 4,   9, 3, 6,   9, 1, 3,  11, 6, 3,                                            -1,-1,-1,  -1,
 8,11, 1,   8, 1, 0,  11, 6, 1,   9, 1, 4,   6, 4, 1,                                            -1,
 3,11, 6,   3, 6, 0,   0, 6, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,
 6, 4, 8,  11, 6, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 7,10, 6,   7, 8,10,   8, 9,10,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 7, 3,   0,10, 7,   0, 9,10,   6, 7,10,                                            -1,-1,-1,  -1,
10, 6, 7,   1,10, 7,   1, 7, 8,   1, 8, 0,                                            -1,-1,-1,  -1,
10, 6, 7,  10, 7, 1,   1, 7, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 2, 6,   1, 6, 8,   1, 8, 9,   8, 6, 7,                                            -1,-1,-1,  -1,
 2, 6, 9,   2, 9, 1,   6, 7, 9,   0, 9, 3,   7, 3, 9,                                            -1,
 7, 8, 0,   7, 0, 6,   6, 0, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,
 7, 3, 2,   6, 7, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 2, 3,11,  10, 6, 8,  10, 8, 9,   8, 6, 7,                                            -1,-1,-1,  -1,
 2, 0, 7,   2, 7,11,   0, 9, 7,   6, 7,10,   9,10, 7,                                            -1,
 1, 8, 0,   1, 7, 8,   1,10, 7,   6, 7,10,   2, 3,11,                                            -1,
11, 2, 1,  11, 1, 7,  10, 6, 1,   6, 7, 1,                                            -1,-1,-1,  -1,
 8, 9, 6,   8, 6, 7,   9, 1, 6,  11, 6, 3,   1, 3, 6,                                            -1,
 0, 9, 1,  11, 6, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 7, 8, 0,   7, 0, 6,   3,11, 0,  11, 6, 0,                                            -1,-1,-1,  -1,
 7,11, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 7, 6,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 3, 0, 8,  11, 7, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 1, 9,  11, 7, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 8, 1, 9,   8, 3, 1,  11, 7, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 1, 2,   6,11, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,10,   3, 0, 8,   6,11, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
 2, 9, 0,   2,10, 9,   6,11, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
 6,11, 7,   2,10, 3,  10, 8, 3,  10, 9, 8,                                            -1,-1,-1,  -1,
 7, 2, 3,   6, 2, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 7, 0, 8,   7, 6, 0,   6, 2, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,
 2, 7, 6,   2, 3, 7,   0, 1, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 6, 2,   1, 8, 6,   1, 9, 8,   8, 7, 6,                                            -1,-1,-1,  -1,
10, 7, 6,  10, 1, 7,   1, 3, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 7, 6,   1, 7,10,   1, 8, 7,   1, 0, 8,                                            -1,-1,-1,  -1,
 0, 3, 7,   0, 7,10,   0,10, 9,   6,10, 7,                                            -1,-1,-1,  -1,
 7, 6,10,   7,10, 8,   8,10, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,
 6, 8, 4,  11, 8, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 3, 6,11,   3, 0, 6,   0, 4, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
 8, 6,11,   8, 4, 6,   9, 0, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 4, 6,   9, 6, 3,   9, 3, 1,  11, 3, 6,                                            -1,-1,-1,  -1,
 6, 8, 4,   6,11, 8,   2,10, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,10,   3, 0,11,   0, 6,11,   0, 4, 6,                                            -1,-1,-1,  -1,
 4,11, 8,   4, 6,11,   0, 2, 9,   2,10, 9,                                            -1,-1,-1,  -1,
10, 9, 3,  10, 3, 2,   9, 4, 3,  11, 3, 6,   4, 6, 3,                                            -1,
 8, 2, 3,   8, 4, 2,   4, 6, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 4, 2,   4, 6, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 9, 0,   2, 3, 4,   2, 4, 6,   4, 3, 8,                                            -1,-1,-1,  -1,
 1, 9, 4,   1, 4, 2,   2, 4, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
 8, 1, 3,   8, 6, 1,   8, 4, 6,   6,10, 1,                                            -1,-1,-1,  -1,
10, 1, 0,  10, 0, 6,   6, 0, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,
 4, 6, 3,   4, 3, 8,   6,10, 3,   0, 3, 9,  10, 9, 3,                                            -1,
10, 9, 4,   6,10, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 9, 5,   7, 6,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 3,   4, 9, 5,  11, 7, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5, 0, 1,   5, 4, 0,   7, 6,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
11, 7, 6,   8, 3, 4,   3, 5, 4,   3, 1, 5,                                            -1,-1,-1,  -1,
 9, 5, 4,  10, 1, 2,   7, 6,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 6,11, 7,   1, 2,10,   0, 8, 3,   4, 9, 5,                                            -1,-1,-1,  -1,
 7, 6,11,   5, 4,10,   4, 2,10,   4, 0, 2,                                            -1,-1,-1,  -1,
 3, 4, 8,   3, 5, 4,   3, 2, 5,  10, 5, 2,  11, 7, 6,                                            -1,
 7, 2, 3,   7, 6, 2,   5, 4, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 5, 4,   0, 8, 6,   0, 6, 2,   6, 8, 7,                                            -1,-1,-1,  -1,
 3, 6, 2,   3, 7, 6,   1, 5, 0,   5, 4, 0,                                            -1,-1,-1,  -1,
 6, 2, 8,   6, 8, 7,   2, 1, 8,   4, 8, 5,   1, 5, 8,                                            -1,
 9, 5, 4,  10, 1, 6,   1, 7, 6,   1, 3, 7,                                            -1,-1,-1,  -1,
 1, 6,10,   1, 7, 6,   1, 0, 7,   8, 7, 0,   9, 5, 4,                                            -1,
 4, 0,10,   4,10, 5,   0, 3,10,   6,10, 7,   3, 7,10,                                            -1,
 7, 6,10,   7,10, 8,   5, 4,10,   4, 8,10,                                            -1,-1,-1,  -1,
 6, 9, 5,   6,11, 9,  11, 8, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3, 6,11,   0, 6, 3,   0, 5, 6,   0, 9, 5,                                            -1,-1,-1,  -1,
 0,11, 8,   0, 5,11,   0, 1, 5,   5, 6,11,                                            -1,-1,-1,  -1,
 6,11, 3,   6, 3, 5,   5, 3, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,10,   9, 5,11,   9,11, 8,  11, 5, 6,                                            -1,-1,-1,  -1,
 0,11, 3,   0, 6,11,   0, 9, 6,   5, 6, 9,   1, 2,10,                                            -1,
11, 8, 5,  11, 5, 6,   8, 0, 5,  10, 5, 2,   0, 2, 5,                                            -1,
 6,11, 3,   6, 3, 5,   2,10, 3,  10, 5, 3,                                            -1,-1,-1,  -1,
 5, 8, 9,   5, 2, 8,   5, 6, 2,   3, 8, 2,                                            -1,-1,-1,  -1,
 9, 5, 6,   9, 6, 0,   0, 6, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,
 1, 5, 8,   1, 8, 0,   5, 6, 8,   3, 8, 2,   6, 2, 8,                                            -1,
 1, 5, 6,   2, 1, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 3, 6,   1, 6,10,   3, 8, 6,   5, 6, 9,   8, 9, 6,                                            -1,
10, 1, 0,  10, 0, 6,   9, 5, 0,   5, 6, 0,                                            -1,-1,-1,  -1,
 0, 3, 8,   5, 6,10,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
10, 5, 6,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
11, 5,10,   7, 5,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
11, 5,10,  11, 7, 5,   8, 3, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5,11, 7,   5,10,11,   1, 9, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,
10, 7, 5,  10,11, 7,   9, 8, 1,   8, 3, 1,                                            -1,-1,-1,  -1,
11, 1, 2,  11, 7, 1,   7, 5, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 3,   1, 2, 7,   1, 7, 5,   7, 2,11,                                            -1,-1,-1,  -1,
 9, 7, 5,   9, 2, 7,   9, 0, 2,   2,11, 7,                                            -1,-1,-1,  -1,
 7, 5, 2,   7, 2,11,   5, 9, 2,   3, 2, 8,   9, 8, 2,                                            -1,
 2, 5,10,   2, 3, 5,   3, 7, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 8, 2, 0,   8, 5, 2,   8, 7, 5,  10, 2, 5,                                            -1,-1,-1,  -1,
 9, 0, 1,   5,10, 3,   5, 3, 7,   3,10, 2,                                            -1,-1,-1,  -1,
 9, 8, 2,   9, 2, 1,   8, 7, 2,  10, 2, 5,   7, 5, 2,                                            -1,
 1, 3, 5,   3, 7, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 7,   0, 7, 1,   1, 7, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 0, 3,   9, 3, 5,   5, 3, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9, 8, 7,   5, 9, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 5, 8, 4,   5,10, 8,  10,11, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,
 5, 0, 4,   5,11, 0,   5,10,11,  11, 3, 0,                                            -1,-1,-1,  -1,
 0, 1, 9,   8, 4,10,   8,10,11,  10, 4, 5,                                            -1,-1,-1,  -1,
10,11, 4,  10, 4, 5,  11, 3, 4,   9, 4, 1,   3, 1, 4,                                            -1,
 2, 5, 1,   2, 8, 5,   2,11, 8,   4, 5, 8,                                            -1,-1,-1,  -1,
 0, 4,11,   0,11, 3,   4, 5,11,   2,11, 1,   5, 1,11,                                            -1,
 0, 2, 5,   0, 5, 9,   2,11, 5,   4, 5, 8,  11, 8, 5,                                            -1,
 9, 4, 5,   2,11, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 2, 5,10,   3, 5, 2,   3, 4, 5,   3, 8, 4,                                            -1,-1,-1,  -1,
 5,10, 2,   5, 2, 4,   4, 2, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3,10, 2,   3, 5,10,   3, 8, 5,   4, 5, 8,   0, 1, 9,                                            -1,
 5,10, 2,   5, 2, 4,   1, 9, 2,   9, 4, 2,                                            -1,-1,-1,  -1,
 8, 4, 5,   8, 5, 3,   3, 5, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 4, 5,   1, 0, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 8, 4, 5,   8, 5, 3,   9, 0, 5,   0, 3, 5,                                            -1,-1,-1,  -1,
 9, 4, 5,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4,11, 7,   4, 9,11,   9,10,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 8, 3,   4, 9, 7,   9,11, 7,   9,10,11,                                            -1,-1,-1,  -1,
 1,10,11,   1,11, 4,   1, 4, 0,   7, 4,11,                                            -1,-1,-1,  -1,
 3, 1, 4,   3, 4, 8,   1,10, 4,   7, 4,11,  10,11, 4,                                            -1,
 4,11, 7,   9,11, 4,   9, 2,11,   9, 1, 2,                                            -1,-1,-1,  -1,
 9, 7, 4,   9,11, 7,   9, 1,11,   2,11, 1,   0, 8, 3,                                            -1,
11, 7, 4,  11, 4, 2,   2, 4, 0,                                            -1,-1,-1,  -1,-1,-1,  -1,
11, 7, 4,  11, 4, 2,   8, 3, 4,   3, 2, 4,                                            -1,-1,-1,  -1,
 2, 9,10,   2, 7, 9,   2, 3, 7,   7, 4, 9,                                            -1,-1,-1,  -1,
 9,10, 7,   9, 7, 4,  10, 2, 7,   8, 7, 0,   2, 0, 7,                                            -1,
 3, 7,10,   3,10, 2,   7, 4,10,   1,10, 0,   4, 0,10,                                            -1,
 1,10, 2,   8, 7, 4,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 9, 1,   4, 1, 7,   7, 1, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,
 4, 9, 1,   4, 1, 7,   0, 8, 1,   8, 7, 1,                                            -1,-1,-1,  -1,
 4, 0, 3,   7, 4, 3,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 4, 8, 7,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 9,10, 8,  10,11, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 3, 0, 9,   3, 9,11,  11, 9,10,                                            -1,-1,-1,  -1,-1,-1,  -1,
 0, 1,10,   0,10, 8,   8,10,11,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3, 1,10,  11, 3,10,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 2,11,   1,11, 9,   9,11, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,
 3, 0, 9,   3, 9,11,   1, 2, 9,   2,11, 9,                                            -1,-1,-1,  -1,
 0, 2,11,   8, 0,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 3, 2,11,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 2, 3, 8,   2, 8,10,  10, 8, 9,                                            -1,-1,-1,  -1,-1,-1,  -1,
 9,10, 2,   0, 9, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 2, 3, 8,   2, 8,10,   0, 1, 8,   1,10, 8,                                            -1,-1,-1,  -1,
 1,10, 2,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 1, 3, 8,   9, 1, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 9, 1,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
 0, 3, 8,                                            -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,
                                         -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,  -1,-1,-1,   -1
};


const bool VERTS_FROM_MAP = true;


typedef struct {
    std::vector<REAL> &verts3;
    std::vector<int> &faces3;
    e3map_t &e3map;
    int& next_unique_vect_counter;
} result_state;


void MarchingCubes::flush_geometry_queue(std::ostream& cout, int& marching_cube_start,

  std::vector<REAL> &verts3, std::vector<int> &faces3, e3map_t &e3map, int& next_unique_vect_counter)
{

    for ( int vert_i = 0; vert_i < this->queue_counter; vert_i++ ) {

        int a = vert_i * 3;
        int b = a + 1;
        int c = a + 2;

        REAL x,y,z;
        x = this->positionQueue[ a ];
        y = this->positionQueue[ b ];
        z = this->positionQueue[ c ];
        if(!VERTS_FROM_MAP){
            verts3.push_back(x);
            verts3.push_back(y);
            verts3.push_back(z);
        }

        if(VERTS_FROM_MAP)
        {
            index3_t  e3_code = this->e3Queue[vert_i];
            std::pair<e3map_t::iterator, bool> e = e3map.emplace(e3_code, next_unique_vect_counter);
            const bool& novel = e.second;


            int overall_vert_index;
            if(novel)
                overall_vert_index = next_unique_vect_counter;
            else
                overall_vert_index = e.first->second;

            if(novel)
                next_unique_vect_counter++;

            if(novel){
                verts3.push_back(x);
                verts3.push_back(y);
                verts3.push_back(z);
            }
            else{
                assert(verts3.size()/3 > overall_vert_index);
                assert(overall_vert_index < next_unique_vect_counter);
            }

            //Loop invariant
            assert(verts3.size()/3 == next_unique_vect_counter);

            faces3.push_back(overall_vert_index);
            int old_overall_vert_index = vert_i + marching_cube_start*3;
        }



    }


    int nfaces = this->queue_counter / 3;

    for ( int face_i = 0; face_i < nfaces; face_i++ ) {

        int a = ( marching_cube_start + face_i ) * 3;
        int b = a + 1;
        int c = a + 2;


        if(!VERTS_FROM_MAP){
            faces3.push_back(a);
            faces3.push_back(b);
            faces3.push_back(c);
        }

    }

    marching_cube_start += nfaces;
    this->queue_counter = 0;

    if(REPORT_STATS){
    std::cout << "flush_geometry_queue(): " ;
    int mapctr = 0;
    for (auto& kv_pair: e3map){
        if(0)
            std::cout << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
    }
    std::cout << " e3Map: " << mapctr;
    std::cout << " Faces: " << faces3.size()/3;
    std::cout << " Verts: " << verts3.size()/3;
    std::cout << std::endl;
    }
}


class MarchingCubesMock {

public:
    MarchingCubesMock( dim_t resolution, bool enableUvs, bool enableColors ) {};
    ~MarchingCubesMock() {};

    void flush_geometry_queue(std::ostream& cout, int& marching_cube_start, std::vector<REAL> &verts3, std::vector<int> &faces3, e3map_t &e3map, int& next_unique_vect_counter)
        {};

    inline int polygonize_cube( REAL fx, REAL fy, REAL fz, index_t q, REAL isol, const callback_t& callback ) {return 0;};

    REAL isolation;

//shape:
    void addBall( REAL ballx, REAL bally, REAL ballz, REAL strength, REAL subtract ) {};
    void addPlaneX( REAL strength, REAL subtract ) {};
    void addPlaneZ( REAL strength, REAL subtract ) {};
    void addPlaneY( REAL strength, REAL subtract ) {};
    void seal_exterior(const REAL exterior_value = -1.) {};
//field
    void reset() {};

//geometry/threejs interface side.
    void render_geometry(const callback_t& renderCallback ) {};
    void sow() {};

// output. filled using sow()
    int resultqueue_faces_start = 0;
    std::vector<REAL> result_verts;
    std::vector<int> result_faces;
};



extern "C" {
    void build_geometry(int resolution, REAL time);
    int get_v_size();
    int get_f_size();
    void get_f(int*, int);
    void get_v(REAL*, int);
    void finish_geometry();
    void* get_f_ptr();
    void* get_v_ptr();
};


typedef struct {
    bool active = 0;
    MarchingCubes* mc = 0;
} state_t;

state_t _state;


void check_state() {
    if(!_state.active) std::cout << "Error: not active.";
}
void check_state_null() {
    if(_state.active)
        std::cout << "Error: should not be active.";
}

void build_geometry(int resolution, REAL time){

    check_state_null();


    bool enableUvs = true;
    bool enableColors = true;

    _state.mc = new MarchingCubes(resolution, enableUvs, enableColors);

    _state.mc -> isolation = 80.0/4;
      // before we had some amazing meatballs! merde a celui qui le lira!

    int min_z = - _state.mc->size2;
    int max_z = _state.mc->size2;
    int min_x = - _state.mc->size;
    int max_x = _state.mc->size;
    int min_y = - _state.mc->size;
    int max_y = _state.mc->size;

    boost::array<int, 2> grid_shape = {{ resolution*resolution*resolution , 3 }};
    boost::multi_array<REAL, 2> grid(grid_shape);

    boost::array<int, 1> implicit_function_shape = {{ resolution*resolution*resolution }};
    boost::multi_array<REAL, 1> implicit_function(implicit_function_shape);

    for (int z = min_z; z < max_z; z++ ) {
        for (int y = min_y; y < max_y; y++ ) {

            for (int x = min_x; x < max_x; x++ ) {
                cout << "Patate" << endl;
                grid[x + y*_state.mc->size + z*_state.mc->size2][0] = (REAL)x/(REAL)_state.mc->size;
                grid[x + y*_state.mc->size + z*_state.mc->size2][1] = (REAL)y/(REAL)_state.mc->size;
                grid[x + y*_state.mc->size + z*_state.mc->size2][2] = (REAL)z/(REAL)_state.mc->size;
            }
        }
    }

    cout << "Ligne 1301!" << endl;
    unit_sphere sphere(2.0);
    sphere.eval_implicit(grid, implicit_function);

    for (int z = min_z; z < max_z; z++ ) {
        for (int y = min_y; y < max_y; y++ ) {
            for (int x = min_x; x < max_x; x++ ) {
              _state.mc->field[x + y*_state.mc->size + z*_state.mc->size2] = implicit_function[x + y*_state.mc->size + z*_state.mc->size2];
            }
        }
    }
    cout << "Ligne 1312!" << endl;
    _state.mc->seal_exterior();

    const callback_t renderCallback;
    _state.mc->render_geometry(renderCallback);

    if(REPORT_STATS){
    int mapctr = 0;
    for (auto& kv_pair: _state.mc->result_e3map){
        if(0)
            std::cout << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
      }

    }

    if(VERBOSE){
        std::cout << resolution << " " << time << std::endl;
        std::cout << _state.mc << std::endl;
    }
    _state.active = true;

    check_state();
}
int get_f_size() {
    check_state();
    return _state.mc->result_faces.size()/3;
}
int get_v_size(){
    check_state();
    return _state.mc->result_verts.size()/3;
}
void get_v(REAL* v_out, int vcount){
    check_state();

    // Vertices
    int ctr = 0;
    for(std::vector<REAL>::iterator it=_state.mc->result_verts.begin(); it < _state.mc->result_verts.end(); it+=3 ){
        for(int di=0; di<3; di++){
            v_out[ctr] = *( it + di );
            ctr++;
        }
    }

    if(vcount*3 != ctr)  std::cout << "sizes dont match: " << (float)ctr/3. << " " << vcount << std::endl;
}

void get_f(int* f_out, int fcount){
    check_state();

    int ctr = 0;
    for(std::vector<int>::iterator it=_state.mc->result_faces.begin(); it < _state.mc->result_faces.end(); it+=3 ){
        for(int di=0; di<3; di++){
            f_out[ctr] = *( it + di );

            ctr++;
        }
    }
    if(fcount*3 != ctr)  std::cout << "sizes dont match: " << (float)ctr/3. << " " << fcount << std::endl;

};

void* get_v_ptr(){
    check_state();
    return (void*)(_state.mc->result_verts.data());
}

void* get_f_ptr(){
    check_state();
    return (void*)(_state.mc->result_faces.data());
}


void finish_geometry() {
    check_state();
    if(_state.mc == 0){
        std::cout << "Error: finish_geometry() before producing the shape()" << std::endl;
    }
    if(!_state.active){

    }
    else{
    }
    delete _state.mc;
    _state.active = false;
    _state.mc = 0;
};


#include "timer.hpp"

int main() {

    std::cout << "main();" << std::endl;
    return 0;
}
