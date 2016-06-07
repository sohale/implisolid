//AlteredQualia's Marching Cubes, C++ version, based on AQ's code.

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <cassert>

/**
 * Based on alteredq's version  https://github.com/mrdoob/three.js/blob/master/examples/js/MarchingCubes.js
 *
 * Port of greggman's ThreeD version of marching cubes to Three.js
 * http://webglsamples.googlecode.com/hg/blob/blob.html
 */


//typedef unsigned short int size_t;
typedef unsigned short int dim_t; //small integers for example the size of one side of the grid
typedef float REAL;
typedef unsigned long int index_t;


REAL lerp( REAL a, REAL b, REAL t )
{
    return a + ( b - a ) * t;
};


class MarchingCubes{
    bool enableUvs, enableColors;
    dim_t resolution;
    REAL isolation;
    index_t size, size2, size3;
    index_t  yd, zd;
    REAL halfsize;
    REAL delta;
    //boost::multi_array<REAL, 1> field;
    boost::multi_array<REAL, 1> field;
    boost::multi_array<REAL, 1> normal_cache;

    //parameters
    const static dim_t maxCount = 4096; // TODO: find the fastest size for this buffer

protected:
    //caching
    boost::multi_array<REAL, 1> vlist;
    boost::multi_array<REAL, 1> nlist;

protected:
    void init( dim_t resolution );
public:
    MarchingCubes( dim_t resolution, bool enableUvs, bool enableColors );
};

//static dim_t MarchingCubes::maxCount = 12;


void MarchingCubes::init( dim_t resolution )
{

        this->resolution = resolution;

        // parameters

        this->isolation = 80.0;

        // size of field, 32 is pushing it in Javascript :)

        this->size = resolution;
        this->size2 = this->size * this->size;
        this->size3 = this->size2 * this->size;
        this->halfsize = ((REAL)this->size) / 2.0;

        // deltas

        this->delta = 2.0 / (REAL)this->size;
        this->yd = this->size;
        this->zd = this->size2;


        //this->field = new Float32Array( this.size3 );
        //this->field = boost::array<REAL, 1>(this->size3); //does not work
        //this->field = std::array<REAL, this->size3>();
        //need a guarantee:
        assert(this->size3 < 10000000); // todo: get available heap. make it an exception.
        assert(this->size3 > 0);
        boost::array<int, 1> field_shape = {{ (int)this->size3, }};
        this->field = boost::multi_array<REAL, 1>( field_shape );
        //this->normal_cache = new Float32Array( this.size3 * 3 );
        boost::array<int, 1> normals_shape = {{ (int)this->size3 * 3, }};
        this->normal_cache = boost::multi_array<REAL, 1>( normals_shape );

        // temp buffers used in polygonize

        //this.vlist = new Float32Array( 12 * 3 );
        //this.nlist = new Float32Array( 12 * 3 );
        boost::array<int, 1> twelve3 = {{ 12 * 3, }};
        this->vlist = boost::multi_array<REAL, 1>( twelve3 );
        this->nlist = boost::multi_array<REAL, 1>( twelve3 );


        // immediate render mode simulator

        //this::maxCount = 4096; // TODO: find the fastest size for this buffer

        this->count = 0;
        /*

        this.hasPositions = false;
        this.hasNormals = false;
        this.hasColors = false;
        this.hasUvs = false;

        this.positionArray = new Float32Array( this.maxCount * 3 );
        this.normalArray   = new Float32Array( this.maxCount * 3 );

        if ( this.enableUvs ) {

            this.uvArray = new Float32Array( this.maxCount * 2 );

        }

        if ( this.enableColors ) {

            this.colorArray   = new Float32Array( this.maxCount * 3 );

        }
        */
    };

MarchingCubes::MarchingCubes( dim_t resolution, bool enableUvs=false, bool enableColors=false )
{

    //THREE.ImmediateRenderObject.call( this, material );

    this->enableUvs = enableUvs;
    this->enableColors = enableColors;

    // functions have to be object properties
    // prototype functions kill performance
    // (tested and it was 4x slower !!!)

//
}

#include "timer.hpp"

int main(){
    timer t;
    MarchingCubes mc( dim_t resolution, bool enableUvs, bool enableColors );
    t.stop();
    return 0;
}
