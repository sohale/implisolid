/* Copyright 2016  @sohale */

/**
 * AlteredQualia's Marching Cubes, C++ version, based on AQ's code based on Henrik Rydgård and @greggman.
 * https://github.com/WebGLSamples/WebGLSamples.github.io/blob/master/blob/marching_cubes.js
 *
 * Based on alteredq's version  https://github.com/mrdoob/three.js/blob/master/examples/js/MarchingCubes.js
 *
 * Port of greggman's ThreeD version of marching cubes to Three.js
 * http://webglsamples.googlecode.com/hg/blob/blob.html
 */

/*
Todo:
- WebWorker (incremental)
- Improve JS data exchange to Geometry
- Adding walls
- WebWorker wrappen in a Geometry
- Geometry in 3js r73
-
*/

#include <cassert>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"


// #include <math.h>

extern "C" {
    void produce_object_old2(float* verts, int *nv, int* faces, int *nf, float param);
    int main(int argc, char **argv);
}

const bool VERBOSE = false;
const bool REPORT_STATS = false;

// typedef unsigned short int size_t;
typedef uint16_t dim_t;  // small integers for example the size of one side of the grid
typedef float REAL;
// typedef unsigned long int index_t;

// boost::array will not work becasue the size of a boost:array has to be known in compile-time (static).
typedef  boost::multi_array<REAL, 1>  array1d;
// typedef array1d::index  array_index_t;
typedef boost::array<array1d::index, 1>  array_shape_t;
// #define array1d  boost::multi_array<REAL, 1>
typedef array1d::index  index_t;

typedef index_t index3_t;  // Range of the element type has to be large enough, larger than (size^3)*3.
typedef boost::multi_array<index3_t, 1>   array1d_e3;
typedef std::map<index3_t, int>  e3map_t;


struct callback_t { void call(void*) const { } callback_t(){} };

/*template<typename Index_Type=int>
boost::array<Index_Type, 1> make_shape_1d(Index_Type size)
{
    // Make a shape to be used in array initialisation
    // fixme: the type of size

    // ASSERT(size>=0);
    boost::array<Index_Type, 1> shape = {{ size, }};
    return shape;
}
*/

#include "mcc2_marching_cubes.hpp"
#include "tests/marching_cubes_mock.hpp"

void meta_balls(MarchingCubes& mc, int num_blobs, REAL time, REAL scale) {
    int numblobs = num_blobs;  // default: 4
    for (int ball_i = 0; ball_i < numblobs; ball_i++) {
        REAL D = 1;
        REAL ballx = sin(ball_i + 1.26 * time * (1.03 + 0.5*cos(0.21 * ball_i))) * 0.27 * D + 0.5;
        REAL bally = std::abs(cos(ball_i + 1.12 * time * cos(1.22 + 0.1424 * ball_i))) * 0.77 * D;  // dip into the floor
        REAL ballz = cos(ball_i + 1.32 * time * 0.1*sin((0.92 + 0.53 * ball_i))) * 0.27 * D + 0.5;
        REAL subtract = 12;
        REAL strength = 1.2 / ((sqrt(numblobs)- 1) / 4 + 1);
        mc.addBall(ballx, bally, ballz, strength, subtract, scale);
    }
}


/*********************************************************
    Old version, legacy
 *********************************************************/

void build_vf(
    // std::vector<REAL>& verts3,
    // std::vector<int>& faces3
    REAL scale
    ) {
    //  Includes allocations.

    dim_t resolution = 28;  // 28;
    bool enableUvs = true;
    bool enableColors = true;
    mp5_implicit::bounding_box box = {0, 3, 0, 3, 0, 3};

    MarchingCubes mc(resolution, box, enableUvs, enableColors);


    REAL time = 0.1;
    meta_balls(mc, 4, time, scale);
    mc.seal_exterior();

    /*
    int numblobs = 4;
    REAL subtract = (REAL)12.;
    REAL strength = (REAL)(1.2 / ( ( sqrt( numblobs ) - 1. ) / 4. + 1. ));
    mc.addBall(0.5, 0.5, 0.5, strength, subtract, 1);
    */
    // MarchingCubes& object = mc;
    // mc.addBall(0.5, 0.5, 0.5, strength, subtract, 1);

    // mc.flush_geometry_queue(std::cout, mc.resultqueue_faces_start, mc.result_normals, verts3, faces3);

    const callback_t renderCallback;
    mc.render_geometry(renderCallback);

    if (VERBOSE)
        std::cout << "MC:: v,f: " << mc.result_verts.size() << " " << mc.result_faces.size() << std::endl;

    // verts3.resize(0);
    // faces3.resize(0);
}



// Not used
void produce_object_old2(REAL* verts, int *nv, int* faces, int *nf, REAL time, REAL scale) {
    // not used

    dim_t resolution = 28;  // 28;
    bool enableUvs = true;
    bool enableColors = true;

    mp5_implicit::bounding_box box = {0, 3, 0, 3, 0, 3};

    if (VERBOSE)
        std::cout << "Leak-free (old version)" << std::endl;


    MarchingCubes mc(resolution, box, enableUvs, enableColors);
    // MarchingCubes* mc0 = new MarchingCubes(resolution, enableUvs, enableColors);
    // MarchingCubes &mc = *mc0;
    // MarchingCubesMock mc(resolution, enableUvs, enableColors);

    // REAL time = 0.1 ;
    meta_balls(mc, 4, time, scale);
    mc.seal_exterior();

    /*
    int numblobs = 4;
    REAL subtract = (REAL)12.;
    REAL strength = (REAL)(1.2 / ( ( sqrt( numblobs ) - 1. ) / 4. + 1. ));
    mc.addBall(0.5, 0.5, 0.5, strength, subtract, scale);
    mc.addBall(0, 0, 0.5, strength, subtract, scale);
    */

    // todo: init-receiver side

    const callback_t renderCallback;
    mc.render_geometry(renderCallback);

    std::cout << "map2" << std::endl;

    // mc.result_faces.resize(100);

    if (VERBOSE)
        std::cout << "MC:: v,f: " << mc.result_verts.size() << " " << mc.result_faces.size() << std::endl;

    *nv = mc.result_verts.size()/3;
    *nf = mc.result_faces.size()/3;

    // Vertices
    int ctr = 0;
    for ( std::vector<REAL>::iterator it=mc.result_verts.begin(); it < mc.result_verts.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            verts[ctr] = *(it + di);
            ctr++;
        }
      }

    // Faces
    ctr = 0;
    for ( std::vector<int>::iterator it=mc.result_faces.begin(); it < mc.result_faces.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            faces[ctr] = *(it + di);
            ctr++;
        }
      }
}


/*
#include <emscripten/bind.h>
using namespace emscripten;
EMSCRIPTEN_BINDINGS(my_module) {
    // produce_object(REAL* verts, int *nv, int* faces, int *nf, REAL time);
    function("produce_object_bind", &produce_object);
}
*/
// #include "mcc1-glue.cpp"


/*void produce_v(){
}*/



/*********************************************************
    public interface for JavaScript
 *********************************************************/


extern "C" {
    void build_geometry(char* shape_parameters_json, char* mc_parameters_json);
    int get_v_size();
    int get_f_size();
    void get_f(int*, int);
    void get_v(REAL*, int);
    void finish_geometry();
    void* get_f_ptr();
    void* get_v_ptr();
    void implicit_value(REAL x, REAL y, REAL z);
    void implicit_grad_value(REAL x, REAL y, REAL z, REAL * x_out, REAL * y_out, REAL * z_out);
    // also: queue, etc.
    // bad: one instance only.
    //     Solution 1:  MarchingCubes* build_geometry();
    //     Solution 2: ids (for workers! ; a statically determined number of them (slots/workers/buckets).).
};


typedef struct {
    bool active = 0;
    MarchingCubes* mc = 0;
} state_t;

state_t _state;

// _state.active = false;
// _state.mc = 0;

bool check_state() {
    if (!_state.active) {
        std::cout << "Error: There are no allocated geometry resources to deallocate.";
        return false;
    }
    return true;
}
bool check_state_null() {
    if (_state.active) {
        std::cout << "Error: There are non-deallocated geometry resources. Call finit_geometry() first. Not doing anything.";
        return false;
    }
    return true;
}


namespace mp5_implicit{
struct mc_settings {
    mp5_implicit::bounding_box box;
    int resolution;
};
}

mp5_implicit::mc_settings parse_mc_properties_json(char* mc_parameters_json) {
    std::stringstream mc_json_stream;
    mc_json_stream << mc_parameters_json;

    namespace pt = boost::property_tree;
    pt::ptree mcparams_dict;


    // TODO(charles): find an alternativ to catch exceptions pt::json_parser::json_parser_error pt::ptree_bad_path
    // try{
    pt::read_json(mc_json_stream, mcparams_dict);

    REAL xmin = mcparams_dict.get<REAL>("box.xmin", NaN);
    REAL xmax = mcparams_dict.get<REAL>("box.xmax", NaN);
    REAL ymin = mcparams_dict.get<REAL>("box.ymin", NaN);
    REAL ymax = mcparams_dict.get<REAL>("box.ymax", NaN);
    REAL zmin = mcparams_dict.get<REAL>("box.zmin", NaN);
    REAL zmax = mcparams_dict.get<REAL>("box.zmax", NaN);

    int resolution = mcparams_dict.get<int>("resolution", -1);

    if ( isNaN(xmin) || isNaN(xmax) || isNaN(ymin) || isNaN(ymax) || isNaN(zmin) || isNaN(zmax) || resolution <= 2 ) {
        std::cout << "Error: missing or incorrect values in mc_parameters_json"<< std::endl;
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        resolution = 28;
    }
    // std::cout << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << " " << resolution << " " << std::endl;


    /*}catch(pt::json_parser::json_parser_error parse_exception){
        std::cout << "parse_exception"<< std::endl ;

    }catch(pt::ptree_bad_data bad_data_exception){
        std::cout << "bad_data_exception"<< std::endl ;

    }catch(pt::ptree_bad_path bad_path_exception){1
        std::cout << "bad_path_exception" << std::endl ;

    }catch(...){
        std::cout << "other_exception" << std::endl ;
    }*/


    mp5_implicit::mc_settings  mc_settings_from_json;  // settings
    mp5_implicit::bounding_box box = {xmin, xmax, ymin, ymax, zmin, zmax};  // {15,20,15,20,15,20};
    mc_settings_from_json.box = box;
    mc_settings_from_json.resolution = resolution;

    return mc_settings_from_json;
}


#include "../js_iteration_2/object_factory.hpp"

// void build_geometry(int resolution, char* mc_parameters_json, char* obj_name, REAL time){
void build_geometry(char* shape_parameters_json, char* mc_parameters_json) {
    if (!check_state_null())
        return;
    // std::cout << "In build_geometry obj_name : " << obj_name << std::endl;
    // std::cout << "Mc_params : " << mc_parameters_json << endl;
    // std::cout << "shape_json : " << shape_parameters_json << endl;

    bool use_metaball;
    std::string shape_parameters_json_str = std::string(shape_parameters_json);
    implicit_function* object = object_factory(shape_parameters_json_str , use_metaball);


    mp5_implicit::mc_settings  mc_settings_from_json = parse_mc_properties_json(mc_parameters_json);

    // std::cout << "Leak-free : new" << std::endl;

    // dim_t resolution = 28;
    bool enableUvs = true;
    bool enableColors = true;

    // MarchingCubes mc(resolution, enableUvs, enableColors);
    _state.mc = new MarchingCubes(mc_settings_from_json.resolution, mc_settings_from_json.box, enableUvs, enableColors);
    // std::cout << "constructor called. " << _state.mc << std::endl;

    _state.mc -> isolation = 80.0/4*0;



    // ****************************
    // Does thismake things slower ?
    boost::multi_array<REAL, 2>  mcgrid_vectorized = _state.mc -> prepare_grid();  // 10.0
    _state.mc -> eval_shape(*object, mcgrid_vectorized);

     delete object;
     object = NULL;
     if (use_metaball) {
         REAL metaball_time = 0;
         meta_balls(*_state.mc, 4, metaball_time, 1.0);
     }

    _state.mc->seal_exterior(-10000000.0);

    // std::cout << "balls added." << std::endl;

    const callback_t renderCallback;
    _state.mc->render_geometry(renderCallback);
    // std::cout << "MC executed" << std::endl;

    // std::cout << "map4" << std::endl;

    if (REPORT_STATS) {
    int mapctr = 0;
    for (auto& kv_pair : _state.mc->result_e3map) {
        if (0)
            std::cout << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
    }
    /*
    std::cout << "build_geometry(): ";
    std::cout << " e3Map: " << mapctr;
    std::cout << " Faces: " << _state.mc->result_faces.size()/3;
    std::cout << " Verts: " << _state.mc->result_verts.size()/3;
    std::cout << std::endl;
    */
    }


    _state.active = true;

    check_state();
    // std::cout << "MC:: v,f: " << _state.mc->result_verts.size() << " " << _state.mc->result_faces.size() << std::endl;
}
int get_f_size() {
    if (!check_state()) return -1;
    return _state.mc->result_faces.size()/3;
}
int get_v_size() {
    if (!check_state()) return -1;
    return _state.mc->result_verts.size()/3;
}
void get_v(REAL* v_out, int vcount) {
    if (!check_state())
        return;

    // int nf = get_f_size();
    // Vertices
    int ctr = 0;
    for ( std::vector<REAL>::iterator it=_state.mc->result_verts.begin(); it < _state.mc->result_verts.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            v_out[ctr] = *(it + di);
            // if(ctr<3*3*3)
            //    std::cout << v_out[ctr] << " ";
            ctr++;
        }
    }
    // std::cout << std::endl;
    // assert nf*3 == ctr;
    if (vcount*3 != ctr)  std::cout << "sizes dont match: " << static_cast<float>(ctr)/3. << " " << vcount << std::endl;
}



void get_f(int* f_out, int fcount) {
    if (!check_state())
        return;
    // int nf = get_f_size();
    int ctr = 0;
    for ( std::vector<int>::iterator it=_state.mc->result_faces.begin(); it < _state.mc->result_faces.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            f_out[ctr] = *(it + di);
            // if(ctr<3*3*3)
            //    std::cout << f_out[ctr] << " ";
            ctr++;
        }
    }
    if (fcount*3 != ctr)  std::cout << "sizes dont match: " << static_cast<float>(ctr)/3. << " " << fcount << std::endl;
    // std::cout << std::endl;
}


/* Data is already there, so why copy it? Also, keep it until next round. */
void* get_v_ptr() {
    if (!check_state())
        return NULL;
    return reinterpret_cast<void*>(_state.mc->result_verts.data());
}

void* get_f_ptr() {
    if (!check_state())
        return NULL;
    return reinterpret_cast<void*>(_state.mc->result_faces.data());
}


// int get_v_size(){};
// int get_f_size(){};
// void get_f(int*){};
// void get_v(REAL*){};
// void* get_f_ptr();
// void* get_v_ptr();
// void finish_geometry();

//  Can cause an exception (but not when null).
void finish_geometry() {
    if (!check_state())
        return;
    if (_state.mc == 0) {
        std::cerr << "Error: finish_geometry() before producing the shape()" << std::endl;
    }
    if (!_state.active) {
        // std::cout << "Cannot finish_geometry() while still active." << std::endl;
    } else {
        // std::cout << "_state.active " << _state.active << "  _state.mc " << _state.mc << std::endl;
    }
    // Dos not cause an exception if null. But it causes exception.
    delete _state.mc;
    _state.active = false;
    _state.mc = 0;
}


#include "timer.hpp"

/*int main() {

    timer t;
    t.stop();
    // MarchingCubes mc( dim_t resolution, bool enableUvs, bool enableColors );
    dim_t resolution = 28;  // 28;
    bool enableUvs = true;
    bool enableColors = true;
    MarchingCubes mc(resolution, enableUvs, enableColors);
    t.stop();

    int numblobs = 4;
    REAL subtract = (REAL)12.;
    REAL strength = (REAL)(1.2 / ( ( sqrt( numblobs ) - 1. ) / 4. + 1. ));

    REAL scale = 1;
    mc.addBall(0.5, 0.5, 0.5, strength, subtract, scale);

    const callback_t renderCallback;
    mc.render_geometry(renderCallback);
    t.stop();


    std::vector<REAL> verts3;
    std::vector<int> faces3;
    MarchingCubes& object = mc;

    // int normals_start = 0;
    mc.flush_geometry_queue(std::cout, mc.resultqueue_faces_start, mc.result_normals, verts3, faces3);

    t.stop();

    cout << resolution << endl;

    cout << endl;


    cout << "verts, faces: ";
    cout << mc.result_verts.size();
    cout << " ";
    cout << mc.result_faces.size();
    cout << endl;

    t.stop();

    // build_vf( verts3, faces3 );  // 21.3 msec using O3
    build_vf(  );  // 26 msec.


    t.stop();

    std::cout << "main();" << std::endl;
    return 0;
}
*/
