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


#include "../js_iteration_2/object_collector.hpp"

#include "../js_iteration_2/pointset_set.hpp"

#include "../js_iteration_2/vertex_resampling.hpp"
#include "../js_iteration_2/apply_v_s_to_mc_buffers.hpp"
#include "../js_iteration_1/centroids_projection.cpp"

#include "../js_iteration_2/faces_verts_algorithms.hpp"

#include "../js_iteration_2/polygoniser_settings.hpp"

#include "timer.hpp"

// #include <math.h>

extern "C" {
    void produce_object_old2(float* verts, int *nv, int* faces, int *nf, float param);
    int main(int argc, char **argv);
}

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

    // mc.flush_geometry_queue(std::clog, mc.resultqueue_faces_start, mc.result_normals, verts3, faces3);

    const callback_t renderCallback;
    mc.render_geometry(renderCallback);

    if (VERBOSE)
        std::clog << "MC:: v,f: " << mc.result_verts.size() << " " << mc.result_faces.size() << std::endl;

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
        std::clog << "Leak-free (old version)" << std::endl;


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

    std::clog << "map2" << std::endl;

    // mc.result_faces.resize(100);

    if (VERBOSE)
        std::clog << "MC:: v,f: " << mc.result_verts.size() << " " << mc.result_faces.size() << std::endl;

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
    void build_geometry(const char* shape_parameters_json, const char* mc_parameters_json);
    int get_v_size();
    int get_f_size();
    void get_f(int*, int);
    void get_v(REAL*, int);
    void finish_geometry();
    void* get_f_ptr();
    void* get_v_ptr();

    // void implicit_value(REAL x, REAL y, REAL z);
    // void implicit_grad_value(REAL x, REAL y, REAL z, REAL * x_out, REAL * y_out, REAL * z_out);

    /********************************************************
        New API for direct evaluation of implicit functions
    *********************************************************/
    //global variables: mp5_implicit::implicit_function* last_object, last_x, v_results, g_results,
    void set_object(const char* shape_parameters_json, bool ignore_root_matrix);  // call object_factory. Sets a global variable   last_object
    void unset_object();

    bool set_x(void* verts, int n);   // sets last_x, a multi_array.
    void unset_x();

    void calculate_implicit_values();  // multi_array or vector that has a .data().  // calculate v_results from  last_object.eval_implicit()
    void* get_values_ptr();  // return v_results.data()
    int get_values_size();   // return v_results

    void calculate_implicit_gradients(bool normalize_and_invert);  // calculate g_results from  last_object.eval_gradient()
    void* get_gradients_ptr(); // return g_results.data()
    int get_gradients_size(); // return g_results


    // query point_set. For example, centroids.
    void* get_pointset_ptr(char* id);
    int get_pointset_size(char* id);


    void about();  // shows build information

    //void calculate_gradients()

    // also: queue, etc.
    // bad: one instance only.
    //     Solution 1:  MarchingCubes* build_geometry();
    //     Solution 2: ids (for workers! ; a statically determined number of them (slots/workers/buckets).).
};



// typedef    std::tuple<std::vector<REAL>, std::vector<vectorized_vect::index>> verts_faces_t;
// typedef    std::tuple<vectorized_vect, std::vector<vectorized_vect::index>> verts_faces_t;

typedef   std::tuple< vectorized_vect, vectorized_faces > vertsfaces_type;

class polygoniser {
    mp5_implicit::implicit_function const *  object;

    vertsfaces_type vertsfaces;  // not good

    /*
    std::vector<REAL> verts;
    std::vector<int> faces;
    const verts_t & centroids,
    */

public:
    polygoniser(mp5_implicit::implicit_function const * ifunc_object)
        :   object(ifunc_object)  // copy constructor
    {
        ;
    }
};


class state_t {
public:
    bool active = 0;
    MarchingCubes* mc = 0;
    polygoniser pgonizer();

public:
    bool check_state() {
        if (!this->active) {
            std::clog << "Error: There are no allocated geometry resources to deallocate.";
            return false;
        }
        return true;
    }
    bool check_state_null() {
        if (this->active) {
            std::clog << "Error: There are non-deallocated geometry resources. Call finit_geometry() first. Not doing anything.";
            return false;
        }
        return true;
    }

};

state_t _state;

// _state.active = false;
// _state.mc = 0;



//polygoniser_settings

// moved up
// std::map< std::string, vectorized_vect > point_set_set;

void* get_pointset_ptr(char* id) {
    void* p = point_set_set[std::string(id)].origin();
    return p;
}
int get_pointset_size(char* id) {
    clog << "int get_pointset_size(char* id) {return point_set_set[std::string(id)].shape()[0];} = "
        << point_set_set[std::string(id)].shape()[0]
        << std::endl;
    return point_set_set[std::string(id)].shape()[0];
}


#include "../js_iteration_2/object_factory.hpp"

// The only usage of marching cubes

std::pair< std::vector<REAL>, std::vector<int>>  mc_start (mp5_implicit::implicit_function* object, dim_t resolution_, const mp5_implicit::bounding_box & box_, const bool use_metaball) {

    bool enableUvs = true;
    bool enableColors = true;

    _state.mc = new MarchingCubes(resolution_, box_, enableUvs, enableColors);

    _state.mc -> isolation = 80.0/4*0;



    // ****************************
    // Does thismake things slower ?
    vectorized_vect   mcgrid_vectorized = _state.mc -> prepare_grid();  // 10.0
    _state.mc -> eval_shape(*object, mcgrid_vectorized);

     if (use_metaball) {
         REAL metaball_time = 0;
         meta_balls(*_state.mc, 4, metaball_time, 1.0);
     }

    _state.mc->seal_exterior(-10000000.0);

    const callback_t renderCallback;
    _state.mc->render_geometry(renderCallback);
    // std::clog << "MC executed" << std::endl;

    // std::clog << "map4" << std::endl;

    if (REPORT_STATS) {
    int mapctr = 0;
    for (auto& kv_pair : _state.mc->result_e3map) {
        if (0)
            std::clog << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
    }
    /*
    std::clog << "build_geometry(): ";
    std::clog << " e3Map: " << mapctr;
    std::clog << " Faces: " << _state.mc->result_faces.size()/3;
    std::clog << " Verts: " << _state.mc->result_verts.size()/3;
    std::clog << std::endl;
    */
    }

    return std::make_pair(_state.mc -> result_verts, _state.mc->result_faces);
}


std::vector<REAL> _state_mc_result_verts;
std::vector<int> _state_mc_result_faces;



// void build_geometry(int resolution, char* mc_parameters_json, char* obj_name, REAL time){
void build_geometry(const char* shape_parameters_json, const char* mc_parameters_json) {
    if (!_state.check_state_null()) {
        clog << "build_geometry() called in a bad state.";
        return;
    }
    // std::clog << "In build_geometry obj_name : " << obj_name << std::endl;
    // std::clog << "Mc_params : " << mc_parameters_json << endl;
    // std::clog << "shape_json : " << shape_parameters_json << endl;

    mp5_implicit::mc_settings  mc_settings_from_json = parse_mc_properties_json(mc_parameters_json);

    std::string shape_parameters_json_str = std::string(shape_parameters_json);
    bool ignore_root_matrix = mc_settings_from_json.ignore_root_matrix;

    //unique_pointer<mp5_implicit::implicit_function> object = ...;

    bool use_metaball;  // output
    mp5_implicit::implicit_function* object = object_factory(shape_parameters_json_str , use_metaball, ignore_root_matrix);

    // std::clog << "Leak-free : new" << std::endl;

    // dim_t resolution = 28;

    auto vertsfaces_pair = mc_start(object, mc_settings_from_json.resolution, mc_settings_from_json.box, use_metaball);
    // std::vector<REAL>, std::vector<int>

    // auto  _state_mc_result_verts = _state.mc -> result_verts;
    // auto  _state_mc_result_faces = _state.mc->result_faces;
    _state_mc_result_verts = std::move(vertsfaces_pair.first);
    _state_mc_result_faces = std::move(vertsfaces_pair.second);


    /*
    int vresamp_iters = 0; //3;
    bool apply_projection = false;
    float c = 1.;
    */

    /*
    // For steps I, II
    mc_settings_from_json.projection.enabled = true;
    mc_settings_from_json.vresampl.c = 1.0;
    mc_settings_from_json.vresampl.iters = 1;
    mc_settings_from_json.qem.enabled = true;
    */



    const bool DISABLE_POSTPROCESSING = false;    // DISABLE ALL MESH POST-PROCESSING (mesh optimisation)
    if (!DISABLE_POSTPROCESSING) {
    int vresamp_iters  =  mc_settings_from_json.vresampl.iters; //10; //3;
    bool apply_projection = true;
    // float c = 1.;
    REAL c =  mc_settings_from_json.vresampl.c;  // 1.0;

    if (VERBOSE) {
        clog << "vresampl.c: " << c << std::endl;
        clog << "vresamp_iters: " << vresamp_iters << std::endl;
        clog << ".projection.enabled: " << mc_settings_from_json.projection.enabled << std::endl;
        clog << ".qem.enabled: " << mc_settings_from_json.qem.enabled << std::endl;
    }


    timer timr;
    timr.report_and_continue("timer started.");

    for (int i=0; i < vresamp_iters; i++) {
        timer t1;
        // result_verts is modified
        apply_vertex_resampling_to_MC_buffers__VMS(object, c, _state_mc_result_verts, _state_mc_result_faces, false );
        t1.stop("vertex resampling");  // 400 -> 200 -> 52 msec  (40--70)
        timr.report_and_continue("vertex resampling");
    }

    if (apply_projection) {
        /*
        if (STORE_POINTSETS)
        {
            verts_t v = convert_vectorverts_to_vectorized_vect( _state_mc_result_verts);
            STORE_POINTSET("pre_p_centroid", v);
        }
        if (STORE_POINTSETS) {
            STORE_POINTSET("pre_p_verts", v);
        }
        */

        if (mc_settings_from_json.projection.enabled) {
            std::clog << "centroids_projection:" << std::endl;
            // Dn't send mc_settings_from_json as an argument
            centroids_projection(object, _state.mc->result_verts, _state_mc_result_faces, mc_settings_from_json.qem.enabled);
            timr.report_and_continue("centroids_projection");
        } else {
            // std::clog << "centroids_projection (& qem) skipped because you asked for it." << std::endl;
        }

        /*
        if (STORE_POINTSETS) {
            verts_t v = convert_vectorverts_to_vectorized_vect( _state_mc_result_verts);
            STORE_POINTSET("post_p_centroids", v);
        }
        */
    }
    }

    /*
    // int REPEATS = 2;
    //int REPEATS_VR = 3;
    int REPEATS = 1;
    int REPEATS_VR = 1;
    for (int repeats = 0; repeats < REPEATS ; ++repeats) {
        float c = 2000.;

        for (int i=0; i < REPEATS_VR; i++){
            apply_vertex_resampling_to_MC_buffers__VMS(object, c,  //    *(_state.mc));
                _state.mc->result_verts, _state_mc_result_faces, ??);
        }

        centroids_projection(object, _state.mc->result_verts, _state_mc_result_faces);
    }
    */

    //delete object;
    object = NULL;
    gc_objects();

    _state.active = true;

    _state.check_state();
    // std::clog << "MC:: v,f: " << _state.mc->result_verts.size() << " " << _state_mc_result_faces.size() << std::endl;
}
int get_f_size() {
    if (!_state.check_state()) return -1;
    return _state_mc_result_faces.size()/3;
}
int get_v_size() {
    if (!_state.check_state()) return -1;
    return _state.mc->result_verts.size()/3;
}
void get_v(REAL* v_out, int vcount) {
    if (!_state.check_state())
        return;

    // int nf = get_f_size();
    // Vertices
    int ctr = 0;
    for ( std::vector<REAL>::iterator it=_state.mc->result_verts.begin(); it < _state.mc->result_verts.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            v_out[ctr] = *(it + di);
            // if(ctr<3*3*3)
            //    std::clog << v_out[ctr] << " ";
            ctr++;
        }
    }
    // std::clog << std::endl;
    // assert nf*3 == ctr;
    if (vcount*3 != ctr)  std::clog << "sizes dont match: " << static_cast<float>(ctr)/3. << " " << vcount << std::endl;
}



void get_f(int* f_out, int fcount) {
    if (!_state.check_state())
        return;
    // int nf = get_f_size();
    int ctr = 0;
    for ( std::vector<int>::iterator it=_state_mc_result_faces.begin(); it < _state_mc_result_faces.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            f_out[ctr] = *(it + di);
            // if(ctr<3*3*3)
            //    std::clog << f_out[ctr] << " ";
            ctr++;
        }
    }
    if (fcount*3 != ctr)  std::clog << "sizes dont match: " << static_cast<float>(ctr)/3. << " " << fcount << std::endl;
    // std::clog << std::endl;
}


/* Data is already there, so why copy it? Also, keep it until next round. */
void* get_v_ptr() {
    if (!_state.check_state())
        return NULL;
    return reinterpret_cast<void*>(_state.mc->result_verts.data());
}

void* get_f_ptr() {
    if (!_state.check_state())
        return NULL;
    return reinterpret_cast<void*>(_state_mc_result_faces.data());
}


// int get_v_size(){};
// int get_f_size(){};
// void get_f(int*){};
// void get_v(REAL*){};
// void* get_f_ptr();
// void* get_v_ptr();
// void finish_geometry();

//  Can cause an exception (but not when null).
/*
    Leave behind and forget the current geometry. Ready for a new one.
*/
void finish_geometry() {
    if (!_state.check_state())
        return;
    if (_state.mc == 0) {
        std::cerr << "Error: finish_geometry() before producing the shape()" << std::endl;
    }
    if (!_state.active) {
        // std::clog << "Cannot finish_geometry() while still active." << std::endl;
    } else {
        // std::clog << "_state.active " << _state.active << "  _state.mc " << _state.mc << std::endl;
    }
    // Dos not cause an exception if null. But it causes exception.
    delete _state.mc;
    _state.active = false;
    _state.mc = 0;
}

// information
void about() {
    std::clog << "Build 1" << std::endl;
    std::clog << __DATE__ << " " << __TIME__ << std::endl;
}











#include "../js_iteration_2/basic_data_structures.hpp"
#include "../js_iteration_2/basic_functions.hpp"
#include "centroids_projection.cpp"



// ************************************************************
// EXPERIMENTAL
/*
typedef  char  iobjid_type;

// solids, solid_set, model, solid workers, etc
class iobject_set {

private:
    // maximum number of objects
    std::vector<mp5_implicit::implicit_function> current_objects{};
public:
    mp5_implicit::implicit_function* get_current_object_ptr(const iobjid_type id) {
        std::clog << "Dont use this yet." << std::endl;
        return NULL;
    }

public:
    iobjid_type set_object(const char* shape_parameters_json, bool ignore_root_matrix) {
        retur -1;
    }

    void set_matrix(const iobjid_type id, void* m12) {
    }

    //void set_eye(iobjid_type id) { }

    bool unset_object(const iobjid_type id) {
        return 0; // failed
    }
};
*/

/*
The plan is to keep the flow of polygonising of the objects in an object.
This object is not the SolidObject3D. It is just the polygonisation flow.
A polygonisation flow can be used using a builder class.
It has a reference to an object.
Each polygnisation flow is a task in a polygonisation flow and they are managed by a scheduler.
The scheduler can be implemented on JavaScript side (clients, not web worker. For now).

*/

typedef boost::array<vectorized_vect::index, 2>  shape_t;

// Query with the implicit object
// A service that can evaluate X points.
// Designed for querying with javascript.
// QueryableObject
// ifunction_service  pointset_service  iobj_service
// ifunction_service
// A service.
// An object for the implicit function, but not as an object necessarily.
// Used for evaluation of f() and gradient values.
// ifunction_service
// But stores the results too (cahches input and output).
// Good for asynchronous evaluation.
// For example: web-workers
// ifunction_point_service

// workers
struct {
    mp5_implicit::implicit_function* current_object = NULL;

    vectorized_vect* current_x = NULL;
    vectorized_vect* current_grad = NULL;
    vectorized_scalar* current_f = NULL;
} ifunction_service;


/*
Sets the object ptr of the ifunction_service from a model
*/
//void set_object_from_model(iobjid_type id) {}

void set_matrix(void* m12) {
    // keeps the same object. Just updates the outer matrix.
    // ifunction_service.current_object->setMatrix(m12);
    std::clog << "Error: Not implemented." << std::endl;
}




// *************************************************************************



void set_object(const char* shape_parameters_json, bool ignore_root_matrix) {
    if(ifunction_service.current_object != NULL){
        std::clog << "Error: You cannot unset() the object before a set_object(json)." << std::endl;
        return;
    }

    //std::clog << "before: current_object " << ifunction_service.current_object << std::endl;

    std::string str = std::string(shape_parameters_json);
    bool dummy;
    ifunction_service.current_object = object_factory(str , dummy, ignore_root_matrix);

    //std::clog << "after: current_object " << ifunction_service.current_object << std::endl;
}
// unset_object(int id)
void unset_object() {
    if(ifunction_service.current_object == NULL){
        std::clog << "Error: You cannot unset() the object before a set_object(json)." << std::endl;
        return;
    }

    //delete ifunction_service.current_object;
    //ifunction_service.current_object = NULL;
    gc_objects();
    ifunction_service.current_object = NULL;
}

bool set_x(void* verts, int n) {
    if(ifunction_service.current_x != NULL){
        std::clog << "Error: You set() before unset()ing the previous set()." << std::endl;
        return false;
    }
    if( n < 0 || n >= 10000) {
        std::clog << "Error: n is outside [0, 10000]." << std::endl;
        return false;
    }

    ifunction_service.current_x = new vectorized_vect( shape_t{n, 3} );
    REAL* real_verts = reinterpret_cast<REAL*>(verts);
    for(int i = 0; i < n; i++) {
        (*ifunction_service.current_x)[i][0] = real_verts[i*3 + 0];
        (*ifunction_service.current_x)[i][1] = real_verts[i*3 + 1];
        (*ifunction_service.current_x)[i][2] = real_verts[i*3 + 2];
        /*
        if(i < 10) {
            std::clog
                << (*(ifunction_service.current_x))[i][0] << " "
                << (*(ifunction_service.current_x))[i][1] << " "
                << (*(ifunction_service.current_x))[i][2] << " "
                << std::endl;
        }
        */
    }
    //std::clog << std::endl;
    ifunction_service.current_f = new vectorized_scalar( shape_t{n}  );  // n x 0 (?)
    //std::clog << "warning: size may be n x 0:  " << ifunction_service.current_f->shape()[0] << "x" << ifunction_service.current_f->shape()[1] << std::endl;
    ifunction_service.current_grad = new vectorized_vect( shape_t{n, 3}  );
    #if ASSERT_USED
    //ifunction_service.current_grad  <- some init value
    #endif
    return true;
}
void unset_x() {
    if(ifunction_service.current_x == NULL){
        std::clog << "Error: You cannot unset() before a set()." << std::endl;
        return;
    }

    delete ifunction_service.current_x;
    delete ifunction_service.current_f;
    delete ifunction_service.current_grad;
    ifunction_service.current_x = NULL;
    ifunction_service.current_f = NULL;
    ifunction_service.current_grad = NULL;
}


void calculate_implicit_values() {
    if(ifunction_service.current_x == NULL || ifunction_service.current_f == NULL || ifunction_service.current_object == NULL) {
        std::clog << "Error: You need to set_x() and set_object() first." << std::endl;
        return;
    }

    ifunction_service.current_object -> eval_implicit(*(ifunction_service.current_x), ifunction_service.current_f);
}
void* get_values_ptr() {
    return ifunction_service.current_f->data();
}
int get_values_size() {
    return ifunction_service.current_f->shape()[0];
}

void calculate_implicit_gradients(bool normalize_and_invert) {
    if(ifunction_service.current_x == NULL || ifunction_service.current_grad == NULL || ifunction_service.current_object == NULL) {
        std::clog << "Error: You need to set_x() and set_object() first." << std::endl;
        return;
    }

    int problems = 0;

    clog << "size consistency" ;
    if (ifunction_service.current_x->shape()[0] != ifunction_service.current_grad->shape()[0]) {
        clog << " " << ifunction_service.current_x->shape()[0] << "  !=  " << ifunction_service.current_grad->shape()[0];
            //<< std::endl;
    }
    clog << std::endl;
    assert(ifunction_service.current_x->shape()[0] == ifunction_service.current_grad->shape()[0]);


    ifunction_service.current_object -> eval_gradient(*(ifunction_service.current_x), ifunction_service.current_grad);
    if(normalize_and_invert) {
        for(auto it = ifunction_service.current_grad->begin(); it < ifunction_service.current_grad->end(); it++) {
            REAL x = (*it)[0];
            REAL y = (*it)[1];
            REAL z = (*it)[2];
            REAL norm = std::sqrt(x*x + y*y + z*z);

            REAL norm_factor;
            if (norm > 0.0001) {
                norm_factor = -1.0 / norm;
            } else {
                norm_factor = -42.0; // how to avoid look black
                problems++;
            }

            (*it)[0] = x * norm_factor;
            (*it)[1] = y * norm_factor;
            (*it)[2] = z * norm_factor;
        }
    }

    if (problems > 0) {
        cerr << " problems " << problems << std::endl;
    }

    /*
    std::clog << "calculated grad: "
        << (*(ifunction_service.current_grad))[0][0] << " "
        << (*(ifunction_service.current_grad))[0][1] << " "
        << (*(ifunction_service.current_grad))[0][2] << " "
        << "  calculated from x ="
        << (*(ifunction_service.current_x))[0][0] << " "
        << (*(ifunction_service.current_x))[0][1] << " "
        << (*(ifunction_service.current_x))[0][2] << " "
        << std::endl;
    */

}
void* get_gradients_ptr() {
    if(ifunction_service.current_x == NULL || ifunction_service.current_f == NULL || ifunction_service.current_object == NULL) {
        std::clog << "Error: You need to set_x() and set_object() first." << std::endl;
        return NULL;
    }
    /*
    std::clog << "current_grad: "
        << (*(ifunction_service.current_grad))[0][0] << " "
        << (*(ifunction_service.current_grad))[0][1] << " "
        << (*(ifunction_service.current_grad))[0][2] << " "
        << std::endl;
    */
    return ifunction_service.current_grad->data();
}
int get_gradients_size() {
    if(ifunction_service.current_x == NULL || ifunction_service.current_f == NULL || ifunction_service.current_object == NULL) {
        std::clog << "Error: You need to set_x() and set_object() first." << std::endl;
        return 0;
    }

    return ifunction_service.current_grad->shape()[0] * 3;
}




int main(int argc, char **argv) {
    std::clog.tie (&cerr);
    clog << "clog" << std::endl;
    cout << "cout" << std::endl;
    cerr << "cerr" << std::endl;
}

/*
#include "timer.hpp"

int main() {

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
    mc.flush_geometry_queue(std::clog, mc.resultqueue_faces_start, mc.result_normals, verts3, faces3);

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

    std::clog << "main();" << std::endl;
    return 0;
}
*/
