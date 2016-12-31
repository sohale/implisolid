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
- Improve JS data exchange to Geometry (done)
- Adding walls (done)
- WebWorker wrapped in a Geometry
- Geometry in 3js r73  (DONE  r79)
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
#include "../js_iteration_2/worker_call_specs.hpp"

#include "timer.hpp"

// #include <math.h>

extern "C" {
    void produce_object_old2(float* verts, int *nv, int* faces, int *nf, float param);
    int main(int argc, char **argv);
}

// Standard Emscripten compile option
#ifdef BUILD_AS_WORKER
    #define AS_WORKER true
#else
    #define AS_WORKER false
#endif

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

#include "../js_iteration_2/marching_cubes.hpp"
#include "tests/marching_cubes_mock.hpp"



/*********************************************************
    public interface for JavaScript
 *********************************************************/


extern "C" {
    void build_geometry(const char* shape_parameters_json, const char* mc_parameters_json);
    void build_geometry_u(const char* shape_parameters_json, const char* mc_parameters_json, const char* call_specs);
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
    int set_object(const char* shape_parameters_json, bool ignore_root_matrix);  // call object_factory. Sets a global variable   last_object
    bool unset_object(int);

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
}

// typedef    std::tuple<std::vector<REAL>, std::vector<vectorized_vect::index>> verts_faces_t;
// typedef    std::tuple<vectorized_vect, std::vector<vectorized_vect::index>> verts_faces_t;

typedef   std::tuple< vectorized_vect, vectorized_faces > vertsfaces_type;

class polygoniser__old {
    mp5_implicit::implicit_function const *  object_;

    vertsfaces_type vertsfaces;  // not good



    /*
    std::vector<REAL> verts;
    std::vector<int> faces;
    const vectorized_vect & centroids,
    */

private:
    polygoniser__old() = delete;
public:
    polygoniser__old(mp5_implicit::implicit_function const * ifunc_object)
        :   object_(ifunc_object)  // copy constructor
    {
        ;
    }
};


class state_t {
public:
    bool active = 0;
    // MarchingCubes* mc = 0;
    polygoniser__old pgonizer();  // why this does not generate an error?
    //todo: move this into this-> pgonizer
    std::vector<REAL> mc_result_verts;
    std::vector<vertexindex_type> mc_result_faces;

public:
    bool check_state() {
        return true;  // fixme. NOT WHEN UPDATING
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




//polygoniser_settings

// moved up
// std::map< std::string, vectorized_vect > point_set_set;

void* get_pointset_ptr(char* id) {
    if (point_set_set.find(std::string(id)) == point_set_set.end()) {
        return static_cast<void*>(0);  // Javascript interface mess
    }
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

std::pair< std::vector<REAL>, std::vector<vertexindex_type>>  mc_start (const mp5_implicit::implicit_function* object, dim_t resolution_, const mp5_implicit::bounding_box & box_/*, const bool use_metaball*/) {

    bool enableUvs = true;
    bool enableColors = true;

    MarchingCubes mc {resolution_, box_, enableUvs, enableColors};
    MarchingCubes * _state_mc = &mc;

    _state_mc -> isolation = 80.0/4*0;



    // ****************************
    // Does thismake things slower ?
    vectorized_vect   mcgrid_vectorized = _state_mc -> prepare_grid();  // 10.0
    _state_mc -> eval_shape(*object, mcgrid_vectorized);

    /*
     if (use_metaball) {
         REAL metaball_time = 0;
         meta_balls(*_state_mc, 4, metaball_time, 1.0);
     }
     */

    _state_mc->seal_exterior(-10000000.0);

    _state_mc->render_geometry();
    // std::clog << "MC executed" << std::endl;

    // std::clog << "map4" << std::endl;

    if (REPORT_STATS) {
    int mapctr = 0;
    for (auto& kv_pair : _state_mc->result_e3map) {
        if (0)
            std::clog << " [" << kv_pair.first << ':' << kv_pair.second << ']';
        mapctr++;
    }
    /*
    std::clog << "build_geometry(): ";
    std::clog << " e3Map: " << mapctr;
    std::clog << " Faces: " << _state_mc->result_faces.size()/3;
    std::clog << " Verts: " << _state_mc->result_verts.size()/3;
    std::clog << std::endl;
    */
    }

    return std::make_pair(_state_mc -> result_verts, _state_mc->result_faces);
}




std::pair< std::vector<REAL>, std::vector<vertexindex_type>> make_a_cube(REAL size) {
    std::vector<REAL> verts = {
        0,0,0, 1,0,0, 1,1,0, 0,1,0,
        0,0,1, 1,0,1, 1,1,1, 0,1,1};

    std::vector<vertexindex_type> faces = {
        0,2,1, 0,3,2,
        0,5,4, 0,1,5,
        1,6,5, 1,2,6,
        2,3,6, 6,3,7,
        4,3,0, 4,7,3,
        4,5,6, 4,6,7};

    for (auto& v : verts) {
        v -= 0.5;
        v *= size;
    }

    return std::make_pair(verts, faces);
}

std::pair< std::vector<REAL>, std::vector<vertexindex_type>> make_a_square(REAL size) {
    std::vector<REAL> verts = {
        0,0,0, 1,0,0, 1,1,0, 0,1,0
    };

    std::vector<vertexindex_type> faces = {
        0,2,1, 0,3,2,
    };

    for (auto& v : verts) {
        v -= 0.5;
        v *= size;
    }

    return std::make_pair(verts, faces);
}

#include "../js_iteration_2/polygonizer_algorithm_ob02.hpp"

inline void grand_algorithm(const std::string shape_parameters_json_str, const mp5_implicit::mc_settings  mc_settings_from_json, const worker_call_sepcs_t  worker_call_sepcs);

void build_geometry_u(const char* shape_parameters_json, const char* mc_parameters_json, const char* call_specs) {
    std::clog << "PROGRESSIVE" << std::endl;
    const mp5_implicit::mc_settings  mc_settings_from_json = parse_mc_properties_json(mc_parameters_json);
    std::string shape_parameters_json_str = std::string(shape_parameters_json);
    worker_call_sepcs_t  parse_worker_call_sepcs {std::string(call_specs)};
    grand_algorithm(shape_parameters_json_str, mc_settings_from_json, parse_worker_call_sepcs);
}

// void build_geometry(int resolution, char* mc_parameters_json, char* obj_name, REAL time){
void build_geometry(const char* shape_parameters_json, const char* mc_parameters_json)
{
    std::clog << "NON-Progressive" << std::endl;
    const mp5_implicit::mc_settings  mc_settings_from_json = parse_mc_properties_json(mc_parameters_json);
    std::string shape_parameters_json_str = std::string(shape_parameters_json);
    worker_call_sepcs_t  worker_call_sepcs {};
    grand_algorithm(shape_parameters_json_str, mc_settings_from_json, worker_call_sepcs);
}

/**
*************************************************
*        The grand algorithm                    *
*************************************************
*/
inline void grand_algorithm(
    const std::string  shape_parameters_json_str,
    const mp5_implicit::mc_settings  mc_settings_from_json,
    //const std::string  mc_settings_from_json,
    const worker_call_sepcs_t  worker_call_sepcs
) {

    if (!_state.check_state_null()) {
        clog << "build_geometry() called in a bad state.";
        return;
    }


    //const mp5_implicit::mc_settings  mc_settings_from_json = parse_mc_properties_json(mc_parameters_json);
    //const std::string shape_parameters_json_str = std::string(shape_parameters_json);
    bool ignore_root_matrix = mc_settings_from_json.ignore_root_matrix;

    //unique_pointer<mp5_implicit::implicit_function> object = ...;

    // bool use_metaball;  // output
    mp5_implicit::implicit_function* object = object_factory(shape_parameters_json_str , /*use_metaball,*/ ignore_root_matrix);

    // timer timr;
    //timr.report_and_continue("timer started.");

    //polygonizer  algorithm(_state, *object, std::string(mc_parameters_json)/*, timr*/);
    polygonizer algorithm(_state, *object, mc_settings_from_json /*, timr*/, worker_call_sepcs);


    // polygonizer::polygonize_init(_state, *object, mc_settings_from_json, use_metaball, algorithm.steps_report, timr);

    // polygonizer::polygonize_step_0(_state, *object, mc_settings_from_json, use_metaball, algorithm.steps_report, timr);
    algorithm.polygonize_step_0(/*use_metaball*/);

    algorithm.send_mesh_back_to_client();

    const bool DISABLE_POSTPROCESSING = false;    // DISABLE ALL MESH POST-PROCESSING (mesh optimisation)
    if (!DISABLE_POSTPROCESSING) {


        const int vresamp_iters  =  mc_settings_from_json.vresampl.iters; //10; //3;
        // disable hard-coded
        const bool apply_projection = true;

        #if NOTQUIET
        if (VERBOSE) mc_settings_from_json.report(clog);
        #endif

        const int overall_repeats = mc_settings_from_json.overall_repeats;
        for (int overall_iter = 0; overall_iter < overall_repeats; ++overall_iter) {
            for (int i=0; i < vresamp_iters; i++) {

                algorithm.polygonize_step_1();

            }
            algorithm.send_mesh_back_to_client();

            if (apply_projection) {
                /*
                if (STORE_POINTSETS)
                {
                    vectorized_vect v = convert_vectorverts_to_vectorized_vect( _state.mc_result_verts);
                    STORE_POINTSET("pre_p_centroid", v);
                }
                if (STORE_POINTSETS) {
                    STORE_POINTSET("pre_p_verts", v);
                }
                */

                if (mc_settings_from_json.projection.enabled) {

                    algorithm.polygonize_step_2();

                    algorithm.send_mesh_back_to_client();

                } else {
                    // std::clog << "centroids_projection (& qem) skipped because you asked for it." << std::endl;
                }

                const bool is_last = overall_iter == overall_repeats-1;
                if (mc_settings_from_json.subdiv.enabled
                    &&
                    (overall_repeats <= 1 || is_last )  // Skip the last one if overall_repeats > 1
                    )
                {
                    std::clog << "Subdivision:" << std::endl;

                    algorithm.polygonizer::polygonize_step_3(is_last);

                    //algorithm.send_mesh_back_to_client();

                } else {
                    std::clog << "subdivisions skipped because you didn't asked for it." << std::endl;
                }


                /*
                if (STORE_POINTSETS) {
                    vectorized_vect v = convert_vectorverts_to_vectorized_vect( _state.mc_result_verts);
                    STORE_POINTSET("post_p_centroids", v);
                }
                */
            }  // if (apply_projection)
        }  // for (overall_iter)
    }  // if (!DISABLE_POSTPROCESSING)

    /*
    // int REPEATS = 2;
    //int REPEATS_VR = 3;
    int REPEATS = 1;
    int REPEATS_VR = 1;
    for (int repeats = 0; repeats < REPEATS ; ++repeats) {
        float c = 2000.;

        for (int i=0; i < REPEATS_VR; i++){
            apply_vertex_resampling_to_MC_buffers__VMS(*object, c,  //    *(_state.mc));
                _state.mc_result_verts, _state.mc_result_faces, ??);
        }

        centroids_projection(object, _state.mc_result_verts, _state.mc_result_faces);
    }
    */

    // polygonizer::polygonize_terminate(_state, *object, mc_settings_from_json, algorithm.steps_report, timr);
    algorithm.polygonize_terminate();
    //delete object;
    object = NULL;
}

int get_f_size() {
    if (!_state.check_state()) return -1;
    return _state.mc_result_faces.size()/3;
}
int get_v_size() {
    if (!_state.check_state()) return -1;
    return _state.mc_result_verts.size()/3;
}
void get_v(REAL* v_out, int vcount) {
    if (!_state.check_state())
        return;

    // int nf = get_f_size();
    // Vertices
    int ctr = 0;
    for ( std::vector<REAL>::iterator it=_state.mc_result_verts.begin(); it < _state.mc_result_verts.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            v_out[ctr] = *(it + di);  // static_cast<REAL>()
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
    for ( std::vector<vertexindex_type>::iterator it=_state.mc_result_faces.begin(); it < _state.mc_result_faces.end(); it+=3 ) {
        for (int di=0; di < 3; di++) {
            f_out[ctr] = static_cast<int>(*(it + di));
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
    return reinterpret_cast<void*>(_state.mc_result_verts.data());
}

void* get_f_ptr() {
    if (!_state.check_state())
        return NULL;
    return reinterpret_cast<void*>(_state.mc_result_faces.data());
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
    if (!_state.check_state()) {
        std::cerr << "Error: finish_geometry() before producing the shape()" << std::endl;
        return;
    }
    assert (_state.active);

    // Dos not cause an exception if null. But it causes exception.
    // delete _state.mc;
    //_state.mc = 0;

    _state.active = false;
}

/* Set the default about string in case it is not provided by the compile script.
 */
#ifndef MORE_ABOUT_INFO

#define MORE_ABOUT_INFO  "-"

#endif


/** Show build information, e.g. when it was compiled.
 * Usage: Type in browser console:
 *      IMPLICIT.about()
 * Shows build time and date.
 * More information regarding the build script and build environment can be passed.
 * In compiler commandline (em++, g++, etc) use:
 *   -DMORE_ABOUT_INFO='"mytext"'
 * Useful for sending the hostname, IP, etc (if compiled within a docker container).
 * Note the extra '' necessary. Inside the '', use another ""
 * is used becasue it has to be a valid C++ expression (a C++ string).
 */
void about() {
    std::clog << "Build Info: " << std::endl;
    std::clog << __DATE__ << " " << __TIME__ << std::endl;

    std::clog << "Macros: "
        #if ASSERT_USED
            << "ASSERT_USED:yes. "
        #else
            << "ASSERT_USED:disabled. "
        #endif
        #if USE_ASSERT
            << "USE_ASSERT:yes. "
        #else
            << "USE_ASSERT:no. "
        #endif
        #if VERBOSE_SUBDIV
            << "VERBOSE_SUBDIV "
        #endif
        #if USE_PSDE
            << "USE_PSDE "
        #endif
        #if DEBUG_VERBOSE
            << "DEBUG_VERBOSE "
        #endif
        #if RANK_MANUALLY
            << "RANK_MANUALLY "
        #endif
        #if NDEBUG
            << "NDEBUG "
        #endif
        #if AS_WORKER
            << "AS_WORKER:yes. "
        #else
            << "AS_WORKER:no. "
        #endif
        << std::endl;

    std::clog << "CONFIG: "
        << "ROOT_TOLERANCE=" << ROOT_TOLERANCE << " "
        /*<< CONFIG*/
        << std::endl;

    std::clog << MORE_ABOUT_INFO << std::endl;
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

// Low-level API to directly query the implicit functions.
// Implicit_function, implicit solid (CSG), implicit_surfaces (BREP, etc)

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



int set_object(const char* shape_parameters_json, bool ignore_root_matrix) {
    if(ifunction_service.current_object != NULL){
        std::clog << "Error: You cannot unset() the object before a set_object(json)." << std::endl;
        return 0;
    }

    //std::clog << "before: current_object " << ifunction_service.current_object << std::endl;

    std::string str = std::string(shape_parameters_json);
    bool dummy;
    ifunction_service.current_object = object_factory(str, /*dummy,*/ ignore_root_matrix);

    //std::clog << "after: current_object " << ifunction_service.current_object << std::endl;
    int new_object_id = 1;
    return new_object_id;
}

bool unset_object(int id) {
    // id is ignored for now
    if(ifunction_service.current_object == NULL){
        std::clog << "Error: You cannot unset() the object before a set_object(json)." << std::endl;
        return false;
    }
    if (id <= 0){
        std::clog << "Incorrect ID: use the same id returned by set_object(json)." << std::endl;
        return false;
    }
    if (id != 1) {
        std::clog << "Incorrect ID. For now, The only id is 1" << std::endl;
        return false;
    }

    //delete ifunction_service.current_object;
    //ifunction_service.current_object = NULL;
    gc_objects();
    ifunction_service.current_object = NULL;
    return true;
}

bool set_x(void* verts, int n) {
    if(ifunction_service.current_x != NULL){
        std::clog << "Error: You set() before unset()ing the previous set()." << std::endl;
        return false;
    }
    if( n < 0 || n >= 10000 * 5) {
        std::clog << "Error: n is outside [0, 50000]." << std::endl;
        return false;
    }

    ifunction_service.current_x = new vectorized_vect( boost::extents[n][3] );
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
    ifunction_service.current_f = new vectorized_scalar( boost::extents[n]  );  // n x 0 (?)
    //std::clog << "warning: size may be n x 0:  " << ifunction_service.current_f->shape()[0] << "x" << ifunction_service.current_f->shape()[1] << std::endl;
    ifunction_service.current_grad = new vectorized_vect( boost::extents[n][3]  );
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

// void calculate_implicit_values_direct(void* start, int size) {}  // size = pointcount * 3
// void calculate_implicit_gradients_direct(void* start, int points_count) {} // size = pointcount * 3

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
    std::clog << "clog" << std::endl;
    std::cout << "cout" << std::endl;
    std::cerr << "cerr" << std::endl;
    #if AS_WORKER
    std::clog << "Compiled as a Web Worker" << std::endl;
    #else
    std::clog << "Not a a Web Worker" << std::endl;
    #endif
}


#if AS_WORKER
#include "mcc2_worker_api.hpp"
#endif
