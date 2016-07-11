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

#include "boost/multi_array.hpp"
#include "boost/array.hpp"

//#include <math.h>

extern "C" {
    void produce_object_old2(float* verts, int *nv, int* faces, int *nf, float param);
    int main();
}

const bool VERBOSE = false;
const bool REPORT_STATS = false;

//typedef unsigned short int size_t;
typedef unsigned short int dim_t; //small integers for example the size of one side of the grid
typedef float REAL;
//typedef unsigned long int index_t;

//boost::array will not work becasue the size of a boost:array has to be known in compile-time (static).
typedef  boost::multi_array<REAL, 1>  array1d;
//typedef array1d::index  array_index_t;
typedef boost::array<array1d::index, 1>  array_shape_t;
//#define array1d  boost::multi_array<REAL, 1>
typedef array1d::index  index_t;

typedef index_t index3_t; //Range of the element type has to be large enough, larger than (size^3)*3.
typedef boost::multi_array<index3_t, 1>   array1d_e3;
typedef std::map<index3_t,int>  e3map_t;


struct callback_t { void call (void*) const { } callback_t(){} };

/*template<typename Index_Type=int>
boost::array<Index_Type, 1> make_shape_1d(Index_Type size)
{
    //Make a shape to be used in array initialisation
    //fixme: the type of size

    //ASSERT(size>=0);
    boost::array<Index_Type, 1> shape = {{ size, }};
    return shape;
}
*/

#include "mcc2_marching_cubes.hpp"



void build_vf(
    //std::vector<REAL>& verts3,
    //std::vector<int>& faces3
    REAL scale
    ){
    // Includes allocations.

    dim_t resolution = 28;  // 28;
    bool enableUvs = true;
    bool enableColors = true;
    mp5_implicit::bounding_box box = {0,3,0,3,0,3};

    MarchingCubes mc(resolution, box, enableUvs, enableColors);


    int numblobs = 4;
    REAL time = 0.1 ;
    for (int ball_i = 0; ball_i < numblobs; ball_i++) {
        REAL ballx = sin(ball_i + 1.26 * time * (1.03 + 0.5*cos(0.21 * ball_i))) * 0.27 + 0.5;
        REAL bally = std::abs(cos(ball_i + 1.12 * time * cos(1.22 + 0.1424 * ball_i))) * 0.77; // dip into the floor
        REAL ballz = cos(ball_i + 1.32 * time * 0.1*sin((0.92 + 0.53 * ball_i))) * 0.27 + 0.5;
        REAL subtract = 12;
        REAL strength = 1.2 / ((sqrt(numblobs)- 1) / 4 + 1);
        mc.addBall(ballx, bally, ballz, strength, subtract, scale);
    }
    mc.seal_exterior();

    /*
    int numblobs = 4;
    REAL subtract = (REAL)12.;
    REAL strength = (REAL)(1.2 / ( ( sqrt( numblobs ) - 1. ) / 4. + 1. ));
    mc.addBall(0.5, 0.5, 0.5, strength, subtract, 1);
    */
    //MarchingCubes& object = mc;
    //mc.addBall(0.5, 0.5, 0.5, strength, subtract, 1);

    //mc.flush_geometry_queue(std::cout, mc.resultqueue_faces_start, mc.result_normals, verts3, faces3);

    const callback_t renderCallback;
    mc.render_geometry(renderCallback);

    if(VERBOSE)
        std::cout << "MC:: v,f: " << mc.result_verts.size() << " " << mc.result_faces.size() << std::endl;

    //verts3.resize(0);
    //faces3.resize(0);
}


class MarchingCubesMock {

public:
    MarchingCubesMock( dim_t resolution, bool enableUvs, bool enableColors ) {};
    ~MarchingCubesMock() {}; //why does this have to be public: ?

    //void flush_geometry_queue(std::ostream&);
    void flush_geometry_queue(std::ostream& cout, int& normals_start, std::vector<REAL> &normals,  std::vector<REAL> &verts3, std::vector<int> &faces3, e3map_t &e3map, int& next_unique_vect_counter)
        {};

    inline int polygonize_cube( REAL fx, REAL fy, REAL fz, index_t q, REAL isol, const callback_t& callback ) {return 0;};

    REAL isolation;

//shape:
    void addBall( REAL ballx, REAL bally, REAL ballz, REAL strength, REAL subtract, REAL scale ) {};
    void addPlaneX( REAL strength, REAL subtract ) {};
    void addPlaneZ( REAL strength, REAL subtract ) {};
    void addPlaneY( REAL strength, REAL subtract ) {};
    void seal_exterior(const REAL exterior_value) {};
//field
    void reset() {};

//geometry/threejs interface side.
    void render_geometry(const callback_t& renderCallback ) {};
    void sow() {};

// output. filled using sow()
    int resultqueue_faces_start = 0;
    std::vector<REAL> result_normals;
    std::vector<REAL> result_verts;
    std::vector<int> result_faces;
};


//Not used
void produce_object_old2(REAL* verts, int *nv, int* faces, int *nf, REAL time, REAL scale){
    //not used

    dim_t resolution = 28;  // 28;
    bool enableUvs = true;
    bool enableColors = true;

    mp5_implicit::bounding_box box = {0,3,0,3,0,3};

    if(VERBOSE)
        std::cout << "Leak-free (old version)" << std::endl;


    MarchingCubes mc(resolution, box, enableUvs, enableColors);
    //MarchingCubes* mc0 = new MarchingCubes(resolution, enableUvs, enableColors);
    //MarchingCubes &mc = *mc0;
    //MarchingCubesMock mc(resolution, enableUvs, enableColors);

    int numblobs = 4;
    //REAL time = 0.1 ;
    for (int ball_i = 0; ball_i < numblobs; ball_i++) {
        REAL D = 1;
        REAL ballx = sin(ball_i + 1.26 * time * (1.03 + 0.5*cos(0.21 * ball_i))) * 0.27 * D + 0.5   ;
        REAL bally = std::abs(cos(ball_i + 1.12 * time * cos(1.22 + 0.1424 * ball_i))) * 0.77 * D; // dip into the floor
        REAL ballz = cos(ball_i + 1.32 * time * 0.1*sin((0.92 + 0.53 * ball_i))) * 0.27 * D+ 0.5;
        REAL subtract = 12;
        REAL strength = 1.2 / ((sqrt(numblobs)- 1) / 4 + 1);
        mc.addBall(ballx, bally, ballz, strength, subtract, scale);
    }
    mc.seal_exterior();

    /*
    int numblobs = 4;
    REAL subtract = (REAL)12.;
    REAL strength = (REAL)(1.2 / ( ( sqrt( numblobs ) - 1. ) / 4. + 1. ));
    mc.addBall(0.5, 0.5, 0.5, strength, subtract, scale);

    mc.addBall(0, 0, 0.5, strength, subtract, scale);
    */

    //todo: init-receiver side

    const callback_t renderCallback;
    mc.render_geometry(renderCallback);

    std::cout << "map2" << std::endl;

    //mc.result_faces.resize(100);

    if(VERBOSE)
        std::cout << "MC:: v,f: " << mc.result_verts.size() << " " << mc.result_faces.size() << std::endl;

    *nv = mc.result_verts.size()/3;
    *nf = mc.result_faces.size()/3;

    // Vertices
    int ctr = 0;
    for(std::vector<REAL>::iterator it=mc.result_verts.begin(); it < mc.result_verts.end(); it+=3 ){
        for(int di=0; di<3; di++){
            verts[ctr] = *( it + di );
            ctr++;
        }
      }

    // Faces
    ctr = 0;
    for(std::vector<int>::iterator it=mc.result_faces.begin(); it < mc.result_faces.end(); it+=3 ){
        for(int di=0; di<3; di++){
            faces[ctr] = *( it + di );
            ctr++;
        }
      }


}


/*
#include <emscripten/bind.h>
using namespace emscripten;
EMSCRIPTEN_BINDINGS(my_module) {
    //produce_object(REAL* verts, int *nv, int* faces, int *nf, REAL time);
    function("produce_object_bind", &produce_object);
}
*/
//#include "mcc1-glue.cpp"


/*void produce_v(){
}*/

extern "C" {
    void build_geometry(int resolution, REAL mc_size, char* obj_name, REAL time);
    int get_v_size();
    int get_f_size();
    void get_f(int*, int);
    void get_v(REAL*, int);
    void finish_geometry();
    void* get_f_ptr();
    void* get_v_ptr();
    void implicit_value(REAL x, REAL y, REAL z);
    void implicit_grad_value(REAL x, REAL y, REAL z, REAL * x_out, REAL * y_out, REAL * z_out);
    //also: queue, etc.
    //bad: one instance only.
    //    Solution 1:  MarchingCubes* build_geometry();
    //    Solution 2: ids (for workers! ; a statically determined number of them (slots/workers/buckets).).
};


typedef struct {
    bool active = 0;
    MarchingCubes* mc = 0;
} state_t;

state_t _state;

//_state.active = false;
//_state.mc = 0;

bool check_state() {
    if(!_state.active){
        std::cout << "Error: There are no allocated geometry resources to deallocate.";
        return false;
    }
    return true;
}
bool check_state_null() {
    if(_state.active){
        std::cout << "Error: There are non-deallocated geometry resources. Call finit_geometry() first. Not doing anything.";
        return false;
    }
    return true;
}


void meta_balls(MarchingCubes& mc, REAL time, REAL scale){

    int numblobs = 4;
    for (int ball_i = 0; ball_i < numblobs; ball_i++) {
        REAL ballx = sin(ball_i + 1.26 * time * (1.03 + 0.5*cos(0.21 * ball_i))) * 0.27 + 0.5;
        REAL bally = std::abs(cos(ball_i + 1.12 * time * cos(1.22 + 0.1424 * ball_i))) * 0.77; // dip into the floor
        REAL ballz = cos(ball_i + 1.32 * time * 0.1*sin((0.92 + 0.53 * ball_i))) * 0.27 + 0.5;
        REAL subtract = 12;
        REAL strength = 1.2 / ((sqrt(numblobs)- 1) / 4 + 1) ;
        mc.addBall(ballx, bally, ballz, strength, subtract, scale);
    }

}
//#include "../js_iteration_2/unit_sphere.hpp"
#include "../js_iteration_2/primitives.cpp"
#include "../js_iteration_2/crisp_subtract.hpp"
//using namespace mp5_implicit;


implicit_function*  object_factory(REAL f_argument, std::string name){
    implicit_function* object;
    if (name == "double_mushroom"){
        object = new mp5_implicit::double_mushroom(0.8, 1/(f_argument+3), 1/(f_argument+3), f_argument+3);
    }
    else if (name == "egg"){
        object = new mp5_implicit::egg(f_argument,f_argument,f_argument);
    }
    else if (name == "sphere"){
        object = new mp5_implicit::unit_sphere((sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3)*10);
    }
    else if (name == "cube"){
        object = new mp5_implicit::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
    }
    else if (name == "super_bowl"){// not working
        object = new mp5_implicit::super_bowl(1.5/(f_argument+3.0));
    }
    else if (name == "scone"){
        object = new mp5_implicit::scone(f_argument +2.5,f_argument +2.5,f_argument +2.5,-0.1);
    }
    else if (name == "scylinder"){
        object = new mp5_implicit::scylinder(f_argument, 1.6); //0.7
    }else if(name == "meta_balls"){
        REAL r = (sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3)*10;
        std::cout << " META BALLS r : " << r << std::endl;
        object = new mp5_implicit::unit_sphere(r);
    }
    else if(name == "sub_spheres"){
        mp5_implicit::unit_sphere * s1 = new mp5_implicit::unit_sphere(2, 1, 1, 1);
        mp5_implicit::unit_sphere * s2 = new mp5_implicit::unit_sphere(1.3);
        object = new mp5_implicit::CrispSubtract(*s1, *s2);
    }
    else {
        std::cout << "Error! You must enter a valid name! So I made a sphere!" << std::endl;
        object = new mp5_implicit::unit_sphere(sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3);
    }
    return object;
}
void build_geometry(int resolution, REAL mc_size, char* obj_name, REAL time){

    if(!check_state_null())
        return;
    std::cout << "In build_geometry obj_name : " << obj_name << std::endl;
    //dim_t resolution = 28;
    bool enableUvs = true;
    bool enableColors = true;
    mp5_implicit::bounding_box box = {0,3,0,3,0,3};//{15,20,15,20,15,20};
    //std::cout << "Leak-free : new" << std::endl;

    //MarchingCubes mc(resolution, enableUvs, enableColors);
    _state.mc = new MarchingCubes(resolution, box, enableUvs, enableColors);
    //std::cout << "constructor called. " << _state.mc << std::endl;

    _state.mc -> isolation = 80.0/4*0;

    /*
    mp5_implicit :: unit_sphere   object(sin(0.033*10 * time * 3.1415*2.)*0.33+0.3);
    // //_state.mc -> prepare_grid(1.0);
    // //object.eval_implicit(grid, implicit_values);
    _state.mc -> eval_shape(object, 1.0);
    */



    std::string name = std::string(obj_name);
    std::cout << "Name : " << name << std::endl;

    REAL f_argument = time;

    implicit_function* object = object_factory(f_argument, name);


    // ****************************
    // Does thismake things slower ?
    boost::multi_array<REAL, 2>  mcgrid_vectorized = _state.mc -> prepare_grid();  // 10.0
    _state.mc -> eval_shape(*object, mcgrid_vectorized);

     delete object;
     object = NULL;
     if(name == "meta_balls"){
         meta_balls(*_state.mc, time, 1.0);
     }

    _state.mc->seal_exterior();

    //std::cout << "balls added." << std::endl;

    const callback_t renderCallback;
    _state.mc->render_geometry(renderCallback);
    //std::cout << "MC executed" << std::endl;

    //std::cout << "map4" << std::endl;

    if(REPORT_STATS){
    int mapctr = 0;
    for (auto& kv_pair: _state.mc->result_e3map){
        if(0)
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


    if(VERBOSE){
        std::cout << resolution << " " << time << std::endl;
        std::cout << _state.mc << std::endl;
    }
    _state.active = true;

    check_state();
    //std::cout << "MC:: v,f: " << _state.mc->result_verts.size() << " " << _state.mc->result_faces.size() << std::endl;
}
int get_f_size() {
    if(!check_state()) return -1;
    return _state.mc->result_faces.size()/3;
}
int get_v_size(){
    if(!check_state()) return -1;
    return _state.mc->result_verts.size()/3;
}
void get_v(REAL* v_out, int vcount){
    if(!check_state())
        return;

    //int nf = get_f_size();
    // Vertices
    int ctr = 0;
    for(std::vector<REAL>::iterator it=_state.mc->result_verts.begin(); it < _state.mc->result_verts.end(); it+=3 ){
        for(int di=0; di<3; di++){
            v_out[ctr] = *( it + di );
            //if(ctr<3*3*3)
            //    std::cout << v_out[ctr] << " ";
            ctr++;
        }
    }
    //std::cout << std::endl;
    //assert nf*3 == ctr;
    if(vcount*3 != ctr)  std::cout << "sizes dont match: " << (float)ctr/3. << " " << vcount << std::endl;
}

void get_f(int* f_out, int fcount){
    if(!check_state())
        return;
    //int nf = get_f_size();
    int ctr = 0;
    for(std::vector<int>::iterator it=_state.mc->result_faces.begin(); it < _state.mc->result_faces.end(); it+=3 ){
        for(int di=0; di<3; di++){
            f_out[ctr] = *( it + di );
            //if(ctr<3*3*3)
            //    std::cout << f_out[ctr] << " ";
            ctr++;
        }
    }
    if(fcount*3 != ctr)  std::cout << "sizes dont match: " << (float)ctr/3. << " " << fcount << std::endl;
    //std::cout << std::endl;
};

/* Data is already there, so why copy it? Also, keep it until next round. */
void* get_v_ptr(){
    if(!check_state())
        return NULL;
    return (void*)(_state.mc->result_verts.data());
}

void* get_f_ptr(){
    if(!check_state())
        return NULL;
    return (void*)(_state.mc->result_faces.data());
}


//int get_v_size(){};
//int get_f_size(){};
//void get_f(int*){};
//void get_v(REAL*){};
//void* get_f_ptr();
//void* get_v_ptr();
//void finish_geometry();

// Can cause an exception (but not when null).
void finish_geometry() {
    if(!check_state())
        return ;
    if(_state.mc == 0){
        std::cerr << "Error: finish_geometry() before producing the shape()" << std::endl;
    }
    if(!_state.active){
        //std::cout << "Cannot finish_geometry() while still active." << std::endl;
    }
    else{
        //std::cout << "_state.active " << _state.active << "  _state.mc " << _state.mc << std::endl;
    }
    //Dos not cause an exception if null. But it causes exception.
    delete _state.mc;
    _state.active = false;
    _state.mc = 0;
};


#include "timer.hpp"

int main() {
    /*
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

    //int normals_start = 0;
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

    //build_vf( verts3, faces3 );  // 21.3 msec using O3
    build_vf(  );  // 26 msec.


    t.stop();
*/
    std::cout << "main();" << std::endl;
    return 0;
}
