// code moved from three files:
//    mcc2.cpp
// object_factory.hpp


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
void meta_balls_old(MarchingCubes& mc, int num_blobs, REAL time, REAL scale) {
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
    meta_balls_old(mc, 4, time, scale);
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
    meta_balls_old(mc, 4, time, scale);
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
    for ( std::vector<vertexindex_type>::iterator it=mc.result_faces.begin(); it < mc.result_faces.end(); it+=3 ) {
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



///////////////////////////////////////////////////////////////////////////////////////////////
// code moved from object_factory.hpp




/*

Non-Curated implicit_function s

implicit_function*  object_factory_simple(REAL f_argument, std::string name){
    std::cout << "This method is deprecated. Never call it" << std::endl;
    implicit_function* object;
    
    if (name == "double_mushroom"){
        object = new mp5_implicit::double_mushroom(0.8, 1/(f_argument+3), 1/(f_argument+3), f_argument+3);
        register_new_object(object);
    }
    else if (name == "egg"){
        object = new mp5_implicit::egg(f_argument,f_argument,f_argument);
        register_new_object(object);
    }
    else if (name == "sphere"){
        object = new mp5_implicit::unit_sphere((sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3)*10);
        register_new_object(object);
    }
    else if (name == "cube"){
        object = new mp5_implicit::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
        register_new_object(object);
    }
    else if (name == "super_bowl"){// not working
        object = new mp5_implicit::super_bowl(1.5/(f_argument+3.0));
        register_new_object(object);
    }
    else if (name == "scone"){
        object = new mp5_implicit::scone(f_argument +2.5,f_argument +2.5,f_argument +2.5,-0.1);
        register_new_object(object);
    }
    //else if (name == "scylinder"){
    //    object = new mp5_implicit::scylinder(f_argument, 1.6); //0.7
    //    register_new_object(object);
    //}
    //else

    if(name == "meta_ball (NOT!)"){
        REAL r = (sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3)*1;
        std::clog << " META BALLS r : " << r << std::endl;
        object = new mp5_implicit::unit_sphere(r);
        register_new_object(object);
    }


    else if(name == "sub_spheres"){
        mp5_implicit::unit_sphere * s1 = new mp5_implicit::unit_sphere(2, 1, 1, 1);
        mp5_implicit::unit_sphere * s2 = new mp5_implicit::unit_sphere(1.3);
        object = new mp5_implicit::CrispSubtract(*s1, *s2);
        register_new_object(s1);
        register_new_object(s2);
        register_new_object(object);
    }

    else {
        std::clog << "Error! You must enter a valid name " <<  name << "! So I made a sphere!" << std::endl;
        //object = new mp5_implicit::unit_sphere(sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3);
        object = new mp5_implicit::unit_sphere(1000.);
        register_new_object(object);
    }
    return object;
}
*/




//REAL xmax = shapeparams_dict.get<REAL>("matrix",NaN);
//std::clog << "############Name : " << name << std::endl;
//REAL zmax = shapeparams_dict.get<REAL>("box.zmax",NaN);
// int resolution = shapeparams_dict.get<int>("resolution",-1);

// if(isNaN(xmin) || isNaN(xmax) || isNaN(ymin) || isNaN(ymax) || isNaN(zmin) || isNaN(zmax) || resolution <= 2 ){
//     std::clog << "Error: missing or incorrect values in mc_parameters_json"<< std::endl;
//     xmin = -1;
//     xmax = 1;
//     ymin = -1;
//     ymax = 1;
//     zmin = -1;
//     zmax = 1;
//     resolution = 28;
// }



bool use_metaball = false;
/*
if(name=="meta_balls"){
    use_metaball = true;
}
else{
    use_metaball = false;
}
*/





//else
        /*
        // if(name=="meta_balls")
        REAL f_argument = shapeparams_dict.get<REAL>("time", NaN);

        std::clog << "otherwise " << "you asked for " << name << std::endl;
        object = object_factory_simple(f_argument, name);
        */




        /*
        mp5_implicit :: unit_sphere   object(sin(0.033*10 * time * 3.1415*2.)*0.33+0.3);
        // //_state.mc -> prepare_grid(1.0);
        // //object.eval_implicit(grid, implicit_values);
        _state.mc -> eval_shape(object, 1.0);
        */


    /*
        std::string name = std::string(obj_name);
        std::clog << "Name : " << name << std::endl;

        REAL f_argument = time;

        implicit_function* object = object_factory_simple(f_argument, name);
    */
        //std::clog << "############################" << shape_parameters_json << std::endl;



        //#include "implicit_function/crisp_subtract.hpp"
        //#include "implicit_function/linearly_transformed.hpp"



struct callback_t { void call(void*) const { } callback_t(){} };

