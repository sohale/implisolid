#pragma once

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include "implicit_function/primitives.hpp"

//#include""  will be moved into sweep.hpp, extrude.hpp, screw.hpp, etc
#include "implicit_function/2d/primitives_2d.hpp"

#include "object_collector.hpp"
//using namespace mp5_implicit;
using mp5_implicit::implicit_function;


namespace pt = boost::property_tree ;

void getMatrix12(REAL * matrix12, const pt::ptree& shapeparams_dict){
        int i = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("matrix")){

            REAL x = element.second.get_value<REAL>();
            //std::clog << "matrix value : " << x << std::endl;
            matrix12[i] = x;
            i++;
        }
}

void getCorners(std::vector<boost::array<REAL,3>>& corners, const pt::ptree& shapeparams_dict){

    int i = 0;
    for (const pt::ptree::value_type &element : shapeparams_dict.get_child("corners")) {
        int j = 0;
        for (const pt::ptree::value_type &cell : element.second)
        {
            corners[i][j] = cell.second.get_value<REAL>();
            j++;
        }
        i++;

    }
}

void copy_eye(REAL matrix12[12]){
    REAL eye[12] = {1,0,0,0,  0,1,0,0,  0,0,1,0 };
    for(int j=0;j<12;j++)
        matrix12[j] = eye[j];
}

// Curated implicit_functions only
implicit_function*  object_factory(pt::ptree shapeparams_dict, bool ignore_root_matrix) {
    // std::clog << "ignore_root_matrix: " << ignore_root_matrix << std::endl;
    std::string name = shapeparams_dict.get<std::string>("type");
    //std::cout << "-------------------------------------" <<std::endl;
    //std::cout << name <<std::endl;


    implicit_function* object;

    // name = "meta_balls";  // meta_balls demo

    if (name == "implicit_double_mushroom"){
        // std::clog << "implicit_double_mushroom case " << std::endl;
        // object = new mp5_implicit::double_mushroom(0.9, 0.4 ,0.4 , 1/0.2 );
        implicit_function* dm = new mp5_implicit::double_mushroom(0.9, 0.4/2, 0.4/2, 1/0.2 );

        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new mp5_implicit::linearly_transformed(dm, matrix12);
        register_new_object(object);
        //object = dm;
    }
    else
    /*if (name == "simple_sphere"){
        //std::clog << "******************* simple_sphere case " << std::endl;
        REAL radius = shapeparams_dict.get<REAL>("radius");
        object = new mp5_implicit::unit_sphere(radius);
        std::clog << "radius " << radius << std::endl;
        register_new_object(object);
    }
    else*/
    if (name == "icube" || name == "cube" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }
        object = new mp5_implicit::cube(matrix12);
       // object = new mp5_implicit::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
        register_new_object(object);
    }
    else
    if (name == "icylinder" || name == "cylinder" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new mp5_implicit::scylinder(matrix12);
        register_new_object(object);
       // object = new mp5_implicit::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
        //register_new_object(object);
    }
    else
    if (name == "iellipsoid" || name == "ellipsoid" ){
        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new mp5_implicit::egg(matrix12);
        register_new_object(object);
    }else
    if (name == "icone" || name == "cone" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new mp5_implicit::scone(matrix12);
        register_new_object(object);
    }else if(name == "iheart" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new mp5_implicit::heart(matrix12);
        register_new_object(object);

    }else if(name == "itorus" ){
        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new mp5_implicit::torus(matrix12);
        register_new_object(object);

    }
    else
    if (name == "tetrahedron"){

        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);

        std::vector<boost::array<REAL,3>> corners(4);
        getCorners(corners, shapeparams_dict);

        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }
        object = new mp5_implicit::tetrahedron(corners, matrix12);
        register_new_object(object);
    }
    else if (name == "screw") {

        Matrix<REAL, 4, 4> transformation_matrix;
        REAL pitch = 0.0;
        std::string profile;
        std::string end_type;
        REAL delta_ratio = 0.0;
        Matrix<REAL, 3, 1> v;

        mp5_implicit::screw::getScrewParameters(
                           transformation_matrix, pitch, profile, 
                           end_type, delta_ratio, v, 
                           shapeparams_dict);

        transformation_matrix << 1, 0, 0, 0,
                                 0, 1, 0, 0,
                                 0, 0, 1, 0,
                                 0, 0, 0, 1;

        object = new mp5_implicit::screw(transformation_matrix,
                                         pitch,
                                         profile,
                                         end_type,
                                         delta_ratio,
                                         v);

        register_new_object(object);
    }
    else if (name == "Union") {
        //todo: Use SimpleUnion if (matrix12 == eye(4))
        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        implicit_function * a = NULL;

        implicit_function * o_matrix;
        implicit_function * o_plain;

        for (pt::ptree::value_type &element : shapeparams_dict.get_child("children")) {
            if (a == NULL){
                a = object_factory(element.second, false);
            } else {
                implicit_function * b = object_factory(element.second, false);

                //The following always prints an empty line:
                //std::clog << "element.second.get_value<string>(\"type\")" << element.second.get_value<string>("type") << std::endl ;

                //a = new mp5_implicit::CrispUnion(*a, *b);
                // register_new_object(a);
                //std::vector<const implicit_function*> versus std::vector<implicit_function*>
                std::vector<implicit_function*> ab = std::vector<implicit_function*>();
                ab.push_back(a);
                ab.push_back(b);
                o_matrix = new mp5_implicit::transformed_union(ab, matrix12);
                register_new_object(o_matrix);
                REAL eye_matrix12[12] = {1,0,0,0,  0,1,0,0,  0,0,1,0};
                o_plain =  new mp5_implicit::transformed_union(ab, eye_matrix12);
                a = o_plain;
                register_new_object(o_plain);
            }

            //std::clog << "##### " << element.second.get_child("type") << std::endl;
            //i++;
        }
        a = o_matrix;
        //std::clog  << "#####" << shapeparams_dict.get<string>("children") << std::endl;

        object = a;
        // register_new_object
    }else if (name == "Intersection") {
        //todo: Use SimpleUnion if (matrix12 == eye(4))
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        implicit_function * a = NULL;
        implicit_function * b = NULL;
        int count = 0;
        for (pt::ptree::value_type &element : shapeparams_dict.get_child("children")) {
            if (a==NULL){
                a = object_factory(element.second, false);
            }else{
                if(count > 1){
                    std::clog << "An CrispIntersection should have only 2 child" << std::endl;
                    break;
                }
                b = object_factory(element.second, false);

            }
            count++;
            if (count > 2) {
                std::clog << "Error: Intersection cannot be applied to more than two objects." << std::endl;
            }
        }

        // object = new mp5_implicit::CrispIntersection(*a, *b);
        //object = new mp5_implicit::transformed_intersection(*a, *b, matrix12);
        std::vector<implicit_function *> children;
        children.push_back(a);
        children.push_back(b);
        object = new mp5_implicit::transformed_intersection(children, matrix12);
        register_new_object(object);

    }else if (name == "Difference") {
        //todo: Use SimpleUnion if (matrix12 == eye(4))
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        implicit_function * a = NULL;
        implicit_function * b = NULL;
        int count = 0;
        for (pt::ptree::value_type &element : shapeparams_dict.get_child("children")) {
            if (a==NULL){
                a = object_factory(element.second, false);
            }else{
                if(count > 1){
                    std::clog << "An CrispSubstraction should have only 2 child" << std::endl;
                    break;
                }
                b = object_factory(element.second, false);

            }
            count++;
            if (count > 2) {
                std::clog << "Error: Intersection cannot be applied to more than two objects." << std::endl;
            }

        }

        // object = new mp5_implicit::CrispSubtract(*a, *b);
        // object = new mp5_implicit::transformed_subtract(*a, *b, matrix12);
        std::vector<implicit_function *> children;
        children.push_back(a);
        children.push_back(b);
        object = new mp5_implicit::transformed_subtract(children, matrix12);
        register_new_object(object);

    } else if(name == "meta_balls") {

        //REAL r = (sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3)*1;
        //std::clog << " META BALLS r : " << r << std::endl;
        // (REAL matrix12[12], int num_blobs, REAL time, REAL scale)

        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        int num_blobs = 4;
        //REAL time = 0.1; //0;  // 0.1
        //REAL time = shapeparams_dict.get_value<REAL>("time", 0.1);   // get_value<REAL>() or get<REAL>()
        REAL time = shapeparams_dict.get<REAL>("time", 0.1);   // get_value<REAL>() or get<REAL>()
        REAL scale = 1.0;
        object = new mp5_implicit::meta_ball_Rydgård(matrix12, num_blobs, time, scale);

        register_new_object(object);
    } else {
        std::cerr << "Invalid object " << "you asked for: \"" << name << "\"" << std::endl;
        abort();
    }

    return object;

}



implicit_function*  object_factory(string shape_parameters_json, bool ignore_root_matrix)
{
    std::stringstream shape_json_stream;
    shape_json_stream << shape_parameters_json ;
    pt::ptree shapeparams_dict;
    pt::read_json(shape_json_stream, shapeparams_dict);
    return object_factory(shapeparams_dict, ignore_root_matrix);

}
