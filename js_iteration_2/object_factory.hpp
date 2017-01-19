#pragma once

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include "implicit_function/primitives.hpp"

//#include""  will be moved into sweep.hpp, extrude.hpp, screw.hpp, etc
#include "implicit_function/2d/primitives_2d.hpp"

#include "object_collector.hpp"
#include "implicit_function/2d/implicit_function_2d.hpp"
#include "implicit_function/2d/GDT/convex_polygon.hpp"
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

namespace mp5_implicit {

// Curated implicit_functions only
implicit_function*  object_factory(pt::ptree shapeparams_dict, bool ignore_root_matrix) {

    /*
    using mp5_implicit::implicit_functions::double_mushroom;
    using mp5_implicit::implicit_functions::linearly_transformed;
    using mp5_implicit::implicit_functions::cube;
    using mp5_implicit::implicit_functions::scylinder;
    using mp5_implicit::implicit_functions::egg;
    using mp5_implicit::implicit_functions::scone;
    using mp5_implicit::implicit_functions::heart;
    using mp5_implicit::implicit_functions::torus;
    using mp5_implicit::implicit_functions::tetrahedron;
    using mp5_implicit::implicit_functions::meta_ball_Rydgård;
    using mp5_implicit::implicit_functions::transformed_subtract;
    using mp5_implicit::implicit_functions::transformed_intersection;
    using mp5_implicit::implicit_functions::transformed_union;
    using mp5_implicit::implicit_functions::screw;
    */

    // std::clog << "ignore_root_matrix: " << ignore_root_matrix << std::endl;
    std::string name = shapeparams_dict.get<std::string>("type");
    //std::cout << "-------------------------------------" <<std::endl;
    //std::cout << name <<std::endl;


    implicit_function* object;

    // name = "meta_balls";  // meta_balls demo

    if (name == "implicit_double_mushroom"){
        // std::clog << "implicit_double_mushroom case " << std::endl;
        // object = new implicit_functions::double_mushroom(0.9, 0.4 ,0.4 , 1/0.2 );
        implicit_function* dm = new implicit_functions::double_mushroom(0.9, 0.4/2, 0.4/2, 1/0.2 );

        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new implicit_functions::linearly_transformed(dm, matrix12);
        register_new_object(object);
        //object = dm;
    }
    else
    /*if (name == "simple_sphere"){
        //std::clog << "******************* simple_sphere case " << std::endl;
        REAL radius = shapeparams_dict.get<REAL>("radius");
        object = new implicit_functions::unit_sphere(radius);
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
        object = new implicit_functions::cube(matrix12);
       // object = new implicit_functions::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
        register_new_object(object);
    }
    else
    if (name == "icylinder" || name == "cylinder" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new implicit_functions::scylinder(matrix12);
        register_new_object(object);
       // object = new implicit_functions::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
        //register_new_object(object);
    }
    else
    if (name == "iellipsoid" || name == "ellipsoid" ){
        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new implicit_functions::egg(matrix12);
        register_new_object(object);
    }else
    if (name == "icone" || name == "cone" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new implicit_functions::scone(matrix12);
        register_new_object(object);
    }else if(name == "iheart" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new implicit_functions::heart(matrix12);
        register_new_object(object);

    }else if(name == "itorus" ){
        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        object = new implicit_functions::torus(matrix12);
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
        object = new implicit_functions::tetrahedron(corners, matrix12);
        register_new_object(object);

    } else if (name == "inf_screw") {

        Matrix<REAL, 3, 4> transformation_matrix;
        REAL pitch = 0.0;
        std::string profile;
        std::string end_type;
        REAL delta_ratio = 0.0;
        Matrix<REAL, 3, 1> v;

        implicit_functions::screw::getScrewParameters(
                           transformation_matrix, pitch, profile,
                           end_type, delta_ratio, v,
                           shapeparams_dict);

        if(ignore_root_matrix) {
            transformation_matrix << 1, 0, 0, 0,
                                     0, 1, 0, 0,
                                     0, 0, 1, 0;
        }        
        object = new implicit_functions::screw(transformation_matrix,
                                         pitch,
                                         profile,
                                         end_type,
                                         delta_ratio,
                                         v);

        register_new_object(object);
    } else if (name == "screw_diff_two_plane") {

        // dangerous
        // ignore_root_matrix = !ignore_root_matrix;


        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }
        REAL eye[12] = {1,0,0,0,  0,1,0,0,  0,0,1,0 };

        Matrix<REAL, 3, 4> transformation_matrix;
        REAL pitch = 0.0;
        std::string profile;
        std::string end_type;
        REAL delta_ratio = 0.0;
        Matrix<REAL, 3, 1> v;

        implicit_functions::screw::getScrewParameters(
                           transformation_matrix, pitch, profile,
                           end_type, delta_ratio, v,
                           shapeparams_dict);

        if(true) {
            transformation_matrix << 1, 0, 0, 0,
                                     0, 1, 0, 0,
                                     0, 0, 1, 0;
        }        
        // object = new implicit_functions::screw(transformation_matrix,
        //                                  pitch,
        //                                  profile,
        //                                  end_type,
        //                                  delta_ratio,
        //                                  v);

        // register_new_object(object);

        // the following code is taken from difference

        implicit_function * a = NULL;
        implicit_function * b = NULL;
        a = new implicit_functions::screw(transformation_matrix,
                                                         pitch,
                                                         profile,
                                                         end_type,
                                                         delta_ratio,
                                                         v);

        Eigen::Matrix<REAL, 3, 1> plane_vector;
        Eigen::Matrix<REAL, 3, 1> plane_point;
        plane_vector << 0.0, 0.0, 1.0;
        plane_point << 0.0, 0.0, 0.25;

        implicit_function * hp = NULL;
        hp = new half_plane(transformation_matrix,
                             plane_vector,
                             plane_point);

        std::vector<implicit_function *> children_first;
        children_first.push_back(a);
        children_first.push_back(hp);

        implicit_function* first_object;

        first_object = new implicit_functions::transformed_subtract(children_first, eye);


        Eigen::Matrix<REAL, 3, 1> plane_vector_bottom;
        Eigen::Matrix<REAL, 3, 1> plane_point_bottom;
        plane_vector_bottom << 0.0, 0.0, -1.0;
        plane_point_bottom << 0.0, 0.0, -0.25;

        implicit_function * hp_bottom = NULL;
        hp_bottom = new half_plane(transformation_matrix,
                             plane_vector_bottom,
                             plane_point_bottom);

        std::vector<implicit_function *> children;
        children.push_back(first_object);
        children.push_back(hp_bottom);

        object = new implicit_functions::transformed_subtract(children, matrix12);


        register_new_object(object);

    } else if (name == "screw") {

        // inf screw subtract top bottom lid
        std::cout << "new _screw " << std::endl;

        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        Matrix<REAL, 3, 4> transformation_matrix;
        REAL pitch = 0.0;
        std::string profile;
        std::string end_type;
        REAL delta_ratio = 0.0;
        Matrix<REAL, 3, 1> v;

        implicit_functions::screw::getScrewParameters(
                           transformation_matrix, pitch, profile,
                           end_type, delta_ratio, v,
                           shapeparams_dict);


        transformation_matrix << 1, 0, 0, 0,
                                 0, 1, 0, 0,
                                 0, 0, 1, 0;

        // // the following code is taken from difference

        implicit_function * a = NULL;
        implicit_function * b = NULL;
        a = new implicit_functions::screw(transformation_matrix,
                                                         pitch,
                                                         profile,
                                                         end_type,
                                                         delta_ratio,
                                                         v);

        b = new top_bottom_lid(transformation_matrix);

        std::vector<implicit_function *> children;
        children.push_back(a);
        children.push_back(b);

        object = new implicit_functions::transformed_subtract(children, matrix12);

        register_new_object(object);

    } else if (name == "half_plane") {

        Eigen::Matrix<REAL, 3, 4> matrix;
        Eigen::Matrix<REAL, 3, 1> plane_vector;
        Eigen::Matrix<REAL, 3, 1> plane_point;

        half_plane::getHalfPlaneParameters(
            matrix,
            plane_vector,
            plane_point,
            shapeparams_dict
        );

        // plane_vector << 0, 0, 1;
        // plane_point << 0, 0, 0.5;
        // half_plane::getHalfPlaneParametersMatrixOnly(
        //     matrix,
        //     shapeparams_dict
        // );

        if(ignore_root_matrix) {
            matrix << 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0;
        };

        std::cout << "--------------ignore_root_matrix------------" << std::endl;
        std::cout << ignore_root_matrix << std::endl;

        std::cout << "--------------transformation_matrix------------" << std::endl;
        std::cout << matrix << std::endl;

        object = new half_plane(matrix,
                                 plane_vector,
                                 plane_point
                                 );

        register_new_object(object);

    } else if (name == "screw_gradient_wrong"){

        std::cout << "screw_tbb\n";

        Eigen::Matrix<REAL, 3, 4> eye;
        eye << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0;


        Matrix<REAL, 3, 4> transformation_matrix;
        REAL pitch = 0.0;
        std::string profile;
        std::string end_type;
        REAL delta_ratio = 0.0;
        Matrix<REAL, 3, 1> v;

        implicit_functions::screw::getScrewParameters(
                           transformation_matrix, pitch, profile,
                           end_type, delta_ratio, v,
                           shapeparams_dict);


        // implicit_functions::screw s  = new implicit_functions::screw(transformation_matrix,
        //                                  pitch,
        //                                  profile,
        //                                  end_type,
        //                                  delta_ratio,
        //                                  v);

//        convex_polygon* square = new convex_polygon(corners_x, corners_y);
//        unique_ptr<implicit_function_2d> polygon {square};


       implicit_function* s = new implicit_functions::screw(transformation_matrix,
                                         pitch,
                                         profile,
                                         end_type,
                                         delta_ratio,
                                         v);
       unique_ptr<implicit_function> screw {s};

        object = new implicit_functions::inf_top_bot_bound(transformation_matrix, screw);
        register_new_object(object);

    } else if (name == "top_bottom_lid") {

        Eigen::Matrix<REAL, 3, 4> matrix;
        Eigen::Matrix<REAL, 3, 1> plane_vector;
        Eigen::Matrix<REAL, 3, 1> plane_point;

        top_bottom_lid::getMatrix(
            matrix,
            shapeparams_dict
        );

        if(ignore_root_matrix) {
            matrix << 1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0;
        };

        std::cout << "--------------ignore_root_matrix------------" << std::endl;
        std::cout << ignore_root_matrix << std::endl;

        std::cout << "--------------transformation_matrix------------" << std::endl;
        std::cout << matrix << std::endl;

        object = new top_bottom_lid(matrix);

        register_new_object(object);

    } else if (name == "asmjscb") {

        REAL param1 = shapeparams_dict.get<REAL>("param1", 0.01);
        //int id = shapeparams_dict.get<int>("id", 94);
        int id = shapeparams_dict.get<int>("index", 94);

        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);

        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }
        object = new implicit_functions::javascript_implicit_function(id, matrix12, param1);
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

                //a = new implicit_functions::CrispUnion(*a, *b);
                // register_new_object(a);
                //std::vector<const implicit_function*> versus std::vector<implicit_function*>
                std::vector<implicit_function*> ab = std::vector<implicit_function*>();
                ab.push_back(a);
                ab.push_back(b);
                o_matrix = new implicit_functions::transformed_union(ab, matrix12);
                register_new_object(o_matrix);
                REAL eye_matrix12[12] = {1,0,0,0,  0,1,0,0,  0,0,1,0};
                o_plain =  new implicit_functions::transformed_union(ab, eye_matrix12);
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

        // object = new implicit_functions::CrispIntersection(*a, *b);
        //object = new implicit_functions::transformed_intersection(*a, *b, matrix12);
        std::vector<implicit_function *> children;
        children.push_back(a);
        children.push_back(b);
        object = new implicit_functions::transformed_intersection(children, matrix12);
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
                std::clog << "Error: Difference cannot be applied to more than two objects." << std::endl;
            }

        }

        // object = new implicit_functions::CrispSubtract(*a, *b);
        // object = new implicit_functions::transformed_subtract(*a, *b, matrix12);
        std::vector<implicit_function *> children;
        children.push_back(a);
        children.push_back(b);
        object = new implicit_functions::transformed_subtract(children, matrix12);
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
        object = new implicit_functions::meta_ball_Rydgård(matrix12, num_blobs, time, scale);

        register_new_object(object);
    } else if(name == "extrusion") {
        
        REAL matrix12[12];
        getMatrix12(matrix12, shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }

        Matrix<REAL, 3, 4> transformation_matrix;
        int size = 0;

        implicit_functions::extrusion::getExtrusionParameters(
                           size, shapeparams_dict);


        transformation_matrix << 1, 0, 0, 0,
                                 0, 1, 0, 0,
                                 0, 0, 1, 0;

        // // the following code is taken from difference

        implicit_function * a = NULL;
        implicit_function * b = NULL;

        REAL eye[12];
        copy_eye(eye);

        a = new implicit_functions::extrusion(eye, size);

        b = new top_bottom_lid(transformation_matrix);

        std::vector<implicit_function *> children;
        children.push_back(a);
        children.push_back(b);

        object = new implicit_functions::transformed_subtract(children, matrix12);

        register_new_object(object);

        /*
        //polygon is a square for now
        std::vector<REAL> corners_x = {0.1, -0.5, -0.5, 0.5};
        std::vector<REAL> corners_y = {0.1, 0.5, -0.5, -0.5};
        convex_polygon* square = new convex_polygon(corners_x, corners_y);
        unique_ptr<implicit_function_2d> polygon {square};
        */
        /*
        // infinite extrusion

        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);
        if(ignore_root_matrix) {
            copy_eye(matrix12);
        }
         object = new implicit_functions::extrusion(matrix12, 12);
        // object = new implicit_functions::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
        register_new_object(object);*/
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

}  // namespace mp5_implicit
