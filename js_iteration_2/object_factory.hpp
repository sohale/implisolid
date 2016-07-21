
//#include "unit_sphere.hpp"
#include "primitives.cpp"
#include "crisp_subtract.hpp"
//using namespace mp5_implicit;


implicit_function*  object_factory_simple(REAL f_argument, std::string name){
    implicit_function* object;
    /*
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
    //else if (name == "scylinder"){
    //    object = new mp5_implicit::scylinder(f_argument, 1.6); //0.7
    //}
    //else
    */
    if(name == "meta_balls"){
        REAL r = (sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3)*1;
        std::cout << " META BALLS r : " << r << std::endl;
        object = new mp5_implicit::unit_sphere(r);
    }
    /*else if(name == "sub_spheres"){
        mp5_implicit::unit_sphere * s1 = new mp5_implicit::unit_sphere(2, 1, 1, 1);
        mp5_implicit::unit_sphere * s2 = new mp5_implicit::unit_sphere(1.3);
        object = new mp5_implicit::CrispSubtract(*s1, *s2);
    }
    */
    else {
        std::cout << "Error! You must enter a valid name " <<  name << "! So I made a sphere!" << std::endl;
        //object = new mp5_implicit::unit_sphere(sin(0.033*10 * f_argument * 3.1415*2.)*0.33+0.3);
        object = new mp5_implicit::unit_sphere(1000.);
    }
    return object;
}


namespace pt = boost::property_tree ;

void getMatrix12(REAL * matrix12, const pt::ptree& shapeparams_dict){
        int i = 0;
        for (const pt::ptree::value_type &element : shapeparams_dict.get_child("matrix")){

            REAL x = element.second.get_value<REAL>();
            //std::cout << "matrix value : " << x << std::endl;
            matrix12[i] = x;
            i++;
        }
}

implicit_function*  object_factory(pt::ptree shapeparams_dict, bool& use_metaball, bool ignore_root_matrix) {
    std::string name = shapeparams_dict.get<std::string>("type");
    //REAL xmax = shapeparams_dict.get<REAL>("matrix",NaN);
    //std::cout << "############Name : " << name << std::endl;
    //REAL zmax = shapeparams_dict.get<REAL>("box.zmax",NaN);
    // int resolution = shapeparams_dict.get<int>("resolution",-1);

    // if(isNaN(xmin) || isNaN(xmax) || isNaN(ymin) || isNaN(ymax) || isNaN(zmin) || isNaN(zmax) || resolution <= 2 ){
    //     std::cout << "Error: missing or incorrect values in mc_parameters_json"<< std::endl;
    //     xmin = -1;
    //     xmax = 1;
    //     ymin = -1;
    //     ymax = 1;
    //     zmin = -1;
    //     zmax = 1;
    //     resolution = 28;
    // }




    if(name=="meta_balls"){
        use_metaball = true;
    }
    else{
        use_metaball = false;
    }

    implicit_function* object;

    if (name == "implicit_double_mushroom"){
        std::cout << "implicit_double_mushroom case " << std::endl;
        // object = new mp5_implicit::double_mushroom(0.9, 0.4 ,0.4 , 1/0.2 );
        object = new mp5_implicit::double_mushroom(0.9, 0.4/2, 0.4/2, 1/0.2 );
    }
    else
    /*if (name == "simple_sphere"){
        //std::cout << "******************* simple_sphere case " << std::endl;
        REAL radius = shapeparams_dict.get<REAL>("radius");
        object = new mp5_implicit::unit_sphere(radius);
        std::cout << "radius " << radius << std::endl;
    }
    else*/
    if (name == "icube" || name == "cube" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        object = new mp5_implicit::cube(matrix12);
       // object = new mp5_implicit::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
    }
    else
    if (name == "icylinder" || name == "cylinder" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        object = new mp5_implicit::scylinder(matrix12);
       // object = new mp5_implicit::cube(f_argument+0.2, f_argument+0.2, f_argument+0.2);
    }
    else
    if (name == "iellipsoid" || name == "ellipsoid" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        object = new mp5_implicit::egg(matrix12);
    }else
    if (name == "icone" || name == "cone" ){
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        object = new mp5_implicit::scone(matrix12);
    }

    else if (name == "Union") {
        //todo: Use SimpleUnion if (matrix12 == eye(4))
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        implicit_function * a = NULL;

        implicit_function * o_matrix;
        implicit_function * o_plain;

        for (pt::ptree::value_type &element : shapeparams_dict.get_child("children")) {
            if (a == NULL){
                a = object_factory(element.second, use_metaball, false);
            } else {
                implicit_function * b = object_factory(element.second, use_metaball, false);

                //The following always prints an empty line:
                //std::cout << "element.second.get_value<string>(\"type\")" << element.second.get_value<string>("type") << std::endl ;

                //a = new mp5_implicit::CrispUnion(*a, *b);
                //std::vector<const implicit_function*> versus std::vector<implicit_function*>
                std::vector<implicit_function*> ab = std::vector<implicit_function*>();
                ab.push_back(a);
                ab.push_back(b);
                o_matrix = new mp5_implicit::transformed_union(ab, matrix12);
                REAL eye_matrix12[12] = {1,0,0,0,  0,1,0,0,  0,0,1,0};
                o_plain =  new mp5_implicit::transformed_union(ab, eye_matrix12);
                a = o_plain;
            }

            //std::cout << "##### " << element.second.get_child("type") << std::endl;
            //i++;
        }
        a = o_matrix;
        //std::cout  << "#####" << shapeparams_dict.get<string>("children") << std::endl;
        //implicit_function * a = object_factory();
        //implicit_function * b;

        object = a;
        //object = new mp5_implicit::egg(matrix12);
    }else if (name == "Intersection") {
        //todo: Use SimpleUnion if (matrix12 == eye(4))
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        implicit_function * a = NULL;
        implicit_function * b = NULL;
        int count = 0;
        for (pt::ptree::value_type &element : shapeparams_dict.get_child("children")) {
            if (a==NULL){
                a = object_factory(element.second, use_metaball, false);
            }else{
                if(count > 1){
                    std::cout << "An CrispIntersection should have only 2 child" << std::endl;
                    break;
                }
                b = object_factory(element.second, use_metaball, false);

            }
            count++;
        }

        object = new mp5_implicit::CrispIntersection(*a, *b);
    }else if (name == "Difference") {
        //todo: Use SimpleUnion if (matrix12 == eye(4))
        REAL matrix12[12];
        getMatrix12(matrix12,shapeparams_dict);

        implicit_function * a = NULL;
        implicit_function * b = NULL;
        int count = 0;
        for (pt::ptree::value_type &element : shapeparams_dict.get_child("children")) {
            if (a==NULL){
                a = object_factory(element.second, use_metaball, false);
            }else{
                if(count > 1){
                    std::cout << "An CrispSubstraction should have only 2 child" << std::endl;
                    break;
                }
                b = object_factory(element.second, use_metaball, false);

            }
            count++;
        }

        object = new mp5_implicit::CrispSubtract(*a, *b);
    }
    else{
        // if(name=="meta_balls")
        REAL f_argument = shapeparams_dict.get<REAL>("time", NaN);

        std::cout << "otherwise " << "you asked for " << name << std::endl;
        object = object_factory_simple(f_argument, name);
    }

    return object;

}




implicit_function*  object_factory(string shape_parameters_json, bool& use_metaball, bool ignore_root_matrix) {


    /*
    mp5_implicit :: unit_sphere   object(sin(0.033*10 * time * 3.1415*2.)*0.33+0.3);
    // //_state.mc -> prepare_grid(1.0);
    // //object.eval_implicit(grid, implicit_values);
    _state.mc -> eval_shape(object, 1.0);
    */


/*
    std::string name = std::string(obj_name);
    std::cout << "Name : " << name << std::endl;

    REAL f_argument = time;

    implicit_function* object = object_factory_simple(f_argument, name);
*/
    //std::cout << "############################" << shape_parameters_json << std::endl;
    std::stringstream shape_json_stream;
    shape_json_stream << shape_parameters_json ;


    pt::ptree shapeparams_dict;

    pt::read_json(shape_json_stream, shapeparams_dict);

    return object_factory(shapeparams_dict, use_metaball, ignore_root_matrix);

}