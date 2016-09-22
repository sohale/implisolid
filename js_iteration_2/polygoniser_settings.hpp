#pragma once

namespace mp5_implicit{

/* Holds all MC information that is transferred through json. Separate from the shape json (i.e. the MP5 Json). */
struct mc_settings {
    mp5_implicit::bounding_box box;
    dim_t resolution;
    bool ignore_root_matrix;
    struct {
        int iters;
        REAL c;
    } vresampl;
    struct {
        bool enabled = true;
    } qem;
    struct {
        // projection of centroids, not vertices
        bool enabled;
    } projection;  // cproj
};

}

/*void print_mc_settings() {

}
*/
mp5_implicit::mc_settings parse_mc_properties_json(const char* mc_parameters_json) {
    std::stringstream mc_json_stream;
    mc_json_stream << mc_parameters_json;

    namespace pt = boost::property_tree;
    pt::ptree mcparams_dict;


    bool needs_abort = false;

    // MC settings:

    // TODO(charles): find an alternativ to catch exceptions pt::json_parser::json_parser_error pt::ptree_bad_path
    // try{
    pt::read_json(mc_json_stream, mcparams_dict);

    REAL xmin = mcparams_dict.get<REAL>("box.xmin", NaN);
    REAL xmax = mcparams_dict.get<REAL>("box.xmax", NaN);
    REAL ymin = mcparams_dict.get<REAL>("box.ymin", NaN);
    REAL ymax = mcparams_dict.get<REAL>("box.ymax", NaN);
    REAL zmin = mcparams_dict.get<REAL>("box.zmin", NaN);
    REAL zmax = mcparams_dict.get<REAL>("box.zmax", NaN);


    REAL resolution_real = mcparams_dict.get<REAL>("resolution", -1);
    int resolution = static_cast<int>(resolution_real);
    if (static_cast<REAL>(resolution) != resolution_real) {
        cerr << "Error: resolution must be integer: " << static_cast<REAL>(resolution) << " != " << resolution_real << std::endl;
        needs_abort = true;
    }


    if ( isNaN(xmin) || isNaN(xmax) || isNaN(ymin) || isNaN(ymax) || isNaN(zmin) || isNaN(zmax) || resolution <= 2 ) {
        std::cerr << "Error: missing or incorrect values in mc_parameters_json"<< std::endl;
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;
        resolution = 28;

        needs_abort = true;
    }
    // std::clog << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << " " << resolution << " " << std::endl;


    /*}catch(pt::json_parser::json_parser_error parse_exception){
        std::clog << "parse_exception"<< std::endl ;

    }catch(pt::ptree_bad_data bad_data_exception){
        std::clog << "bad_data_exception"<< std::endl ;

    }catch(pt::ptree_bad_path bad_path_exception){1
        std::clog << "bad_path_exception" << std::endl ;

    }catch(...){
        std::clog << "other_exception" << std::endl ;
    }*/


    // post-MC steps (I, II, III):


    // bool apply_projection = false;
    // bool apply_qem = false;
    // REAL qem_reject_maxlen = ...;
    // REAL qem_tau = ...;
    // REAL projection_maxlen = ...;
    // mc_settings_from_json.vertexresampleing.c
    // mc_settings_from_json.vertexresampleing.iterations

    // default: dont do resampling
    REAL c = mcparams_dict.get<REAL>("vresampl.c", 1.0);
    // int default_vresampl_iters = (c == 0.0)? 0 : 1;  // if c is zero, ...
    int default_vresampl_iters = 0;
    int vresamp_iters = mcparams_dict.get<int>("vresampl.iters", default_vresampl_iters);

    // handling the default cases is frustrating. Also it's very difficult to use bool and int.

    // default is false. dont use boolean in Json for *.enabled
    int qem_enabled_int = mcparams_dict.get<int>("qem.enabled", -1);
    if (VERBOSE)
        clog << "json: qem.enabled : " << qem_enabled_int << std::endl;
    bool qem_enabled = !!qem_enabled_int;
    if (qem_enabled_int == -1)   // If you write 'false' or 'true' it is an error, in both cases it will be false.
        qem_enabled = false;
    // todo: interpret as <bool>
    int projection_enabled_int = mcparams_dict.get<int>("projection.enabled", -1);  // Why this goes to the default when I use "false" or "true"?
    if (VERBOSE)
        clog << "projection_enabled : " << projection_enabled_int << std::endl;
    bool projection_enabled = !!projection_enabled_int;
    if (projection_enabled_int == -1)  // to handle the default cases seaprately. "true", "false" in Json both mean error.
        projection_enabled = false;

    // Check if the wrong spelling is used.
    if(mcparams_dict.get<int>("projection.enable", 9123456) != 9123456) {
        std::cerr << "Error: Use projection.enabled instead of use projection.enable" << std::endl;
        needs_abort = true;
    }
    if(mcparams_dict.get<int>("qem.enable", 9123456) != 9123456) {
        std::cerr << "Error: Use qem.enabled instead of use qem.enable" << std::endl;
        needs_abort = true;
    }

    if (needs_abort) {
        abort();
    }


    mp5_implicit::mc_settings  mc_settings_from_json;  // settings
    mp5_implicit::bounding_box box = {xmin, xmax, ymin, ymax, zmin, zmax};  // {15,20,15,20,15,20};
    mc_settings_from_json.box = box;
    mc_settings_from_json.resolution = resolution;

    mc_settings_from_json.vresampl.c = c;
    mc_settings_from_json.vresampl.iters = vresamp_iters;

    mc_settings_from_json.qem.enabled = qem_enabled;
    mc_settings_from_json.projection.enabled = projection_enabled;

    if (mc_settings_from_json.qem.enabled && !mc_settings_from_json.projection.enabled) {
        std::clog << "Warning: QEM will not be applied if centroid projection is disabled" << std::endl;
    }

    // Shape settings

    mc_settings_from_json.ignore_root_matrix = mcparams_dict.get<bool>("ignore_root_matrix", false);

    return mc_settings_from_json;
}

