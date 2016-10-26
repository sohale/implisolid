#pragma once

namespace mp5_implicit{

namespace prtree = boost::property_tree;

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
    struct {
        bool enabled;
    } subdiv;


    struct {
        // not used:
        // bool enabled_pointsets
        REAL post_subdiv_noise;
    } debug;

};

}

/**
 * @brief      Reads a bool from a json. If can tolerate integers as bool.
 * Without this function, using <bool>, it must be true/false, otherwise, it will fail.
 *
 * @param[in]  mcparams_dict    The json dictionary
 * @param[in]  fieldname        The fieldname, e.g. "qem.enabled"
 * @param[in]  default_value    The default value, e.g. false
 * @param      needs_abort      The needs abort. A bool which is set to true if we cannot continue and an exception should be raised or the program abort()ed.
 * @param[in]  forbidden_names  The forbidden field names, the firlds that are likely to be mistaken with our fiel, or the deprecated field.
 *
 * @return     { a bool }
 */
bool read_bool_from_json(
    const prtree::ptree& mcparams_dict,
    std::string fieldname,
    bool default_value,
    bool *needs_abort,
    std::vector<std::string> forbidden_names
        = std::vector<std::string>()  //
) {
    // default is false. dont use boolean in Json for *.enabled

    // todo: Also try to interpret as <bool> to begin with.
    // This goes to the default when I use "false" or "true"?
    int MAGICAL_NUMBER_DEFAULT = -1;
    int field_int = mcparams_dict.get<int>(fieldname, MAGICAL_NUMBER_DEFAULT);
    if (VERBOSE)
        clog << "json: " << fieldname << " : " << field_int << std::endl;

    bool fieldval = !!field_int;

    // If you write 'false' or 'true' it is an error, in both cases it will be false.
    // to handle the default cases seaprately. "true", "false" in Json both mean error.
    if (field_int == MAGICAL_NUMBER_DEFAULT) {
        fieldval = default_value;
        if (VERBOSE) {
            clog << "Value not specified: " << fieldname << " : Using " << fieldval << std::endl;
        }
    }

    // Used for correcting invalide bool numbers. Now we just use int.
    for (std::string  invalid_name : forbidden_names) {
        int MAGICALNUMBER = 9123456;
        // Check if the wrong spelling is used.
        if(mcparams_dict.get<int>(invalid_name, MAGICALNUMBER) != MAGICALNUMBER) {
            std::cerr << "Error: Use "<< fieldname << " instead of using" << invalid_name << ". Aborting! " << std::endl;
            *needs_abort = true;
        }
    }
    return fieldval;
}

/*void print_mc_settings() {

}
*/
mp5_implicit::mc_settings parse_mc_properties_json(const char* mc_parameters_json) {
    std::stringstream mc_json_stream;
    mc_json_stream << mc_parameters_json;

    prtree::ptree mcparams_dict;


    bool needs_abort = false;

    // MC settings:

    // TODO(charles): find an alternativ to catch exceptions prtree::json_parser::json_parser_error prtree::ptree_bad_path
    // try{
    prtree::read_json(mc_json_stream, mcparams_dict);

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


    /*}catch(prtree::json_parser::json_parser_error parse_exception){
        std::clog << "parse_exception"<< std::endl ;

    }catch(prtree::ptree_bad_data bad_data_exception){
        std::clog << "bad_data_exception"<< std::endl ;

    }catch(prtree::ptree_bad_path bad_path_exception){1
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

    /*
    // default is false. dont use boolean in Json for *.enabled
    int qem_enabled_int = mcparams_dict.get<int>("qem.enabled", -1);
    if (VERBOSE)
        clog << "json: qem.enabled : " << qem_enabled_int << std::endl;
    bool qem_enabled = !!qem_enabled_int;
    if (qem_enabled_int == -1) {   // If you write 'false' or 'true' it is an error, in both cases it will be false.
        qem_enabled = false;
    }
    */

    bool projection_enabled = read_bool_from_json(mcparams_dict, "projection.enabled", false, &needs_abort, std::vector<std::string>{"projection.enable"});
    bool qem_enabled = read_bool_from_json(mcparams_dict, "qem.enabled", false, &needs_abort, std::vector<std::string>{"qem.enable"});

    // Used for correcting invalide bool numbers. Now we just use int.

    /*
    int MAGICALNUMBER = 9123456;
    // Check if the wrong spelling is used.
    if(mcparams_dict.get<int>("projection.enable", MAGICALNUMBER) != MAGICALNUMBER) {
        std::cerr << "Error: Use projection.enabled instead of use projection.enable" << std::endl;
        needs_abort = true;
    }
    if(mcparams_dict.get<int>("qem.enable", MAGICALNUMBER) != MAGICALNUMBER) {
        std::cerr << "Error: Use qem.enabled instead of use qem.enable" << std::endl;
        needs_abort = true;
    }
    */


    bool subdiv_enabled = read_bool_from_json(mcparams_dict,
        "subdiv.enabled",
        true,
        &needs_abort, std::vector<std::string>{"subdiv.enable"}
    );

    /*
    // todo: make this a macro
    // todo: interpret as <bool>
    int subdiv_enabled_int = mcparams_dict.get<int>("projection.enabled", -1);  // Why this goes to the default when I use "false" or "true"?
    if (VERBOSE)
        clog << "subdiv_enabled : " << subdiv_enabled_int << std::endl;
    bool subdiv_enabled = !!subdiv_enabled_int;
    if (subdiv_enabled_int == -1) {
        subdiv_enabled = false;
    }
    if(mcparams_dict.get<int>("subdiv.enabled", MAGICALNUMBER) != MAGICALNUMBER) {
        std::cerr << "Error: Use subdiv.enabled instead of use subdiv.enable" << std::endl;
        needs_abort = true;
    }
    */


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

    mc_settings_from_json.subdiv.enabled = subdiv_enabled;

    REAL DEFAULT_SURFACE_NOISE = 0.01;
    mc_settings_from_json.debug.post_subdiv_noise = mcparams_dict.get<REAL>("debug.post_subdiv_noise", DEFAULT_SURFACE_NOISE);
    cout << "mc_settings_from_json.debug.post_subdiv_noise" << mc_settings_from_json.debug.post_subdiv_noise << std::endl;


    // Shape settings

    mc_settings_from_json.ignore_root_matrix = mcparams_dict.get<bool>("ignore_root_matrix", false);

    if (needs_abort) {
        abort();
    }


    return mc_settings_from_json;
}

