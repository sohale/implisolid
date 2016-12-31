#pragma once

namespace mp5_implicit{

namespace prtree = boost::property_tree;

/* Holds all MC information that is transferred through json. Separate from the shape json (i.e. the MP5 Json). */
struct mc_settings {
    mp5_implicit::bounding_box  box;
    dim_t resolution;
    bool ignore_root_matrix;
    int overall_repeats;

    struct {
        int iters;
        REAL c;
    } vresampl;
    struct {
        bool enabled;
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

    public:
    void report(ostream& os) {
        os <<   "ignore_root_matrix: " << ignore_root_matrix;
        os << ", box: ("
            << box.xmin << ","
            << box.ymin << ","
            << box.zmin << ") -- ("
            << box.xmax << ","
            << box.ymax << ","
            << box.zmax <<  ")";
        os << ", resolution: " << resolution;
        os << ", vresampl.iters: " << vresampl.iters;
        os << ", vresampl.c: " << vresampl.c;
        os << ", projection.enabled: " << projection.enabled;
        os << ", qem.enabled: " << qem.enabled;
        os << ", subdiv.enabled: " << subdiv.enabled;
        os << ", overall_repeats: " << overall_repeats;
        os << ", debug.post_subdiv_noise: " << debug.post_subdiv_noise;
        os << ".";

        if (
            this->qem.enabled &&
            !this->projection.enabled
        ) {
            os << "(Invalid settings: QEM will not be applied if centroid projection is disabled)";
        }
        if (
            !this->qem.enabled &&
            this->projection.enabled
        ) {
            os << "(Warning: projection without QEM has no visible effect. Vertices remain unchanged.)";
        }

        os << std::endl;
    }
    static const mc_settings default_settings() {
        mc_settings DEFAULT_SETTINGS;
        DEFAULT_SETTINGS.ignore_root_matrix = false;
        DEFAULT_SETTINGS.resolution = 28;
        // default: dont do resampling
        DEFAULT_SETTINGS.vresampl.iters = 0;
        DEFAULT_SETTINGS.vresampl.c = 1.0;
        DEFAULT_SETTINGS.projection.enabled = false;
        DEFAULT_SETTINGS.qem.enabled = false;
        DEFAULT_SETTINGS.overall_repeats = 1;
        DEFAULT_SETTINGS.debug.post_subdiv_noise = 0.01;
        DEFAULT_SETTINGS.subdiv.enabled = true;
        // DEFAULT_SETTINGS.box.xmin = NaN;
        return DEFAULT_SETTINGS;
    }
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
    std::string & abort_reason,
    std::vector<std::string> forbidden_names
        = std::vector<std::string>()
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
            abort_reason += "  field name instead of another one";
        }
    }
    return fieldval;
}

/*void print_mc_settings() {

}
*/
const mp5_implicit::mc_settings parse_mc_properties_json(const std::string & mc_parameters_json) {

    mp5_implicit::mc_settings DEFAULT_SETTINGS = mp5_implicit::mc_settings::default_settings();

    bool needs_abort = false;
    std::string abort_reason = "";

    // MC settings:
    mp5_implicit::mc_settings  mc_settings_from_json;  // settings

    // TODO(charles): find an alternativ to catch exceptions prtree::json_parser::json_parser_error prtree::ptree_bad_path
    // try{
    prtree::ptree mcparams_dict;
    std::stringstream mc_json_stream;
    mc_json_stream << mc_parameters_json;
    prtree::read_json(mc_json_stream, mcparams_dict);

    REAL xmin = mcparams_dict.get<REAL>("box.xmin", NaN);
    REAL xmax = mcparams_dict.get<REAL>("box.xmax", NaN);
    REAL ymin = mcparams_dict.get<REAL>("box.ymin", NaN);
    REAL ymax = mcparams_dict.get<REAL>("box.ymax", NaN);
    REAL zmin = mcparams_dict.get<REAL>("box.zmin", NaN);
    REAL zmax = mcparams_dict.get<REAL>("box.zmax", NaN);
    if ( isNaN(xmin) || isNaN(xmax) || isNaN(ymin) || isNaN(ymax) || isNaN(zmin) || isNaN(zmax)) {
        std::cerr << "Error: missing or incorrect values in mc_parameters_json"<< std::endl;
        xmin = -1;
        xmax = 1;
        ymin = -1;
        ymax = 1;
        zmin = -1;
        zmax = 1;

        needs_abort = true;
        abort_reason += "  err2";

    }
    mp5_implicit::bounding_box box = {xmin, xmax, ymin, ymax, zmin, zmax};  // {15,20,15,20,15,20};
    mc_settings_from_json.box = box;

    {
        REAL resolution_real = mcparams_dict.get<REAL>("resolution", -1);
        int resolution_int = static_cast<int>(resolution_real);
        if (resolution_int == -1) {
            resolution_int = DEFAULT_SETTINGS.resolution;
        }
        if (static_cast<REAL>(resolution_int) != resolution_real) {
            cerr << "Error: resolution must be integer: " << static_cast<REAL>(resolution_int) << " != " << resolution_real << std::endl;
            needs_abort = true;
            abort_reason += "  err3";

        }
        mc_settings_from_json.resolution = resolution_int;
    }
    if ( mc_settings_from_json.resolution <= 2 ) {
        std::cerr << "Error: missing or incorrect values in mc_parameters_json"<< std::endl;
        // resolution = 28;
        needs_abort = true;
        abort_reason += "  err4";
    }





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

    mc_settings_from_json.vresampl.c = mcparams_dict.get<REAL>(
        "vresampl.c", DEFAULT_SETTINGS.vresampl.c);
    // int default_vresampl_iters = (c == 0.0)? 0 : 1;  // if c is zero, ...
    mc_settings_from_json.vresampl.iters = mcparams_dict.get<int>(
        "vresampl.iters", DEFAULT_SETTINGS.vresampl.iters);


    // handling the default cases is frustrating. Also it's very difficult to use bool and int.

    mc_settings_from_json.projection.enabled = read_bool_from_json(mcparams_dict,
        "projection.enabled", DEFAULT_SETTINGS.projection.enabled, &needs_abort,
        abort_reason,
        std::vector<std::string>{"projection.enable"}
    );

    mc_settings_from_json.qem.enabled = read_bool_from_json(mcparams_dict,
        "qem.enabled", DEFAULT_SETTINGS.qem.enabled, &needs_abort,
        abort_reason,
        std::vector<std::string>{"qem.enable"}
    );

    // Used for correcting invalide bool numbers. Now we just use int.

    mc_settings_from_json.subdiv.enabled =
        read_bool_from_json(mcparams_dict,
            "subdiv.enabled",
            DEFAULT_SETTINGS.subdiv.enabled,
            &needs_abort,
            abort_reason,
            std::vector<std::string>{"subdiv.enable"}
        );


    mc_settings_from_json.debug.post_subdiv_noise = mcparams_dict.get<REAL>(
        "debug.post_subdiv_noise",
        DEFAULT_SETTINGS.debug.post_subdiv_noise);

    cout << "mc_settings_from_json.debug.post_subdiv_noise" << mc_settings_from_json.debug.post_subdiv_noise << std::endl;


    mc_settings_from_json.overall_repeats = mcparams_dict.get<int>(
        "overall_repeats", DEFAULT_SETTINGS.overall_repeats);


    // Shape settings
    mc_settings_from_json.ignore_root_matrix = mcparams_dict.get<bool>(
        "ignore_root_matrix", DEFAULT_SETTINGS.ignore_root_matrix);

    // Verification

    if (
        mc_settings_from_json.qem.enabled &&
        !mc_settings_from_json.projection.enabled
    ) {
        std::clog << "Invalid settings: QEM will not be applied if centroid projection is disabled" << std::endl;
        needs_abort = true;
        abort_reason += "  err7";
    }

    mc_settings_from_json.report(cout);

    if (needs_abort) {
        std::cout << "Aborting because of a problem in mc_settings_from_json. Abort reasons: "
            << abort_reason << std::endl;
        abort();
    }


    return mc_settings_from_json;
}

