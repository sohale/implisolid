#pragma once

/**
    A class that contains the internal variables in the loop of the main algorithm.
    The idea is to call the steps separately and externally (via web-worker, etc), while keeping the state.
    Note that _state and polygoniser__old are different classes.
    The class polygoniser__old was a failed attampt.
    
    Local variables in the algorithm are now member variables in this class.

    polygonizer is a state-machine, an automaton with instructions: init, ste0, step1, etc.

*/
struct polygonizer {

public:
// accessible
    std::string steps_report;
    timer timr;

private:
// shared states
    state_t & _state;  // shared state: output
    const mp5_implicit::implicit_function & object;  // shared object: input

// copy (input argument)
    const mp5_implicit::mc_settings mc_settings_from_json;
    const bool report_back_enabled = false;
    const worker_call_sepcs_t worker_call_sepcs;

public:
    polygonizer(
        state_t & _state,
        mp5_implicit::implicit_function& object,
        //const std::string & mc_parameters_json
        const mp5_implicit::mc_settings mc_settings_from_json_,
        const worker_call_sepcs_t &  worker_call_sepcs
    )
    :
        _state(_state),  // Does this assign a reference? or copies it into the referenced?
        object(object),
        // mc_settings_from_json(parse_mc_properties_json(mc_parameters_json)),
        mc_settings_from_json(mc_settings_from_json_),
        steps_report(""),
        timr(),
        report_back_enabled(true),
        worker_call_sepcs(worker_call_sepcs)
    {
        timr.report_and_continue("timer started.");

        // this->report_back_enabled = true;
    }
/*
shape_parameters_json, const char* mc_parameters_json
*/
    //static void polygonize_init(state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr,*/ bool use_metaball);
    void polygonize_step_0(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr,*/ /*bool use_metaball*/);
    void polygonize_step_1(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr*/);
    void polygonize_step_2(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr*/);
    void polygonize_step_3(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr,*/ bool is_last);

    void polygonize_terminate(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr*/);
    void send_mesh_back_to_client();
};


/*
void polygonizer::polygonize_init(state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, bool use_metaball, std::string& steps_report, timer & timr) {
    ;
}
*/

// todo: move into _state
void polygonizer::polygonize_step_0(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr, */ /*bool use_metaball*/) {

    auto vertsfaces_pair = mc_start(&object, mc_settings_from_json.resolution, mc_settings_from_json.box /*, use_metaball*/);
    // std::vector<REAL>, std::vector<int>
    steps_report = steps_report + "MC. ";

    /*
    TEST a SQUARE
    //auto
    vertsfaces_pair = make_a_square(10.0);
    */

    // auto  _state.mc_result_verts = _state.mc -> result_verts;
    // auto  _state.mc_result_faces = _state.mc->result_faces;
    _state.mc_result_verts = std::move(vertsfaces_pair.first);
    _state.mc_result_faces = std::move(vertsfaces_pair.second);

    timr.report_and_continue("marching cubes");
}

void polygonizer::polygonize_step_1(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr*/) {
    timer t1;

    const REAL c =  mc_settings_from_json.vresampl.c;  // 1.0;

    // result_verts is modified
    apply_vertex_resampling_to_MC_buffers__VMS(object, c, _state.mc_result_verts, _state.mc_result_faces, false );
    steps_report = steps_report + "V ";
    t1.stop("vertex resampling");  // 400 -> 200 -> 52 msec  (40--70)

    timr.report_and_continue("vertex resampling");
}


void polygonizer::polygonize_step_2(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr*/) {

    std::clog << "centroids_projection:" << std::endl;
    // Never send mc_settings_from_json as an argument
    centroids_projection(&object, _state.mc_result_verts, _state.mc_result_faces, mc_settings_from_json.qem.enabled);
    steps_report = steps_report + "P ";

    timr.report_and_continue("centroids_projection");
}


void polygonizer::polygonize_step_3(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr, */ bool is_last) {

    // Incorrect logic:
    const REAL scale_noise_according_to_matrix = 1.0 * 10.0;
    /*
    if (ignore_root_matrix) {
        scale_noise_according_to_matrix *= 1.0 / 10.0;
    } else {
        scale_noise_according_to_matrix *= 1.0;
    }
    */

    // For a realistic noise:
    // get_actual_matrix (even if ignored, dont ignore this one)
    // invert-it

    // Better logic: pass this as an argument:
    // object->get_noise_generator(matrix, ignore);

    const REAL actual_noise = is_last ?
            mc_settings_from_json.debug.post_subdiv_noise * scale_noise_according_to_matrix
        :
            0;
    if (!is_last) {
        std::clog << "noise skipped." << std::endl;
    } else {
        std::clog << "noise applied." << std::endl;
    }

    cout << "mc_settings_from_json.debug.post_subdiv_noise" << mc_settings_from_json.debug.post_subdiv_noise
        << "  actual noise: " << actual_noise << std::endl;
    my_subdiv_ ( _state.mc_result_verts, _state.mc_result_faces,
        actual_noise);
    std::clog << "outisde my_subdiv_." << std::endl;
    steps_report = steps_report + "Subdiv("+std::to_string(actual_noise)+")";

    timr.report_and_continue("subdivisions");
    
}
/*
    stores the result into _state.
    _state holds the result, and it will be available for later reference.
    Results are left in _state, available in _state.
*/
void polygonizer::polygonize_terminate(/*state_t & _state, const mp5_implicit::implicit_function& object, const mp5_implicit::mc_settings & mc_settings_from_json, std::string& steps_report, timer & timr*/) {
    //delete object;
    //object = NULL;
    gc_objects();

    _state.active = true;

    _state.check_state();
    // std::clog << "MC:: v,f: " << _state.mc_result_verts.size() << " " << _state.mc_result_faces.size() << std::endl;

    steps_report = steps_report + "finished.";
    std::clog << "algorithm.steps_report: " << steps_report << std::endl;

}

#include <emscripten.h>

void polygonizer::send_mesh_back_to_client() {
    if (!this->report_back_enabled) return;
    const worker_call_sepcs_t& wcs = this->worker_call_sepcs;
    // _state.mc_result_verts, _state.mc_result_faces
    int result = EM_ASM_INT(
        {
            /* This Javascript code is executed by browser (js engine) on WebWorker side */

            console.log(
                'Polyg. progress callback: verts (ptr, count):', $0, $1, " faces: ", $2, $3,
                " progr callback-id:", $4, "shape_id:", $5, " call_id", $6
            );

            /* `wwapi` is an object on WebWorker space */
            // global
            if (typeof wwapi === "undefined") return;

            // Translate those integer numbers to JS objects (verts_a, faces_a) in the scope of WW js engine.
            const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT; // 4
            const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT; // 4
            var nverts = Module. _get_v_size();
            var nfaces = Module. _get_f_size();
            assert(nverts>=0);
            assert(nfaces>=0);
            var verts_address = Module. _get_v_ptr();
            var faces_address = Module. _get_f_ptr();
            var verts_a = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
            var faces_a = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);
            const _progress_update_callback_id = $4; const _shape_id = $5; const _call_id = $6;
            wwapi.send_progress_update(verts_a, faces_a, _progress_update_callback_id, _shape_id, _call_id);

            return 3;
        },
        _state.mc_result_verts.data(),  // $0
        _state.mc_result_verts.size()/3,  // $1
        _state.mc_result_faces.data(),  // $2
        _state.mc_result_faces.size()/3,  // $3
        wcs.progress_callback_id,  // $4
        wcs.shape_id,  // $5
        wcs.call_id  // $6
    );

}

