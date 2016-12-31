typedef int shapeid_t;

/**
  Maintains a single int. But can be used for more complicated call patterns in future.
  Used when multiple results are sent back to the client from one single request.
*/
struct worker_call_sepcs_t {
// protected:
    int progress_callback_id = -1;
    int call_id = -1;
    shapeid_t shape_id = -1;

 public:
    bool is_progressive() {
        return progress_callback_id != -1;
    }

    /**
       Specified the follow-up code when the result of the progress of polygonisation is sent back to the client or the non-worker side.
    */
    explicit  worker_call_sepcs_t (const std::string & callspecs_dict_json) {
        prtree::ptree callspecs_dict;
        std::stringstream mc_json_stream;
        mc_json_stream << callspecs_dict_json;
        prtree::read_json(mc_json_stream, callspecs_dict);

        // parse_worker_call_sepcs_t specs;
        this->progress_callback_id =
            callspecs_dict.get<int>("progressCallback_id", -1);
        this->call_id =
            callspecs_dict.get<int>("call_id", -1);
        this->shape_id =
            (shapeid_t)(callspecs_dict.get<int>("shape_id", -1));
    }

    /**
    When this constructor used, there is no "progressive geometry update", in which ther are multiple postMessage() from a single worker task.
    todo: Still the call_id and shape_id might be relevant
    This contruvtoir is used when progress_update_specs (i.e. worker_call_specs) is undefined on the js side.
    */
    explicit  worker_call_sepcs_t () {
        this->progress_callback_id = -1;
        this->call_id = -1;
        this->shape_id = -1;
    }

};

/**
    progress_update_specs is either undefined or {progressCallback_id: .., call_id:.., shape_id: ...} on JS side.
*/

/*
const worker_call_specs_t  parse_worker_call_sepcs (
    const std::string & callspecs_dict_json
) {
    return worker_call_sepcs_t(callspecs_dict_json);
}
*/

