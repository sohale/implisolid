
'use strict';

/******************************************************
 * Low level: Level 1
 ******************************************************/


var wapi1 = {};  // todo: rename

wapi1.worker = new Worker('../js_iteration_2/js/worker_api.js');
wapi1.worker_response_callbacks_table = {};
const WORKER_DEBUG = true;
wapi1.w_api_debug_info = {};
wapi1.w_api_info2 = {};  // w_api_info2.READABLE_COMPILE_TIME_NAME = ...  // such as function name or a name that is similar to a funciton name
wapi1.call_counter = 0;  // a clobal counter that assigns a unique id to each individual call, to make it easy to trace correspondence of call and callback.

var wcodes1 = {};
wcodes1.__wapi_query_implicit_values = 4;
wcodes1.__wapi_get_latest_vf = 6;
wcodes1.__wapi_make_geometry = 5;
wcodes1.__wapi_make_geometry_progress_update = 7;

var callstamps = {};

wapi1.worker.addEventListener('message', function(event) {
    if (WORKER_DEBUG)
        console.info("Message received from worker:", event, "data=", event.data);

    var is_progress_update = event.data.is_progress_update;
    var isFinal = !is_progress_update;

    var _callbackId = event.data.return_callback_id;
    var mycallback = wapi1.worker_response_callbacks_table[ _callbackId ];
    
    if (mycallback) {
        if (isFinal) {
            // delete wapi1.worker_response_callbacks_table[ _callbackId ];
            // don’t delete array elements. It would make the array transition to a slower internal representation.
            wapi1.worker_response_callbacks_table[ _callbackId ] = undefined;  // per-call
        }
        if (WORKER_DEBUG) {
            // tracking calls, measuring the time, detecting orphan calls (calls with no return), etc.
            var call_id = event.data.call_id;
            var start_time = callstamps[call_id];
            console.info("TIME: ", new Date() - start_time , "msec");
        }
    }

    //if (_callbackId === 4) {
    //    console.info("The result of 'query_implicit_values()' was: ", event.data.returned_data); //event.data.result_allpositive
    //}
    if (!mycallback) {
        console.error("Returned data from worked did not match any registered callbackId: ", event.data, event);
        return;
    }
    if (WORKER_DEBUG)
        console.info("calling the callback ", mycallback.name);

    var call_identification = event.data.call_id;
    var shape_id = event.data.shape_id;
    //var is_progress_update = event.data.is_progress_update;
    if (is_progress_update) {  //bad
        // {return_callback_id:progres_update_callback_id, returned_data: result, call_id: call_id, shape_id:shape_id, is_progress_update: true}
        // actually same format is used (same side args):
        mycallback (event.data.returned_data, call_identification, shape_id);
    }

    mycallback(event.data.returned_data, call_identification, shape_id);
    // todo: mycallback(event.data.returned_data, call_identification, event.data.shape_id);

 });


function worker_call_preparation(_callbackId, result_callback) {
    wapi1.call_counter ++ ;
    if (WORKER_DEBUG) {
        var temp = wapi1.worker_response_callbacks_table[_callbackId];
        if (temp) {
            console.warn("callback already set.");
            if (!(temp === result_callback))
                console.error("Warning: and callback is different!: ", temp, "!==", result_callback);
        }
    }

    wapi1.worker_response_callbacks_table[_callbackId] = result_callback;
    if (WORKER_DEBUG) {
        wapi1.w_api_debug_info[_callbackId] = result_callback;
        // also you can store call timestamp, etc
    }


    //console.info(callstamps, wapi1.call_counter);
    callstamps[wapi1.call_counter] = new Date();
    //console.info(callstamps);

    return _callbackId;
}

/*******************************************************************************
 * API Level 2
*******************************************************************************/

// Level 2
var wapi2 = {};
/**
* instead of result_callback --> completion_callback, updated_callback.
*/
wapi2.wapi_make_geometry = function (shape_id0, mp5_str, polygonization_settings_json, result_callback, progress_update_callback) {
    // based on implisolid_main.js
    //var startTime = new Date();
    const usage_info = "Usage: wapi_make_geometry(shape_id, shape_json, polyg_json, callback_func, progress_callback_func)";

    var progressive_completion = !!progress_update_callback;
    assert(!progressive_completion || typeof progress_update_callback === 'function', "" + usage_info);

    my_assert(typeof result_callback === 'function', "A completion callback is not provided. " + usage_info);
    if (progressive_completion)
    my_assert(typeof progress_update_callback === 'function', "The specified progress callback is not a function. " + usage_info);


    var _callbackId = worker_call_preparation(wcodes1.__wapi_make_geometry, result_callback);
    var _progressCallbackId = undefined;
    if (progressive_completion)
        _progressCallbackId = worker_call_preparation(wcodes1.__wapi_make_geometry_progress_update, progress_update_callback);
    var wreq = {
            funcName: 'make_geometry', //progressive_completion? 'make_geometry' : 'make_geometry_u',
            callbackId: _callbackId,
            progressCallbackId: progressive_completion? _progressCallbackId : 0,
            call_id: wapi1.call_counter,
            // call_timestampe: new Date(),

            mp5_str: null, polygonization_settings: null, obj_req_id: 0
        };
    wreq.mp5_str = mp5_str;
    wreq.polygonization_settings = polygonization_settings_json;  // mc_params;
    wreq.obj_req_id = shape_id0; //wapi1.call_counter;
    wapi1.worker.postMessage(wreq);

    /*
    const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
    const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT

    if(impli1.needs_deallocation) {
        impli1.finish_geometry();
        impli1.needs_deallocation = false;
    }
    impli1.build_geometry(mp5_str, mc_params_json);
    impli1.needs_deallocation = true;

    var nverts = impli1.get_v_size();
    var nfaces = impli1.get_f_size();

    var verts_address = impli1.get_v_ptr();
    var faces_address = impli1.get_f_ptr();

    var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
    var faces = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);

    // first iteration: the callback
    var geom = result_callback(verts, faces);

    var endTime = new Date();
    var timeDiff = endTime - startTime;

    //report_time(timeDiff, function(){hist();});

    return geom;
    */

}

//     w_impli2.query_implicit_values_old = function(mp5_str, points, reduce_callback, call_id, shape_index)
/**
@param client_result_callback is needed for all worker API functions. client_result_callback is a fuction that receives the and called asynchroniously. Even if it sdoesn't have eany argument, it notifies completion f the worker's job.
*/
wapi2.wapi_query_implicit_values_old = function(shape_index, shape_json, points, reduce_type_str, epsilon, client_result_callback) {
    my_assert(typeof reduce_type_str === 'string');
    if(shape_index  === undefined ) shape_index = "test(undefiend)";

    // registers the function client_result_callback to listen to the messages posted from the webworker
    var _callbackId = worker_call_preparation(wcodes1.__wapi_query_implicit_values_old, client_result_callback);
    var wreq = {
            funcName: 'query_implicit_values_old',
            callbackId: _callbackId,
            call_id: wapi1.call_counter,

            mp5_str: null, points: null, reduce_type: "", epsilon: 0.0000,
            obj_req_id: -1
        };
    wreq.mp5_str = shape_json;
    wreq.points = points;  // the points are set to the worker via the postMessage, as the Float32Array they are.
    wreq.reduce_type = reduce_type_str;
    wreq.epsilon = epsilon;
    wreq.shape_id = shape_index;

    wapi1.worker.postMessage(wreq);
}

/** called by query_a_normal() 
@param: clientside_result_callback  = function((normals, call_id, shape_id)){output_vect3.set(-r[0], -r[1], -r[2]);}
*/
wapi2.wapi_query_a_normal = function() {
    //static
    var pointbuffer = new Float32Array(3);

    return function wapi_query_a_normal(shape_index, mp5_shape_json, point, clientside_result_callback) {
        my_assert(typeof mp5_shape_json === 'string');

        // registers the worker callback
        var _callbackId = worker_call_preparation(wcodes1.__wapi_query_a_normal, clientside_result_callback);

        if (point instanceof Float32Array) {
            pointbuffer.set(point);
        } else {
            pointbuffer[0] = point.x;
            pointbuffer[1] = point.y;
            pointbuffer[2] = point.z;
        }

        var wreq = {
                funcName: 'query_normals',
                callbackId: _callbackId,
                call_id: wapi1.call_counter,

                mp5_str: null, points: null,
                obj_req_id: -1
            };
        wreq.mp5_str = mp5_shape_json;
        wreq.points = pointbuffer;
        wreq.shape_id = shape_index;

        wapi1.worker.postMessage(wreq);
    };
}();



wapi2.wapi_get_latest_vf = function (shape_id, result_callback) {
    my_assert(typeof result_callback === 'function');

    var _callbackId = worker_call_preparation(wcodes1.__wapi_get_latest_vf, result_callback);
    var wreq = {
            funcName: 'get_latest_vf',
            callbackId: _callbackId,
            call_id: wapi1.call_counter,

            obj_req_id: 0
        };
    wreq.shape_id = shape_id, // wreq.obj_req_id = obj_id, //wapi1.call_counter;
    wapi1.worker.postMessage(wreq);
}


/*******************************************************************************
 * API Level 3
 ******************************************************************************/

var wapi3 = {};  // not worker side, but uses worker

// wapi3.custom_mc_settings = {};

// Copied from implisolid_main.js
var SET_ME_TO_TRUE = false;
/*
wapi3.getLiveGeometry  = function(shape_properties, bbox, ignore_root_matrix, geom_callback) {
    var mc_properties = IMPLICIT.make_polygonization_settings(bbox, ignore_root_matrix);

    var shape_json_str = JSON.stringify(shape_properties);
    var polygonization_properties_json = JSON.stringify(mc_properties);
    this.getLiveGeometry_from_json(shape_json_str, polygonization_properties_json, geom_callback);
}
*/
wapi3.getLiveGeometry_from_json  = function(shape_id00, shape_json_str, polygonization_setttings_json_str, geom_callback) {

    assert(typeof shape_id00 === "number");
    assert(typeof shape_json_str === "string");
    assert(typeof polygonization_setttings_json_str === "string");
    assert(typeof geom_callback === "function");

    // wapi2.make_geometry
    //var geom =   // no geom returned here anymore
    var shape_id0 = shape_id00;
    wapi2.wapi_make_geometry(shape_id0, shape_json_str, polygonization_setttings_json_str,
        function (vf_dict, call_id, shapeid_) {
            var verts = vf_dict.verts;
            var faces = vf_dict.faces;
            assert(shapeid_ == shape_id0);
            // ThreeJS-specific code

            // var ignore_root_matrix: Does not need other (MC-related) arguments.

            var allocate_buffers = true;
            var geom = new LiveBufferGeometry79(verts, faces, allocate_buffers);

            // Set the normals
            // var ignore_root_matrix = mc_params.ignore_root_matrix;  // Does not need other (MC-related) arguments.
            //geom.update_normals(this_, verts, shape_json_str, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.
            //
            if (!shape_json_str)
                console.error(shape_json_str);
            if (SET_ME_TO_TRUE)
            wapi3.make_normals_into_geometry(geom, shapeid_, shape_json_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

            //this_.aaaaaaaaaA(verts);

            // no geom anymore
            // return geom;
            //var shape_id = 99; // vf_dict.shape_id
            geom_callback(geom, shape_id0);
        }
    );

    // return geom;
};


wapi3.wapi_query_implicit_values = function (obj_id, mp5_str, points, result_callback) {
    var _callbackId = worker_call_preparation(wcodes1.__wapi_query_implicit_values, result_callback);

    var wreq = {
            funcName: 'query_implicit_values',
            callbackId: _callbackId,
            call_id: wapi1.call_counter,
            //args: {
            mp5_str: null, points: null, reduce_callback: "", obj_id: 0
            //},
        };
    wreq.mp5_str = mp5_str;
    wreq.points = points;
    wreq.reduce_callback = 'any-positive';  // all outside or on
    wreq.obj_id = obj_id;
    wapi1.worker.postMessage(wreq);
}





wapi3.receive_mesh = function(shape_id_, result_callback) {
    var vf = {faces: null, verts: null};
    /*
    var nonempty = wapi2.wapi_get_latest_vf(vf);
    if (nonempty) {
        // vf already contains the output of get_latest_vf()
    } else{
        console.log("empty implicit. Using a default shape.");
        vf['verts'] = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 ]);
        vf['faces'] = new Uint32Array([0,1,2, 0,2,3, 0,4,5, 0,5,1, 1,5,6, 1,6,2, 2,6,3, 3,6,7, 4,5,6, 5,6,7]);
    }
    assert(vf['verts']);
    assert(vf['faces']);
    */
    // should be no logic here

    wapi2.wapi_get_latest_vf(shape_id_, function c3(data_received_from_worker, call_id, shape_id){
        data_received_from_worker.verts;
        data_received_from_worker.faces;
        data_received_from_worker.nonempty;
        var vf = data_received_from_worker;
        var nonempty = data_received_from_worker.nonempty;

        if (nonempty) {
            // vf already contains the output of get_latest_vf()
        } else{
            console.log("empty implicit. Using a default shape.");
            vf['verts'] = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 ]);
            vf['faces'] = new Uint32Array([0,1,2, 0,2,3, 0,4,5, 0,5,1, 1,5,6, 1,6,2, 2,6,3, 3,6,7, 4,5,6, 5,6,7]);
        }
        assert(vf['verts']);
        assert(vf['faces']);
        //return vf;
        result_callback(shape_id_, vf);
    });
}


/**
 ver 1 is a test for a recursive - like function in side another function (function composition; promise chaining).
*/
wapi3.update_geometry_from_json_ver1 = function(geometry, shape_id, shape_json, polygonization_json, done_callback, _extra) {
    //var startTime = new Date();

    assert(_extra === undefined);
    if (typeof shape_json !== "string") {
        shape_json = JSON.stringify(shape_json);
    }
    assert(typeof done_callback !== "number" && typeof done_callback !== "boolean");
    assert(typeof done_callback === "function");

    // polygonization_json should not be jsonified.
    if (typeof polygonization_json !== "string") {
        polygonization_json = JSON.stringify(polygonization_json);
    }

    wapi2.wapi_make_geometry (shape_id, shape_json, polygonization_json, //result_callback,
        //function (verts, faces) {
        function c3(vf_result, call_id1, shape_id1) {
            //var allocate_buffer=false;
            vf_result.faces;
            vf_result.verts;
            // vf_result already contains it!!

            // todo: make it the same query. dont call get_vf again.
            wapi3.receive_mesh(shape_id1, function(shape_id__, vf) {
                vf.faces;
                vf.verts;
                vf.nonempty;
                
                // todo: match the object_id to the request call_id, etc

                var ignoreDefaultNormals = true;
                var bool_reallocated = geometry.update_geometry1(vf['verts'], vf['faces'], ignoreDefaultNormals, false);
                geometry.use_default_normals_from_vertices();

                /*
                var startTime = callstamps[call_id];
                var endTime = new Date();
                var timeDiff = endTime - startTime;
                console.log(nv3 + " , " + nf3);
                console.log("Time: "+timeDiff+ " msec.");
                report_time(timeDiff);
                */

                //geometry is now updated
                //var call_id = 10000; 
                //done_callback(null, call_id, shape_id__);
                done_callback(shape_id__, null);
            });
        }
    );

};

/**
ver 2 is faster (simpler), doing it in a single request. 
*/
wapi3.update_geometry_from_json_ver2 = function(geometry, shape_id, shape_json, polygonization_json, done_callback, _extra) {
    if (typeof shape_json !== "string") shape_json = JSON.stringify(shape_json);
    assert(_extra === undefined);
    assert(typeof done_callback !== "number" && typeof done_callback !== "boolean");
    assert(typeof done_callback === "function");
    if (typeof polygonization_json !== "string") polygonization_json = JSON.stringify(polygonization_json);

    wapi2.wapi_make_geometry (shape_id, shape_json, polygonization_json, //result_callback,
        function c3(vf_result, call_id1, shape_id1) {
            vf_result.faces; vf_result.verts;
            assert(shape_id1 === shape_id);
            var ignoreDefaultNormals = true;
            var bool_reallocated = geometry.update_geometry1(vf_result.verts, vf_result.faces, ignoreDefaultNormals, false);
            geometry.use_default_normals_from_vertices();
            done_callback(shape_id1, null);
        }
    );
};

/* based on update_geometry_from_json_ver2() */
wapi3.update_geometry_from_json_progressive = function(geometry, shape_id, shape_json, polygonization_json, done_callback, progress_update_callback) {
    if (typeof shape_json !== "string") shape_json = JSON.stringify(shape_json);
    assert(typeof done_callback !== "number" && typeof done_callback !== "boolean");
    assert(typeof done_callback === "function");
    assert(typeof progress_update_callback === "function");
    if (typeof polygonization_json !== "string") polygonization_json = JSON.stringify(polygonization_json);

    wapi2.wapi_make_geometry (shape_id, shape_json, polygonization_json, //result_callback,
        function finished__c3(vf_result, call_id1, shape_id1) {
            vf_result.faces; vf_result.verts;
            assert(shape_id1 === shape_id);
            var ignoreDefaultNormals = true;
            var bool_reallocated = geometry.update_geometry1(vf_result.verts, vf_result.faces, ignoreDefaultNormals, false);
            geometry.use_default_normals_from_vertices();
            done_callback(shape_id1, null);
        },
        function progressed_update__c3(vf_result, call_id1, shape_id1) {
            console.error("THE UPDATE: ", vf_result, call_id1, shape_id1);
            //                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              assert(is_progress_update);
            vf_result.faces; vf_result.verts;
            assert(shape_id1 === shape_id);
            var ignoreDefaultNormals = true;
            var bool_reallocated = geometry.update_geometry1(vf_result.verts, vf_result.faces, ignoreDefaultNormals, false);
            geometry.use_default_normals_from_vertices();
            //done_callback(shape_id1, null);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            progress_update_callback(shape_id1, null);
        }
    );
};
wapi3.query_implicit_values_old = wapi2.wapi_query_implicit_values_old;

wapi3.wapi_query_a_normal = wapi2.wapi_query_a_normal;


// End of Worker-based API
// ******************************************************************************


var IMPLICIT_WORKER = wapi3;
// wapi_*

console.log("IMPLICIT_WORKER ---------------");
for(var a in IMPLICIT_WORKER) console.log("IMPLICIT_WORKER.",a);
