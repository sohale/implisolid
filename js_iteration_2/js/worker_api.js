'use strict';
// importScripts('../build/mcc2.compiled.js');
// GET ... .js.mem 404 (File not found)
// mcc2.compiled.js:96101 Uncaught could not load memory initializer mcc2.compiled.js.mem
// Module.filePackagePrefixURL  (prefix only)
// Module.locateFile  (fiull filename path; callback)
var Module={};
//Module.filePackagePrefixURL = ... // this is ignored
Module.locateFile = function(){return '../../build/mcc2.compiled.js.mem';};
importScripts('../../build/mcc2.compiled.js');

'use strict';

onmessage = function(event) {
    var api = wwapi;
    // console.log('Message received from main script. event.data=',event.data);
    var data = event.data;
    var _call_id = data.call_id;  // call_identification
    var _callback_id = event.data.callbackId;

    var _progressCallback_id = event.data.progressCallbackId;
    var progressive_completion = !!_progressCallback_id;  // whether it is undefined

    assert(!progressive_completion || typeof _progressCallback_id == 'number', "specified progress callback id sould be an integer");
    switch (data.funcName) {
        case "query_implicit_values":
            var shape_id = data.obj_id;  // data.obj_req_id
            var result =
                api.query_implicit_values(shape_id, data.mp5_str, data.points, data.reduce_callback);
            postMessage({return_callback_id:_callback_id, returned_data: result, call_id: _call_id, shape_id:shape_id});  // {result_allpositive:}
            break;

        case "make_geometry":
        //case "make_geometry_u":
            // the only progressive
            var shape_id = data.obj_req_id;
            var result =
                api.make_geometry__workerside(shape_id, data.mp5_str, data.polygonization_settings, _call_id, {progressCallback_id: _progressCallback_id, call_id: _call_id, shape_id: shape_id});  // _call_id not needed here actually
            postMessage({return_callback_id:_callback_id, returned_data: result, call_id: _call_id, shape_id:shape_id});  // {result_allpositive:}
            break;

        case "get_latest_vf":
            var shape_id = data.shape_id;
            var result_vf_output = {faces: null, verts: null};

            var nonempty = api.get_latest_vf(shape_id, result_vf_output, _call_id);  // _call_id not used here actually
            result_vf_output.nonempty = nonempty;
            // if not nonempty, returns nonempty==false
            postMessage({return_callback_id:_callback_id, returned_data: result_vf_output, call_id: _call_id, shape_id: shape_id});  // {result_allpositive:}
            // definitely box the dict
            break;

        case "query_implicit_values_old":
            var shape_id = data.shape_id;
            assert(
                //data.reduce_type === "all-outside" ||
                //data.reduce_type === "all-non-positive" ||
                //data.reduce_type === "<=0"
                data.reduce_type === "any-inside" ||
                data.reduce_type === "any-positive" ||
                data.reduce_type === ">0"
            );

            var epsilon = data.epsilon;

            // old-style. The new style should get the string on the C++ side and do the reduction there.

            var reduce_callback = function(values_tarray) {
                'use asm';
                const minus_eps = +(-epsilon);
                var found_positive_val = false|0;
                for (var i = 0|0; i < values_tarray.length|0; i++) {
                    if (values_tarray[i] >= +minus_eps ) {
                        found_positive_val = true|0;
                        break;
                    }
                }
                return found_positive_val;
            };

            var result_values = api.query_implicit_values_old(shape_id, data.mp5_str, data.points, reduce_callback, _call_id);
            // result_values.shape_id = shape_id;  // where this should be put?

            postMessage({return_callback_id:_callback_id, returned_data: result_values, call_id: _call_id, shape_id: shape_id});  // {result_allpositive:}
            break;

        case "query_normals":
            var shape_id = data.shape_id;
            var result_normal = new Float32Array(3);
            result_normal[0] = undefined;  // FOR DEBUG
            var is_definitely_called = false;
            var result_callback = function(gradient_vector) {
                result_normal.set(gradient_vector);
                is_definitely_called = true;
            }
            api.query_normals(shape_id, data.mp5_str, data.points, result_callback, _call_id);
            assert(result_normal[0] === result_normal[0]);  // check it's not NaN
            assert(is_definitely_called === true);

            postMessage({return_callback_id:_callback_id, returned_data: result_normal, call_id: _call_id, shape_id: shape_id});  // {result_allpositive:}
            break;

        default:
            console.error("Unknown worker request for unknown function call: \""+(data.funcName)+"\" via event data:", event.data);
            throw new Error("Unrecognised worker request (unknown function call)");
    }
    //postMessage({result_allpositive: result});
    // postMessage({result_data: event.data});
/*
           {
                funcName: 'query_implicit_values',
                callbackId: callbackId,
                data: {mp5_str: '', points: new Floar32Array(), reduce_callback: 'allpositive'},
            }
*/
}

var assert__ = function(x, m){if(!x) {console.error(m,x);throw m;}}

// contains those functions in "Module" that have string arguments.
var Module_cwrapped = {
    //Based on API Version 1:
    // only functions that receive 'string' arguments need to be cwrap()ed.
    set_object: Module.cwrap('set_object','number',['string','number']),
    build_geometry: Module.cwrap('build_geometry', null, [ 'string', 'string']),
    build_geometry_u: Module.cwrap('build_geometry_u', null, [ 'string', 'string', 'string']),

    // API Version 2:
};

// Level 2
var wwapi = {};
wwapi._module = Module;

wwapi.needs_deallocation = false;  // state

    wwapi.query_implicit_values = function(shape_id, mp5_str, points, reduce_type, call_id_arg)
    {
        assert__(points instanceof Float32Array);
        assert__(points.length % 3 === 0);
        assert__(typeof reduce_type === 'string');
        // assert__(['allpositive',].indexOf(redureduce_typereduce_typece_type) > -1);             /* the callback that produces the return value of the function */
        const ENUM_ALL_NON_POSITIVE = 2;
        const ENUM_ALL_NEGATIVE = 3;
        const ENUM_ERROR = 123456;
        var reduce_type_enum =
        (    (reduce_type === "any-inside" ||
                reduce_type === "any-positive" ||
                reduce_type === ">0" ) ? ENUM_ALL_NON_POSITIVE
            : (
                (reduce_type === "any-negative" || reduce_type === ">=0") ? ENUM_ALL_NEGATIVE
                : (
                    ENUM_ERROR
                )
            )
        );

        assert(reduce_type_enum === ENUM_ALL_NON_POSITIVE);

        var _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;                    /* We ll need the _FLOAT_SIZE in bytes, when we deal with allocations and HEAPF32*/
        var nverts = points.length / 3;                                       /* the number of vertices for which we want to calculate the implicit value */
        /* This allocates enough space for vertices (points) in the C++ memory. */
        var verts_space_address = Module._malloc(_FLOAT_SIZE * 3 * nverts);
        /* save the given values into C++ memory. The index in that memory will be used for calling the actual function. */
        /* Note:
            This seems to be too many steps: malloc, copy the Float32Array (which is already made by the caller), then_set_x() copies it into a temporary memory location. The _set_x() could have been avoided by sending a pointer into _calculate_implicit_values(). But it would need to also allocate the resulting values. If the Reduce (see below) is done on C++, this second solution will be more relevant, because we won't need to allocate space for the results. Because they don't need to be communicated back. It would just need to send back the reduces value (a single value) which is a boolean.
            We must copy (using .set() ), because, although both are of type Float32Array, but points cannot be alraedy on C++ memory. The reason is, in that case it would need to have called Module._malloc() before the current function.
        */
        Module.HEAPF32.subarray(verts_space_address / _FLOAT_SIZE, verts_space_address / _FLOAT_SIZE + 3 * nverts
            ).set(points);

        var ignore_root_matrix = false;
        /* obj_id is an integer that will match to the function _unset_object(obj_id). This integer will be used to identify one object among multiple objects that are available on the C++ side. */
        var obj_id = /* Module. */ Module_cwrapped.set_object(mp5_str, ignore_root_matrix);

        /* Consume the given vertices. The vertices are given as an integer, the index of their start in C++ mamory. */
        var success = Module._set_x(verts_space_address, nverts);
        /* If there are a problem, such as memory roblems, the result is 0. */
        //todo: rename "success"
        if (!success) {
            console.log("set_x returned false . Probably an error in memory allocation. nverts was: " + nverts);
            Module._unset_object(obj_id);
            return false;
        }
        Module._calculate_implicit_values();                                /* The actual implicit value calculation*/

        var ptr = Module._get_values_ptr();           /* retrieve a pointer to the position in C++ memory, containing the results of the calculated array of implicit_values. */
        var ptr_len = Module._get_values_size();   /* We need the size too. This is a separate function call because we cannot return more than one value, or a struct in Emscripten (There may be a better way such as using a struct). */
        var values_tarray = Module.HEAPF32.subarray(ptr / _FLOAT_SIZE , ptr / _FLOAT_SIZE + ptr_len ); /* get a 'view' of that part of C++ memory as a Float32Array object. */

        /* Generate one boolean from the array. It resembles reduce() in MapReduce. */
        /* The following code can be moved into the C++ code. In this case, we don't need to query get_values_ptr(), etc. It will also be faster on C++ side. */
        var result = null;
        switch (reduce_type_enum) {
            case ENUM_ALL_NON_POSITIVE:

    	    var reduce_callback = function(values_tarray) {
                const CONFIG__collision_detection__epsilon = 0.001;
        		var found_positive_val = false;
        		for (var i=0; i < values_tarray.length; i++)
        		{
        		    if (values_tarray[i] >= - CONFIG__collision_detection__epsilon)
        		    {
        		        found_positive_val = true;
        		        break;
        		    }
        		}
        		var result =  found_positive_val;
        		return result;
    	    };

            result = reduce_callback(values_tarray);

            break;
           otherwise:
               console.error("unrecognised reduced type: ", reduce_type);
        }

        Module._free( verts_space_address );  /* Free space allocated by JS (here) on C++ memory */
        Module._unset_object(obj_id);    /* Ask C++ to free its resources; i.e. the implicit_function object. */
        Module._unset_x(); /* Ask C++ to free the memory it used to pass the verts. */
        return  result;
    };

  // in Style of an older version, although it was written later.    
    wwapi.query_implicit_values_old = function(shape_id, mp5_str, points, reduce_callback, call_id)
    {
        //todo: inside of this function needs refactoring.
        assert(points instanceof Float32Array);
        assert(typeof reduce_callback === 'function');

        assert(typeof mp5_str === "string");

        /* Declarations */

        var _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var nverts = points.length / 3;
        var verts_space_address = Module._malloc(_FLOAT_SIZE * 3 * nverts);

        /* body */
        Module.HEAPF32.subarray(
            verts_space_address / _FLOAT_SIZE, verts_space_address / _FLOAT_SIZE + 3 * nverts
        ).set(points);

        const ignore_root_matrix = false;

        var obj_id = Module_cwrapped.set_object(mp5_str, ignore_root_matrix);
        var success = Module. _set_x(verts_space_address, nverts);
        //todo: rename "success"
        if (!success){
            console.error("error allocating x: nverts is: " + nverts);
            Module. _unset_object(objid);
            Module._free( verts_space_address );
            return false;
        }

        Module. _calculate_implicit_values();

        var ptr = Module. _get_values_ptr();
        var ptr_len = Module. _get_values_size()
        var values_tarray = Module.HEAPF32.subarray(ptr / _FLOAT_SIZE , ptr / _FLOAT_SIZE + ptr_len );

        var result = reduce_callback(values_tarray);

        Module._free( verts_space_address );
        Module. _unset_object(obj_id);
        Module. _unset_x();
        return  result;
    }

    /** called by query_a_normal() */
    wwapi.query_normals = function(shape_id, mp5_shape_json, input_verts, result_callback, call_id) {
        // Always takes the root matrix into account. Similar to query_implicit_values()
        assert(input_verts instanceof Float32Array);
        var nverts = input_verts.length / 3;             

        const ignore_root_matrix = false;
        var objid = Module_cwrapped.set_object(mp5_shape_json, ignore_root_matrix);

        // var input_verts = new Float32Array([point.x, point.y, point.z]);

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var verts_space_address = Module._malloc(_FLOAT_SIZE * 3 * nverts);  // This allocates space in the C++ side, in order to create the array of vertices.
        Module.HEAPF32.subarray(verts_space_address / _FLOAT_SIZE, verts_space_address / _FLOAT_SIZE + 3 * nverts).set(input_verts);

        var allocation_success = Module. _set_x(verts_space_address, nverts);
        if (!allocation_success) {
            console.error("error allocating x: nverts is: " + nverts);
            Module. _unset_object(objid);
            Module._free( verts_space_address );
            return false;
        } else {
            // console.log("OK! set_c successful");
        }

        Module. _calculate_implicit_gradients(true);  // true = normalise and invert
        var ptr = Module. _get_gradients_ptr();                                 // retrieve a pointer to the position in memory of the calculated array
        var ptr_len = Module. _get_gradients_size()
        var r = Module.HEAPF32.subarray(ptr / _FLOAT_SIZE , ptr / _FLOAT_SIZE + ptr_len );

        // console.log(r[0], r[1], r[2]);
        //output_vect3.set(-r[0], -r[1], -r[2]);
        result_callback(r);

        // dos this also free the gradients_ptr ?
        Module. _unset_x();
        Module. _unset_object(objid);
        Module._free( verts_space_address );
    };




//based on  implisolid.js
    wwapi.make_geometry__workerside = function (shape_id, mp5_str, polygonization_settings_json, /*result_callback,*/ call_id_arg, progress_update_specs) {
        // assert(typeof result_callback !== 'undefined');
        // assert(result_callback);
        assert(typeof call_id_arg !== 'undefined');
        assert(polygonization_settings_json);
        assert(typeof polygonization_settings_json === "string");

        //progress_update_specs: {progressCallback_id: .., call_id:.., shape_id: ...}
        var progressive_completion = !!progress_update_specs.progressCallback_id;  // whether it is undefined
        // the only place call_id is usedful here. But we have it through progress_update_specs, so probably we don't need it as a call_id_arg arg, unless for debugging purposes.
        //if (!progressive_completion)  progress_update_specs = {};
        //var call_specs = JSON.stringify(progress_update_specs);  // string
        var call_specs = undefined;
        if (progressive_completion) call_specs = JSON.stringify(progress_update_specs);

        var startTime = new Date();
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT
        if(wwapi.needs_deallocation) {
            /*Module_cwrapped.*/ Module. _finish_geometry();
            wwapi. needs_deallocation = false;
        }
        //console.log("mc_params.resolution " + mc_params.resolution);
        //mc_params.resolution = 40;
        //var mp5_str = JSON.stringify(shape_params);
        //var mp5_str = JSON.stringify(shape_params);
        if (progressive_completion) {
            // should call_specs be mixed with polygonization_settings_json or kept separate?
            /*Module.*/ Module_cwrapped. build_geometry_u(mp5_str, polygonization_settings_json, call_specs); 
        } else {
            // should call_specs be mixed with polygonization_settings_json or kept separate?
            assert(call_specs === undefined);
            /*Module.*/ Module_cwrapped. build_geometry(mp5_str, polygonization_settings_json, call_specs);  // undefined!
        }
        wwapi. needs_deallocation = true;
        var nverts = Module. _get_v_size();
        var nfaces = Module. _get_f_size();
        var verts_address = Module. _get_v_ptr();
        var faces_address = Module. _get_f_ptr();
        var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);
        // first iteration: using the callback
        //var geom = result_callback(verts, faces);
        var endTime = new Date();
        var timeDiff = endTime - startTime;
        //report_time(timeDiff, function(){hist();});
        //return geom;

        return {verts:verts, faces:faces};
    };

    wwapi.get_latest_vf = function (shape_id, output_vf_dict, call_id) {
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;
        const POINTS_PER_FACE = 3;

        // todo: use shape_id
        var nverts = Module. _get_v_size();
        var nfaces = Module. _get_f_size();

        if(nfaces > 0){
            var verts_address = Module. _get_v_ptr();
            var faces_address = Module. _get_f_ptr();

            var verts = Module.HEAPF32.subarray(
                verts_address/_FLOAT_SIZE,
                verts_address/_FLOAT_SIZE + 3*nverts);

            var faces = Module.HEAPU32.subarray(
                faces_address/_INT_SIZE,
                faces_address/_INT_SIZE + nfaces * POINTS_PER_FACE);
            // return {faces: faces, verts: verts};
            output_vf_dict['faces'] = faces;
            output_vf_dict['verts'] = verts;
            return true;
        } else {
            // empty
            output_vf_dict['faces'] = null;
            output_vf_dict['verts'] = null;
            return false;
        }
    }

/** called by the Emscripten */
wwapi.send_progress_update = function (verts, faces, progres_update_callback_id, shape_id, call_id) {
    var endTime = new Date();
    //var timeDiff = endTime - startTime;
    //var result = {verts:verts, faces:faces};
    
    console.error("COPY CONSTR");
    var verts_a = new Float32Array(verts);  // suggested by Brion Vibber
    var faces_a = new Uint32Array(faces);  // suggested by Brion Vibber
    var result = {verts:verts_a, faces:faces_a};

    wwapi = 0;

    postMessage(
        {return_callback_id:progres_update_callback_id, returned_data: result, call_id: call_id, shape_id:shape_id, is_progress_update: true}
        //, [verts_a, faces_a]   // transfer list
        //, [result]
    );
}

// worker side
console.info("worker js loaded");
