// it will not  deploy/
// sometimes ignore_matrix is necessary
// still some refactoring needed
// unused functions to be deleted
// unnecessary console.logs to be removed (bbox)

/*
 Worker version of implisolid.js
 Also the  Node version.
 The new callback-based structure is created, which is much more complicated, but allows separated ThreeJS-related code from the rest (core) of the code.
 The code has three levels, each one has a different service:

Level I:  IMPLICIT_LOW
 1-low level MCC2 API, which directly calls C++ functions. Can run in NodeJS and as a Worker.
Fully incode the Worker.
 Level II:
 2- ImpliSolid API functions: which call the MCC2 API functions. They also receive callbacks to receive front-end logic code (ThreeJS). Can also run in Node or as a Worker. This part however most likrely will remain on client side. This mediaes calls to level I. e.g. ideally, no direct call to Level I should be necessary.
Between the worker and front-end. NodeJS compatible.

IMPLISOLID
 Level III:
 3- The ThreeJS-related code, which runs only on browser. The methods for manipulating ThreeJS objects (LiveGeometry, setting of the normals in the relevant ThreeJS attributes, etc). This logic is defines are a few functions as well and some callbacks defined within them, which are sent to level II functions.

*/
'use strict';

// LiveBufferGeometry71  MyBufferGeometry77 LiveBufferGeometry79

var ImplicitService = (function () {
'use strict';

function init(service) {
    'use strict';

    // Low level API
    //main = Module.cwrap('main', 'number', []);
    //var service={}; //= newProducer //is an interface
    service.build_geometry = Module.cwrap('build_geometry', null, [ 'string', 'string']);
    service.get_v_size = Module.cwrap('get_v_size', 'number', []);
    service.get_f_size = Module.cwrap('get_f_size', 'number', []);
    service.get_v = Module.cwrap('get_v', null, ['number']);
    service.get_f = Module.cwrap('get_f', null, ['number']);
    service.get_v_ptr = Module.cwrap('get_v_ptr', 'number', []);
    service.get_f_ptr = Module.cwrap('get_f_ptr', 'number', []);
    service.finish_geometry = Module.cwrap('finish_geometry', null, []);

    service.set_object = Module.cwrap('set_object', 'number', ['string', 'number']);
    service.unset_object = Module.cwrap('unset_object', 'number', ['number']);
    service.set_x = Module.cwrap('set_x', 'number', ['number', 'number']);  // sends the x points for evaluation of implicit or gradient
    service.unset_x = Module.cwrap('unset_x', null, []);
    service.calculate_implicit_values = Module.cwrap('calculate_implicit_values', null, []);
    service.get_values_ptr = Module.cwrap('get_values_ptr', 'number', []);
    service.get_values_size = Module.cwrap('get_values_size', 'number', []);
    service.calculate_implicit_gradients = Module.cwrap('calculate_implicit_gradients', null, ['number']);  // boolean argument to normalize and reverse the vevtors, suitable for rendering.
    service.get_gradients_ptr = Module.cwrap('get_gradients_ptr', 'number', []);
    service.get_gradients_size = Module.cwrap('get_gradients_size', 'number', []);

    service.get_pointset_ptr = Module.cwrap('get_pointset_ptr', 'number', ['string']);
    service.get_pointset_size = Module.cwrap('get_pointset_size', 'number', ['string']);

    if (Module["_about"]) {
        service.about = Module.cwrap('about', null, []);
    } else {
        service.about = "C++ method problem: no about()";
        console.error("C++ method problem: a correct version not found");
    }
    /*
    try {
        service.about = Module.cwrap('about', null, []);
    } catch(err) {
        service.about = "C++ method problem";
        console.log("C++ method problem: not found");
    }
    */

    service.init_ = function(){
        service.needs_deallocation = false;
    }

    service.finish_with = function (){
        //after the last round.
        if(!service.needs_deallocation){
            console.error("cannot `finish_geometry()`. Geometry not produced.");
        }
        service.finish_geometry();
        service.needs_deallocation = false;
    }

    service.set_vect = function (float32Array) {
        // Accesses module
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        if (float32Array.length % 3 != 0) {console.error("bad input array");};
        var nverts = float32Array.length / 3;
        var verts_space = Module._malloc(_FLOAT_SIZE*3*nverts);
        Module.HEAPF32.subarray(verts_space/_FLOAT_SIZE, verts_space/_FLOAT_SIZE + 3*nverts).set(float32Array);
        var result = service.set_x(verts_space, nverts);
        console.log("result: "+result);
        if (!result) {
            console.error("Something went wrong: ", result);
        }
        Module._free( verts_space );
    };

    // move to service2
    /*
    service.get_emc_array = function (verts_address, point_count) {

        // todo: should work with point_count==0

        // var verts_address = implisolid_.get_pointset_ptr(name);
        // var point_count = implisolid_.get_pointset_size(name);

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var float_array = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 1 * point_count);
        if (point_count == 0) {
            console.warn("Error: zero points", point_count);
            //return // null; // make_empty_lines_mesh();
        }

        return float_array;
    }
    */

    service.init_();

    return service;
}



/**
 * Service 2:
 * Technical:
 * Wraps all accesses to Module. Suitable for weparating Worker from the main one.
 */
function init2(impli2, impli1) {

    // replaced by query_normals()?
    function prepare_gradients1(implicit_service, shape_json, verts) {
        var ignore_root_matrix = true;

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var result1 =
            implicit_service.set_object(shape_json, ignore_root_matrix);
        console.log("SET_OBJECT() RESULTS:", result1);
        implicit_service.set_vect(verts);  // overhead
        implicit_service.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
        var ptr = implicit_service.get_gradients_ptr();
        var ptr_len = implicit_service.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);

        for( var i = 0 ; i < ptr_len; i++) {
            // gradients[i]  += (Math.random() *2 -1.0 )* 0.2;
        }
       // geom.update_normals_from_array(gradients);
       return gradients;
    }
    function release_gradients1(implicit_service) {
        implicit_service.unset_x();
        var result1 =
            implicit_service.unset_object(1);
        console.log("UNSET_OBJ() RESULTS:", result1);
    }


    /**
     * Callback receives the outcome of the normals. Instead of returning the normals, they are used in the callback and then the resources are freed.
     */
    impli2.query_normals = function(x, callback) {
        console.error("what?", mp5_str, ignore_root_matrix);
        impli1.set_object(mp5_str, ignore_root_matrix);
        impli2.set_vect(x);  // overhead
        impli1.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
        var ptr = impli1.get_gradients_ptr();
        var ptr_len = impli1.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        //console.log("grad len = " +  ptr_len+ "  grad = " + gradients);  // x 4

        //for( var i = 0 ; i < ptr_len; i++) {
        //    gradients[i] += Math.random() * 0.2;
        //}


        //geom.update_normals_from_array(gradients);
        callback(gradients);

        impli1.unset_x();
        impli1.unset_object();
    };


    // mid-level API
    //
    
    /** Note that nomals are NOT queried & applied in this function.
    The inputs should be jsonified already. Although for bavkward compatibilty it is automatically converted.
    Note that in the similar "update" method, the polygonization_params must be non-JSONified. 
    the update methos is: update_geometry_from_json() or update_geometry()
    // todo: remove allocate_buffer, because it will be set in the callback
    */
    impli2.make_geometry = function (mp5_str, polygonization_params_str, geometry_callback, allocate_buffer) {
        assert(typeof geometry_callback !== 'undefined');
        assert(geometry_callback);

        if (typeof polygonization_params_str !== "string") {
            polygonization_params_str = JSON.stringify(polygonization_params_str);
            console.error("Use a string, a JSONified settings:", polygonization_params_str);
        }

        var startTime = new Date();
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT

        if(impli1.needs_deallocation) {
            impli1.finish_geometry();
            impli1.needs_deallocation = false;
        }
        
        // todo: surround build_geometry() by try{}catch{}

        //console.log("polygonization_params_str.resolution " + polygonization_params_str.resolution);
        //polygonization_params_str.resolution = 40;

        //var mp5_str = JSON.stringify(shape_params);
        //var mp5_str = JSON.stringify(shape_params);
        assert(typeof polygonization_params_str === "string");
        console.log('mp5_str', mp5_str);
        console.log('polygonization_params_str', polygonization_params_str);
        impli1.build_geometry(mp5_str, polygonization_params_str);
        impli1.needs_deallocation = true;

        var nverts = impli1.get_v_size();
        var nfaces = impli1.get_f_size();

        var verts_address = impli1.get_v_ptr();
        var faces_address = impli1.get_f_ptr();

        var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);

        // first iteration: the geometry_callback
        var geom = geometry_callback(verts, faces, allocate_buffer);


        var endTime = new Date();
        var timeDiff = endTime - startTime;

        //report_time(timeDiff, function(){hist();});

        return geom;
    };

    /**
    Function: query_implicit_values(shape_str, points, reduce_callback)
    Evaluates the points **points** in the implicit function. The function is specified as a mp5-fragment (json) string **shape_str**.
    The values are not directly returned. Instead, a callback is used to prepare a value.
    The reason a callback is used is that the resources need to be freed.

    Returned value:
    The value returned us the value returned by the callback reduce_callback.
    In examples below, the return value is:
    - Example1: returns Float32Array of implicit values in the following usage.
    - Example2: returns   true for collision, false for no collision,

    Example 1:
        var my_shape = mainModel.root.sons[#nr]
        var verts_array = new Float32Array([x1,y1,z1, x2,y2,z2, ..., xn,yn,zn])
        var implicit_vals = query_implicit_values(my_shape, points,  // false
            reduce_callback = function (values_tarray) {
                // here we create a view on the HEAP from position
                // ( ptr / _FLOAT_SIZE ) to position (ptr / _FLOAT_SIZE + ptr_len)
                var result = new Float32Array(values_tarray);  // clones (makes a copy of) data
                return result;
            });
    Example 2:
        var my_shape = mainModel.root.sons[#nr]
        var verts_array = new Float32Array([x1,y1,z1, x2,y2,z2, ..., xn,yn,zn])
        var has_positive_vals = query_implicit_values(my_shape, verts_array,  // true
            function (values_tarray) {
                var found_positive_val = false;
                for (var i=0; i< values_tarray.length; i++)
                {
                    if (values_tarray[i] >= - CONFIG.collision_detection.epsilon)
                    {
                        found_positive_val = true;
                        break;
                    }
                }
                var result =  found_positive_val;
                return result;
            });
    */
    impli2.query_implicit_values = function(mp5_str, points, reduce_callback)
    {
        assert(points instanceof Float32Array);                   /* simple check */
        assert(typeof reduce_callback === 'function');             /* the callback that produces the return value of the function */

        /* Declarations */

        var _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;                    /* We ll need the _FLOAT_SIZE in bytes, when we deal with allocations and HEAPF32*/
        var nverts = points.length / 3;                                       /* the number of vertices for which we want to calculate the implicit value */
        var verts_space_address = Module._malloc(_FLOAT_SIZE * 3 * nverts);  /* This allocates space in the C++ side, in order to create the array of vertices. */

        /* body */
        Module.HEAPF32.subarray(
            verts_space_address / _FLOAT_SIZE, verts_space_address / _FLOAT_SIZE + 3 * nverts
        ).set(points);

        // Never ignore the root matrix. (Always takes the root matrix into account).
        // Hence, query_implicit_values() does not need an `ignore_root_matrix` argument.
        const ignore_root_matrix = false;

        var obj_id = impli1.set_object(mp5_str, ignore_root_matrix);                    /* Allocations on the C++ side */
        var success = impli1.set_x(verts_space_address, nverts);
        //todo: rename "success"
        if (!success){
            console.log("set_x returned false . Probably an error in memory allocation. nverts was: " + nverts);
            impli1.unset_object(obj_id);
            return false;
        }
        // impli1.set_x_with_matrix(verts_space_address, nverts, matrix);     /* can create a function like this so we dont have the matrix issue */
        impli1.calculate_implicit_values();                                /* The actual implicit value calculation*/

        var ptr = impli1.get_values_ptr();                                 /* retrieve a pointer to the position in memory of the calculated array */
        var ptr_len = impli1.get_values_size()
        var values_tarray = Module.HEAPF32.subarray(ptr / _FLOAT_SIZE , ptr / _FLOAT_SIZE + ptr_len );

        var result = reduce_callback(values_tarray);

        //Bug! Forgot to FREE!!!
        Module._free( verts_space_address );
        impli1.unset_object(obj_id);    /* Free allocated C++ memory */
        impli1.unset_x();
        return  result;
    }

    impli2.query_a_normal = function(mp5_shape_json, point, result_callback) {
        // Always takes the root matrix into account. Similar to query_implicit_values()
        const ignore_root_matrix = false;
        var objid = impli1.set_object(mp5_shape_json, ignore_root_matrix);

        var input_verts = new Float32Array([point.x, point.y, point.z]);

        var nverts = 1;                                       // the number of vertices for which we want to calculate the implicit value
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var verts_space_address = Module._malloc(_FLOAT_SIZE * 3 * nverts);  // This allocates space in the C++ side, in order to create the array of vertices.
        Module.HEAPF32.subarray(verts_space_address / _FLOAT_SIZE, verts_space_address / _FLOAT_SIZE + 3 * nverts).set(input_verts);

        var setx_resut__address = impli1.set_x(verts_space_address, nverts);  // !
        if (!setx_resut__address){
            console.log("set_x returned false: nverts is: " + nverts);
            impli1.unset_object(objid);
            return false;
        } else {
            // console.log("OK! set_c successful");
        }

        impli1.calculate_implicit_gradients(true);  // true = normalise and invert
        var ptr = impli1.get_gradients_ptr();                                 // retrieve a pointer to the position in memory of the calculated array
        var ptr_len = impli1.get_gradients_size()
        var r = Module.HEAPF32.subarray(ptr / _FLOAT_SIZE , ptr / _FLOAT_SIZE + ptr_len );

        // console.log(r[0], r[1], r[2]);
        //output_vect3.set(-r[0], -r[1], -r[2]);
        result_callback(r);

        impli1.unset_x();
        impli1.unset_object(objid);
        Module._free( verts_space_address );
    };

    /*
     * The output is written into vf_dict, which is a dictionary with keys 'verts' and 'faces'.
     * Returns true on success, i.e. if the shape is not empty. Otherwise, returns false;
     */
    impli2.get_latest_vf = function (output_vf_dict) {
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;
        const POINTS_PER_FACE = 3;

        var nverts = impli1.get_v_size();
        var nfaces = impli1.get_f_size();

        if(nfaces > 0){
            var verts_address = impli1.get_v_ptr();
            var faces_address = impli1.get_f_ptr();

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

    /*
    get_pointset_() would belong to impli2 if it was defined in this file
    impli2.get_pointset_ = function (pointset_label) {
        var n = impli1.get_pointset_size(pointset_label);
        if (n == 0) {
            return null;
        }
        var verts = impli1.get_emc_array(
            impli1.get_pointset_ptr(pointset_label),
            n * 3 );
        return verts;
    }
    */

    impli2.about = impli1.about;
    // Enable access to lower level functions. not recommended except for debugging.
    impli2.service1 = impli1;
}




/**
 * High level API. Works with ThreeJS objects (THREE.Geometry)
 * Dependency: ThreeJS (not-explicit)
 * This may need to go inside liveGeomtry or a subclass.
 * move to front-end:  implisolid_front.js
 * This is part of WeDesign, not ImpliSolid back-end.
 * There is no mention of 'Module' here.
 * No access to impli1.
 */

function init3(service3, service2) {

    service3.use_II = true;
    service3.use_II_qem = true;
    service3.use_III = false;
    service3.use_noise = false;
    service3.repeats = 0;
    service3.custom_mc_settings = null;  // can be manually set in browser's console


    service3.receive_mesh = function (vf) {

        var vf = {faces: null, verts: null};
        var nonempty = service2.get_latest_vf(vf);
        if (nonempty) {
            // vf already contains the output of get_latest_vf()
        } else{
            console.log("empty implicit. Using a default shape.");
            vf['verts'] = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 ]);
            vf['faces'] = new Uint32Array([0,1,2, 0,2,3, 0,4,5, 0,5,1, 1,5,6, 1,6,2, 2,6,3, 3,6,7, 4,5,6, 5,6,7]);
        }
        assert(vf['verts']);
        
        assert(vf['faces']);
        return vf;
    }

    // update_geometry_from_json  <-->  update_geometry()
    /**
    useless. to be removed.
    */
    service3.update_geometry = function(geometry, ignoreDefaultNormals) {
        // call after build_geometry()

        /*
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;
        const POINTS_PER_FACE = 3;

        var nverts = service1__.get_v_size();
        var nfaces = service1__.get_f_size();

        // non-empty
        if(nfaces > 0) {
            var verts_address = service1__.get_v_ptr();
            var faces_address = service1__.get_f_ptr();

            var verts = Module.HEAPF32.subarray(
                verts_address/_FLOAT_SIZE,
                verts_address/_FLOAT_SIZE + 3*nverts);

            var faces = Module.HEAPU32.subarray(
                faces_address/_INT_SIZE,
                faces_address/_INT_SIZE + nfaces * POINTS_PER_FACE);
        }
        else{
            console.log("empty implicit");
            var verts = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 ]);
            var faces = new Uint32Array([0,1,2, 0,2,3, 0,4,5, 0,5,1, 1,5,6, 1,6,2, 2,6,3, 3,6,7, 4,5,6, 5,6,7]);
        }

        assert(verts);
        assert(faces);

        return geometry.update_geometry1(verts, faces, ignoreDefaultNormals, false);
       */

        //var implicit_service1 = impli1;
        
        // make it a callback        
        var vf = service3.receive_mesh()

        var bool = geometry.update_geometry1(vf['verts'], vf['faces'], ignoreDefaultNormals, false);
        return bool;

    };


    // update_geometry_from_json  <-->  update_geometry()
    /**
        @param ignoreDefaultNormals=true always use `true`.
        
        Some confusion here:
        update_geometry_from_json() versus update_geometry()
        
        
        todo: fixme: can be simplified a lot. Also remove update_geometry().  See implisolid_worker.js
    */
    service3.update_geometry_from_json_experimental = function(geometry, shape_json, mc_params, ignoreDefaultNormals) {
        // see __update_reused_geometry()
        var startTime = new Date();
        if (!shape_json) {
            console.error(shape_json);
        }

        if (typeof shape_json !== "string") {
            shape_json = JSON.stringify(shape_json);
        }

        // mc_params must nt be jsonified, because we need its `ignore_root_matrix` field for querying the normals.
        assert(typeof mc_params !== "string");
        mc_params = JSON.stringify(mc_params);

        // fixme: this creates a new geometry? No. But compare with update_geometry()
        var allocate_buffer;  // =?
        service2.make_geometry (shape_json, mc_params,
            function (verts, faces, allocate_buffers) {
                /*if (allocate_buffers === undefined) {
                    allocate_buffers = true;
                    console.error("allocate_buffers not specified:", allocate_buffers);
                }*/
                assert(allocate_buffers == false);

                // Based on the latest vf_pair in the most recent step (i.e. in make_geometry() ) which accept a (shape,polyg) pair.
                // Does Not query the normals. It's queries in a separate request. todo: do the normals in the same request but it must remain optional. The normals are likely to be useless in future, so we should not waste the CPU for it. In that case, we can have it in the same polygonization_setttings_json string.
                service3.update_geometry (geometry, ignoreDefaultNormals);

                var ignore_root_matrix = mc_params.ignore_root_matrix;

                // Why normals are not already taken care of in update_geometry()? Is the following line redundant?
                service3.make_normals_into_geometry(geometry, shape_json, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

                var endTime = new Date();
                var timeDiff = endTime - startTime;
                //console.log(nv3 + " , " + nf3);
                //console.log("Time: "+timeDiff+ " msec.");

                report_time(timeDiff);

                return geometry;
            },
        allocate_buffer=false);

        // call after build_geometry()

        // Why? Note that the geometry may be  not prepared yet.
        // This might run after the geometry is finished (because the callback of the make_geometry() is not really async currently).
        // But this means update_geometry() is called twice. why?
        // If it was async, it was there to have an immediate output?
        service3.update_geometry(geometry, ignoreDefaultNormals);
    };

    service3.update_geometry_from_json = service3.update_geometry_from_json_experimental;

    //was: geom. update_normals (service ..)
    // todo: move the callback into LiveGeometry
    service3.make_normals_into_geometry = function(geom, mp5_str, x, ignore_root_matrix) {
        geom.use_default_normals_from_vertices();
        return;

        // console.error(x);
        /*
        for( var i = 0 ; i < x.length; i++) {
            x[i] += Math.random() * 0.9;
        }
        */

        // why query_normals does not use ignore_root_matrix?
        service2.query_normals(x, function(gradients) {
            //for( var i = 0 ; i < ptr_len; i++) {
            //    gradients[i] += Math.random() * 0.2;
            //}
            geom.update_normals_from_array(gradients);
        });

        /*
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        service1.set_object(mp5_str, ignore_root_matrix);
        service2.set_vect(x);  // overhead
        service1.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
        var ptr = service.get_gradients_ptr();
        var ptr_len = service.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        //console.log("grad len = " +  ptr_len+ "  grad = " + gradients);  // x 4

        //for( var i = 0 ; i < ptr_len; i++) {
        //    gradients[i] += Math.random() * 0.2;
        //}

        geom.update_normals_from_array(gradients);

        service.unset_x();
        service.unset_object();
        */

    };

    /**
       This method generates a polygonisation_settings based on the given (bbox, ignore_root_matrix) 
       @param `ignore_root_matrix` is embeded in the result.
       This method must move into designer code. (naturally). Also custom_mc_settings, use_II, etc
    */
    service3.make_polygonization_settings = function(bbox, ignore_root_matrix) {

        assert(bbox, "You need to specify the bounding box");
        assert("min" in bbox, "You need to specify *.min.x");
        assert("max" in bbox, "You need to specify *.max.x");
        assert("x" in bbox["min"], "You need to specify *.max.x");

        var bb ={};
        var sc = 1.0;
        var test = 0.;
        bb["xmin"] = bbox.min.x * sc + test;
        bb["xmax"] = bbox.max.x * sc - test;

        bb["ymin"] = bbox.min.y * sc + test;
        bb["ymax"] = bbox.max.y * sc - test;

        bb["zmin"] = bbox.min.z * sc + test;
        bb["zmax"] = bbox.max.z * sc - test;

        _expect(bb["xmin"], "boundingbox is null");
        _expect(bb["xmax"], "boundingbox is null");
        _expect(bb["ymin"], "boundingbox is null");
        _expect(bb["ymax"], "boundingbox is null");
        _expect(bb["zmin"], "boundingbox is null");
        _expect(bb["zmax"], "boundingbox is null");

        // Designer-specific
        var getResolution  = function(bb) {
            return CONFIG.implisolid.default_mc_resolution;
            const max_value = 40;
            const min_value = 14;
            const factor = CONFIG.implisolid.default_mc_resolution;
            var max_length = Math.max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"], bb["zmax"] - bb["zmin"]);
            var tmp =  Math.min(max_value,max_length*factor);
            return  Math.floor(Math.max(tmp,min_value));

            // 1 -> 28
            // 2 -> 48
        };


        var mc_res = CONFIG.implisolid.default_mc_resolution;
        /*
        var mc_properties = {
            resolution: getResolution(bb),
            box: bb,
            ignore_root_matrix: ignore_root_matrix,

            vresampl: {iters: this3_.use_II? 1:0, c: 1.0},
            projection: {enabled: this_.use_II? 1:0},
            qem: {enabled: (this3_.use_II && this3_.use_II_qem)?1:0},

            subdiv: {enabled: 1},

            overall_repeats: 1,

            debug: {
                post_subdiv_noise: 0.01,
            }

        };
        */

        var mc_properties = {
            resolution: getResolution(bb),
            box: bb,
            ignore_root_matrix: ignore_root_matrix,

            vresampl: service3.use_II? {iters: 1, c: 0.4} : {iters: 0, c: 1.0},
            projection: {enabled: service3.use_II? 1 : 0},
            qem: {enabled: service3.use_II_qem? 1 : 0},
            //subdiv: {enabled: 1},
            //overall_repeats: 2,
            subdiv: {enabled: service3.use_III? 1 : 0},
            overall_repeats: service3.repeats,

            debug: {
                enabled_pointsets: 0,
                post_subdiv_noise: service3.use_noise? 0.01 : 0.0,
            },
        };

        if (service3.custom_mc_settings) {
            // nonrecursive is enough, because we want to replace everything except for "box".
            // note: this is not tested.
            // Example usage:
            //    IMPLICIT.custom_mc_settings = {vresampl: {iters: 1, c: 1} };
            mc_properties = merge_dicts_nonrecursive(mc_properties, service3.custom_mc_settings);
        }
        return mc_properties;
    }

    /** Similar to getLiveGeometry_from_json(), but the full polygonization_setttings_json_str is proovided. 
    */
    service3.getLiveGeometry_from_json = function(mp5_str, polygonization_setttings_json_str, ignore_root_matrix) {
        // contradiction: needs ignore_root_matrix but already has it in polygonization_setttings_json_str??
        assert(typeof mp5_str === "string");
        assert(typeof polygonization_setttings_json_str === "string");

        return service2.make_geometry(mp5_str, polygonization_setttings_json_str,

            function (verts, faces, allocate_buffers) {
                // ThreeJS-specific code

                // var ignore_root_matrix: Does not need other (MC-related) arguments.

                if (allocate_buffers === undefined) {
                    allocate_buffers = true;
                    console.error("allocate_buffers not specified:", allocate_buffers);
                }

                var geom = new LiveBufferGeometry79(verts, faces, allocate_buffers);

                // Set the normals
                // var ignore_root_matrix = mc_params.ignore_root_matrix;  // Does not need other (MC-related) arguments.
                //geom.update_normals(this_, verts, mp5_str, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.
                //
                if (!mp5_str)
                    console.error(mp5_str);

                if (true) {
                    service3.update_geometry(geom, true);
                }

                // should not be used I think.
                service3.make_normals_into_geometry(geom, mp5_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

                //this_.aaaaaaaaaA(verts);

                return geom;
            }
        );
    }

    // Example usage of implisolid
    /**
    This method is called by the designer to obtain the geometry from the ImplicitService
    Usage: service.getLiveGeometry = function(dict, bbox, ignore_root_matrix)
    The getLiveGeometry() := a getLiveGeometry_from_json() + a make_polygonization_settings()
    */
    service3.getLiveGeometry = function(shape_properties, bbox, ignore_root_matrix) {

        var mc_properties = service3.make_polygonization_settings(bbox, ignore_root_matrix);

        // console.log (" mc properties : " + JSON.stringify(mc_properties));
        var mp5_str;
        if (typeof shape_properties === "string")
            mp5_str = shape_properties;
        else
            mp5_str = JSON.stringify(shape_properties);

        return service3.getLiveGeometry_from_json(mp5_str,  JSON.stringify(mc_properties));
    };

    if (typeof PS_UTILS !== "undefined")
        service3.PS_UTILS = PS_UTILS;
    else
        service3.PS_UTILS = {};


    service3.query_implicit_values = service2.query_implicit_values;

    service3.query_a_normal = service2.query_a_normal;

    service3.about = service2.about;

    // Enable access to lower level functions. not recommended except for debugging.
    // Makes sense because service3 is useful (we might want to use it with other callbacks)
    service3.service2 = service2;
}

var ImplicitService = function() {

    var impli1 = {};
    // adds the low-level API to 'this'
    init(impli1);

    var impli2 = {};
    // adds the mid-level API to 'this'
    init2(impli2, impli1);

    var impli3 = this;
    // adds the high-level API to 'this'
    init3(impli3, impli2);

    // API interface functions:
    console.log("impli1 ---------------");
    for(var a in impli1) console.log("API Level1.",a);
    console.log("impli2 ---------------");
    for(var a in impli2) console.log("API Level2.", a);
    console.log("impli3 ---------------");
    for(var a in impli3) console.log("IMPLICIT.",a);
};


// return IMPLICIT;  // is null!
//return _on_cpp_loaded;
// return null;
return ImplicitService;
}());

/* most times you need to write the following function */
var IMPLICIT = null;  // is assigned to at _on_cpp_loaded();
function _on_cpp_loaded() {
    console.log("C++ ready.");
    IMPLICIT = new ImplicitService();
    // IMPLICIT = new ImplicitWorkerService();

    // combines the IMPLICIT as a Worker/Node/npm library with the ThreeJS part. Not good! Solution: divide into two classes.
    // ImpliSolid.js, ...ImpliSolid3js.js (frontend)
    // Alternative: Globally: Module._on_cpp_loaded = ...;

    // ugly
    if (typeof _assert_000 !== "undefined")
        assert = _assert_000;

    return IMPLICIT;  // another usage
};

/**
 * Usage:
 * console.log(merge({}, {key1: 1,    key2: "test",    key3: [5, 2, 76, 21]}, { key1: 1 }));
 * http://stackoverflow.com/a/8625261/4374258
 */
var merge_dicts_nonrecursive = function() {
    var merged_dict = {};
    var nargs = arguments.length;
    for (var arg_i = 0; arg_i < nargs; arg_i++) {
        var newarg = arguments[arg_i];
        for (var key in newarg) {
            if (newarg.hasOwnProperty(key)) {
                merged_dict[key] = newarg[key];
            }
        }
    }
    return merged_dict;
};

/**
 * Confusion:
 * IMPLICIT is an object
 * of type ImplicitService
 * ImplicitService is global
 * IMPLICIT is global too.
 * Do we really need an instancing a class for IMPLICIT ? i.e. is there really a need for separation of class and object?
 * Instantiation is done in _on_cpp_loaded(). It cannot have any parameter. It is called automatically as a callback by Emscripten. In which, a global variable IMPLICIT is set.
 * No code should/can be executed prior to _on_cpp_loaded. todo: make ((function () {..}());) a constructor.
 * Rename? IMPLISOLID.
 */



/*

impli1: Low level API:
=====================
build_geometry
get_v_size
get_f_size
get_v
get_f
get_v_ptr
get_f_ptr
finish_geometry

set_object
unset_object
set_x
unset_x
calculate_implicit_values
get_values_ptr
get_values_size
calculate_implicit_gradients
get_gradients_ptr
get_gradients_size
get_pointset_ptr
get_pointset_size
about
init_
finish_with
set_vect
needs_deallocation


impli2: Mid level API:  (Can run in NodeJS)
=====================
query_implicit_values  (string, array, callback_function)
query_normals (array, function)
query_a_normal (json_str, Float32Array, result_callback_function)

make_geometry (string, string, function, bool)
get_latest_vf ({array,array})


impli3: High level API: (Cannot run on NodeJS)
=====================
IMPLICIT.

getLiveGeometry_from_json(string_json, string_json, bool) -> Geometry
update_geometry(geometry, ignoreDefaultNormals) -> void
update_geometry_from_json(THREE.Geometry, json_string, json_string, bool) -> void

make_normals_into_geometry(THREE.Geometry, json_string, Float32Array, bool)

query_implicit_values = inherit.
query_a_normal = inherit
about = inherit

//todo: raycast.

Remove from IMPLICIT: (move to designer)
    getLiveGeometry(dict, dict, bool) -> Geometry
         (only keep getLiveGeometry_from_json() )

    make_polygonization_settings()

    The configs properties: (which override the polygonized_settings)
        .use_II
        .use_II_qem
        .use_III
        .use_noise
        .repeats
        .custom_mc_settings  // this should replace the other ones

    merge_dicts_nonrecursive()


IMPLICIT.PS_UTILS
IMPLICIT.service2


--------------------
API todo:
rename   getLiveGeometry() --> newLiveGeometry
rename   getLiveGeometry_from_json() --> newLiveGeometry_from_json
replace {    .use_II, .use_II_qem, .use_III, .use_noise, .repeats} by custom_mc_settings.
*/
