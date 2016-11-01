/*
 Worker version of implisolid.js
 Also the  Node version.
 The new callback-based structure is created, which is much more complicated, but allows separated ThreeJS-related code from the rest (core) of the code.
 The code has three levels, each one has a different service:
 Level I:
 1-low level MCC2 API, which directly calls C++ functions. Can run in NodeJS and as a Worker.
 Level II:
 2- ImpliSolid API functions: which call the MCC2 API functions. They also receive callbacks to receive front-end logic code (ThreeJS). Can also run in Node or as a Worker. This part however most likrely will remain on client side. This mediaes calls to level I. e.g. ideally, no direct call to Level I should be necessary.
 Level III:
 3- The ThreeJS-related code, which runs only on browser. The methods for manipulating ThreeJS objects (LiveGeometry, setting of the normals in the relevant ThreeJS attributes, etc). This logic is defines are a few functions as well and some callbacks defined within them, which are sent to level II functions.
*/
'use strict';

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

    service.init = function(){
        service.needs_deallocation = false;

        this.use_II = true;
        this.use_II_qem = true;
    }

    service.finish_with = function (){
        //after the last round.
        if(!this.needs_deallocation){
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
        var result = this.set_x(verts_space, nverts);
        console.log("result: "+result);
        if (!result) {
            console.error("Something went wrong: ", result);
        }
        Module._free( verts_space );

    };

    service.update_geometry = function(geometry, ignoreNormals) {

        var implicit_service = this;
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;
        const POINTS_PER_FACE = 3;

        var nverts = implicit_service.get_v_size();
        var nfaces = implicit_service.get_f_size();

        if(nfaces > 0){
            var verts_address = implicit_service.get_v_ptr();
            var faces_address = implicit_service.get_f_ptr();

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

        return geometry.update_geometry1(verts, faces, ignoreNormals, false);
    };

    /**
     * Callback receives the outcome of the normals. Instead of returning the normals, they are used in the callback and then the resources are freed.
     */
    service.query_normals = function(callback) {
        service.set_object(mp5_str, ignore_root_matrix);
        service.set_vect(x);  // overhead
        service.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
        var ptr = service.get_gradients_ptr();
        var ptr_len = service.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        //console.log("grad len = " +  ptr_len+ "  grad = " + gradients);  // x 4

        //for( var i = 0 ; i < ptr_len; i++) {
        //    gradients[i] += Math.random() * 0.2;
        //}


        //geom.update_normals_from_array(gradients);
        callback(gradients);

        service.unset_x();
        service.unset_object();
    };

    //was: geom. update_normals (service ..)
    service.make_normals_into_geometry = function(geom, mp5_str, x, ignore_root_matrix) {
        geom.removeAttribute('normal');
        geom.computeVertexNormals();
        geom.__set_needsUpdate_flag(false);
        return;


        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        // console.error(x);
        /*
        for( var i = 0 ; i < x.length; i++) {
            x[i] += Math.random() * 0.9;
        }
        */

        service.query_normals(function(gradients) {
            //for( var i = 0 ; i < ptr_len; i++) {
            //    gradients[i] += Math.random() * 0.2;
            //}
            geom.update_normals_from_array(gradients);
        });

        /*
        service.set_object(mp5_str, ignore_root_matrix);
        service.set_vect(x);  // overhead
        service.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
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


    service.init();

    return service;
}

var ImplicitService = function() {

    // adds the low-level API to 'this'
    init(this);

    // High-level API
    this.make_geometry = function (mp5_str, mc_params, callback) {
        var this_ = this;
        if (typeof callback === 'undefined') {
            callback = function (verts, faces) {
                var allocate_buffers = true;
                var geom = new LiveBufferGeometry79(verts, faces, allocate_buffers);

                // Set the normals
                var ignore_root_matrix = mc_params.ignore_root_matrix;  // Does not need other (MC-related) arguments.
                //geom.update_normals(this, verts, mp5_str, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

                this_.make_normals_into_geometry(geom, mp5_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

                //this.aaaaaaaaa(verts);

                return geom;
            }
        }
        var startTime = new Date();
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT

        if(this.needs_deallocation) {
            this.finish_geometry();
            this.needs_deallocation = false;
        }

        //console.log("mc_params.resolution " + mc_params.resolution);
        //mc_params.resolution = 40;

        //var mp5_str = JSON.stringify(shape_params);
        //var mp5_str = JSON.stringify(shape_params);
        this.build_geometry(mp5_str, JSON.stringify(mc_params));
        this.needs_deallocation = true;

        var nverts = this.get_v_size();
        var nfaces = this.get_f_size();

        var verts_address = this.get_v_ptr();
        var faces_address = this.get_f_ptr();

        var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);

        // first iteration: the callback
        var geom = callback(verts, faces);


        var endTime = new Date();
        var timeDiff = endTime - startTime;

        //report_time(timeDiff, function(){hist();});

        return geom;
    };


};


// return IMPLICIT;  // is null!
//return _on_cpp_loaded;
// return null;
return ImplicitService;
}());


var IMPLICIT = null;  // is assigned to at _on_cpp_loaded();
function _on_cpp_loaded() {
    console.log("C++ ready.");
    //IMPLISOLID.
    // global
    IMPLICIT = new ImplicitService();
    // combines the IMPLICIT as a Worker/Node/npm library with the ThreeJS part. Not good! Solution: divide into two classes.
    // ImpliSolid.js, ...ImpliSolid3js.js (frontend)
    // Alternative: Globally: Module._on_cpp_loaded = ...;
    IMPLICIT.getLiveGeometry = IMPLICIT_getLiveGeometry;

    assert = _assert_000;
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
 */

/*
Todo: move to front-end:  implisolid_front.js
IMPLICIT.getLiveGeometry = function(dict, bbox, ignore_root_matrix)
*/



// Example usage of implisolid
// This method is called by the designer to obtain the geometry from the ImplicitService

var IMPLICIT_getLiveGeometry = function(dict, bbox, ignore_root_matrix) {
    var shape_properties = dict;
    var s = 1;

    assert(bbox, "You need to specify the bounding box");
    assert("min" in bbox, "You need to specify *.min.x");
    assert("max" in bbox, "You need to specify *.max.x");
    assert("x" in bbox["min"], "You need to specify *.max.x");

    // var this_ = IMPLICIT;  // an implicit argument

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

        vresampl: {iters: this_.use_II? 1:0, c: 1.0},
        projection: {enabled: this_.use_II? 1:0},
        qem: {enabled: (this_.use_II && this_.use_II_qem)?1:0},

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

        vresampl: {iters: 1, c: 0.4},
        projection: {enabled: 1},
        qem: {enabled: 1},
        //subdiv: {enabled: 1},
        //overall_repeats: 2,
        subdiv: {enabled: 0},
        overall_repeats: 1,

        debug: {
            enabled_pointsets: 0,
            post_subdiv_noise: 0.01,
        },
    };


    console.log("*********************");


    console.log (" mc properties : " + JSON.stringify(mc_properties));
    var mp5_str = JSON.stringify(shape_properties);
    var geom = IMPLICIT.make_geometry(mp5_str, mc_properties,
        function (verts, faces) {
            // ThreeJS-specific code

            // var ignore_root_matrix: Does not need other (MC-related) arguments.

            var allocate_buffers = true;
            var geom = new LiveBufferGeometry79(verts, faces, allocate_buffers);

            // Set the normals
            // var ignore_root_matrix = mc_params.ignore_root_matrix;  // Does not need other (MC-related) arguments.
            //geom.update_normals(this_, verts, mp5_str, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.
            //
            console.error(mp5_str);

            IMPLICIT.make_normals_into_geometry(geom, mp5_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

            //this_.aaaaaaaaa(verts);

            return geom;
        }
    );

    return geom;
}
