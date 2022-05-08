'use strict';

// LiveBufferGeometry71  MyBufferGeometry77 LiveBufferGeometry79

/** @param emscrModuleName The compiled .js file name */
function wait_for_full_reload(emscrModuleName) {
  // A possibility is to allow emscrModuleName to be the `Module`
  let Module;
  if (typeof emscrModuleName == 'string') {
    Module = require(emscrModuleName);
  } else {
    throw new Error('specify the emscripten-compiled .js filename: ' + emscrModuleName);
  }
  // called after require()
  return new Promise( (acc, rej)=> {
    Module.onRuntimeInitialized = () => {
     Service1.emscriptenModule = Module;
     acc(Service1);
    }
  });
}

// see implisolid_main.js
// getInterfaceFucntions, assign_functoins
function assign_emscripten_functions(service, Module) {

    // Low level API v1.0.0
    const low_level_api = {
        build_geometry: Module.cwrap('build_geometry', null, [ 'string', 'string']),
        get_v_size: Module.cwrap('get_v_size', 'number', []),
        get_f_size: Module.cwrap('get_f_size', 'number', []),
        get_v: Module.cwrap('get_v', null, ['number']),
        get_f: Module.cwrap('get_f', null, ['number']),
        get_v_ptr: Module.cwrap('get_v_ptr', 'number', []),
        get_f_ptr: Module.cwrap('get_f_ptr', 'number', []),
        finish_geometry: Module.cwrap('finish_geometry', null, []),

        set_object: Module.cwrap('set_object', 'number', ['string', 'boolean']),
        unset_object: Module.cwrap('unset_object', 'number', ['number']),
        // sends the x points for evaluation of implicit or gradient:
        set_x: Module.cwrap('set_x', 'number', ['number', 'number']),
        unset_x: Module.cwrap('unset_x', null, []),
        calculate_implicit_values: Module.cwrap('calculate_implicit_values', null, []),
        get_values_ptr: Module.cwrap('get_values_ptr', 'number', []),
        get_values_size: Module.cwrap('get_values_size', 'number', []),
        // boolean argument to normalize and reverse the vevtors, suitable for rendering:
        calculate_implicit_gradients: Module.cwrap('calculate_implicit_gradients', null, ['number']),
        get_gradients_ptr: Module.cwrap('get_gradients_ptr', 'number', []),
        get_gradients_size: Module.cwrap('get_gradients_size', 'number', []),

        get_pointset_ptr: Module.cwrap('get_pointset_ptr', 'number', ['string']),
        get_pointset_size: Module.cwrap('get_pointset_size', 'number', ['string']),
        about: Module.cwrap('about', null, []),
    };
    Object.keys(low_level_api).forEach(funcname => {
      service[funcname] = low_level_api[funcname];
    });

    if (!Module["_about"]) {
        // Don;t continue, otherwise, it will be very difficult to track the cause of the issue.
        throw new Error('function `about()` is not `extern "c"`, or not specified in EXPORTED_FUNCTIONS compiler options');
    }

    return service;
}
/*
const ImplicitService = (function () {
  'use strict';
}
*/
/*
{

  service.init_ = function() {
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

/*
// todo:move
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
* /

service.init_();
}
*/

class Service1 {
  constructor() {
    // expect(Service1.emscriptenModule).to.not.be.undefined;

    //this.serviceObj = {};
    const serviceObj = this;
    assign_emscripten_functions(serviceObj, Service1.emscriptenModule);
  }
}


// wait_for_full_reload(Module, callback)

// module.exports = { Service1, wait_for_full_reload };
module.exports = wait_for_full_reload;
