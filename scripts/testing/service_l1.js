'use strict';

// LiveBufferGeometry71  MyBufferGeometry77 LiveBufferGeometry79





// function wait_for_full_reload(Module, callback) {

function wait_for_full_reload(moduleOrModuleName) {
  let Module;
  if (typeof moduleOrModuleName == 'string') {
    Module = require(moduleOrModuleName);
  } else {
    // not recommended
    //Module = moduleOrModuleName;
    throw new Error('specify the emscripten-compiled .js filename: ' + moduleOrModuleName);
  }
  // called after require()
  return new Promise( (acc, rej)=> {
    /*var result =*/
    Module.onRuntimeInitialized = () => {
      /*
      Module.ccall('myFunction', // name of C function 
          null, // return type
          null, // argument types
          null // arguments
     );*/
     // callback(Module);

     Service1.emscriptenModule = Module;
     acc(Service1);
    }
  });
}

// see implisolid_main.js
// getInterfaceFucntions, assign_functoins
function assign_emscripten_functions(service, Module) {
    'use strict';

    // Low level API
    service.build_geometry = Module.cwrap('build_geometry', null, [ 'string', 'string']);
    service.get_v_size = Module.cwrap('get_v_size', 'number', []);
    service.get_f_size = Module.cwrap('get_f_size', 'number', []);
    service.get_v = Module.cwrap('get_v', null, ['number']);
    service.get_f = Module.cwrap('get_f', null, ['number']);
    service.get_v_ptr = Module.cwrap('get_v_ptr', 'number', []);
    service.get_f_ptr = Module.cwrap('get_f_ptr', 'number', []);
    service.finish_geometry = Module.cwrap('finish_geometry', null, []);

    service.set_object = Module.cwrap('set_object', 'number', ['string', 'boolean']);
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
        service.about = "C++ method problem: no .about()";
        console.error("C++ method problem: a correct version not found");
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
