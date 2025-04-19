
// NodeJS code
// A CLI demo file to test the compiled Emscripten module on NodeJS.
// Usage: executed via: scripts/testing/sanity-test-1.bash
// Polygonises (meshify/vectorised renders(!)) a Cone.
// 3D geomeric content: Uses `example_objects.js` 's `provide_input()` 's DEFAULT_OBJ_SELECTOR, which is a cone.


const chai = require('chai');
const expect = chai.expect;

const [, , compiled_js_filename] = process.argv;
if (!compiled_js_filename) {
  throw new Error('Usage: node sanity1.js "./compiled-escripten-filename.js"');
}

// todo: remove comment
/*
const mcc = require(compiled_js_filename);
console.log(mcc);
// mcc.battery()
const {
  _get_pointset_ptr,
  _get_pointset_size,
  _build_geometry_u,
  _build_geometry,
  _get_f_size,
  _get_v_size,
  _get_v,
  _get_f,
  _get_v_ptr,
  _get_f_ptr,
  _finish_geometry,
  _about,
  _set_object,
  _unset_object,
  _set_x,
  _unset_x,
  _calculate_implicit_values,
  _get_values_ptr,
  _get_values_size,
  _calculate_implicit_gradients,
  _get_gradients_ptr,
  _get_gradients_size,
  _main,
} = mcc;
*/
/*
_main(); // native function `main` called before runtime initialization
*/

/*
var Module = mcc;
var result = Module.onRuntimeInitialized = () => {
    Module.ccall('_get_v', // name of C function
        null, // return type
        null, // argument types
        null // arguments
   );
}
*/
// todo: clean up
// const {Service1, wait_for_full_reload} = require('./service_l1');
// s1.set_object();

// todo: ...

/*
https://javascript.info/mixins

Mixing for multiple subsets of Module functions
  for l1
  for l2
  for l3
  for arrow_utils
*/

async function old_pattern_deprecated() {
  const mcc2 = await wait_for_full_reload(mcc);
  console.log('loaded');
  const s1 = new Service1(mcc2);
  s1.about();
  console.log('ok');
}

async function run2() {
  console.log('Runtime-loading of module:');
  const wait_for_full_reload = require('./service_l1');
  const Service1 = await wait_for_full_reload(compiled_js_filename);
  const s1 = new Service1();
  s1.about();
  console.log('(JS Module loaded.)\n')
  // "type":"sdf_3d"
  const example_objects = require('../../examples/js-lib/example_objects.js');

  console.log('\nExample 3D object: two pars: task-options and shape:')
  const {shape_json, polygonization_json} = example_objects.provide_input(0.0, 0, {}, {});
  console.log('3D object generated (above are logs and debug prints).');

  console.log('\noutput to be fed into ImpliSolid Polygoniser:');
  console.log({shape_json, polygonization_json})


  const {
    _on_cpp_loaded,
  } = require('../../js_iteration_2/implisolid_main.js');

  function assert(cond, message) {
    if (!cond) {
        message = message || "Assertion failed for unspecified reason";
        console.error(message);
        console.error(message.stack);
        throw new Error("assert ", message);
    }
  }

  console.log('\nInterface: javascript object `IMPLICIT` (has 3 layers):');
  console.log('Runtime-loading of module, again:');
  const IMPLICIT = _on_cpp_loaded(Service1.emscriptenModule);
  console.log(IMPLICIT);
  console.log(IMPLICIT.about());
  console.log(IMPLICIT.service2);



  // Two polygonisation tests:
  // build_geometry() versus make_geometry():
  // make_geometry wraps around build_geometry

  console.log("\n1. make_geometry:");
  // should not have dependency on threejs. IMPLICIT needs to be generatd separately from service2.
  const q1 = IMPLICIT.service2.make_geometry(shape_json, polygonization_json,
    (verts, faces, allocate_buffer)=>{
      console.log('made');
      console.log('async end.. todo: async or promise');
      console.log({verts, faces, allocate_buffer});
      chai.expect(verts).to.be.an.instanceof(Float32Array);
      chai.expect(faces).to.be.an.instanceof(Uint32Array);
      // chai.expect(verts).to.be.an.instanceof(TypedArray);
      console.log(typeof verts, typeof faces, typeof allocate_buffer);

      console.log('made');
      console.log('');
  }, 'qq');

  console.log('make_geometry returned:', q1);
  console.log();

  console.log("\n2. build_geometry:");
  const q2 = IMPLICIT.service2.service1.build_geometry(JSON.stringify(shape_json), JSON.stringify(polygonization_json));
  console.log('build_geometry returned:', q2);
  console.log();

}

run2();

