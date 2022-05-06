
// NodeJS code

const [, , compiled_js_filename] = process.argv;
if (!compiled_js_filename) {
  throw new Error('Usage: node sanity1.js "./compiled-escripten-filename.js"');
}

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

async function run1() {
  /*
  wait_for_full_reload(mcc)
  .then(mcc =>{
    s1 = new Service1(mcc);
    s1.about();
  })
  */
  const mcc2 = await wait_for_full_reload(mcc);
  console.log('loaded');
  const s1 = new Service1(mcc);
  console.log('about:');
  s1.about();
  console.log('ok');
}

async function run2() {
  const wait_for_full_reload = require('./service_l1');
  const Service1 = await wait_for_full_reload(compiled_js_filename);
  console.log('loaded');
  const s1 = new Service1();
  console.log('about:');
  s1.about();
  console.log('ok');

  // "type":"sdf_3d"
  const example_objects = require('../../js_iteration_2/examples/js/example_objects.js');
  const {shape_json, polygonization_json} = example_objects.provide_input(0.0, 0, {}, {});
  console.log('xx11')
  console.log({shape_json, polygonization_json})
  console.log('222')

  global.Module = Service1.emscriptenModule;
  const {
    _on_cpp_loaded,
  } = require('../../js_iteration_2/implisolid_main.js');

  console.log('1')
  IMPLICIT = _on_cpp_loaded(); // implicit argument: global.Module
  console.log(IMPLICIT);
  console.log(IMPLICIT.about());
  console.log(IMPLICIT.service2);
  // should not have dependency on threejs. IMPLICIT needs to be generatd separately from service2.
  const q = IMPLICIT.service2.service1.build_geometry(JSON.stringify(shape_json), JSON.stringify(polygonization_json));
  console.log(q);

  //console.log(IMPLICIT);
}

run2();

