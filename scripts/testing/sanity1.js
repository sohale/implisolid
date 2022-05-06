
// NodeJS code

console.log(process.argv);
const [, , compiled_js_filename] = process.argv;
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
  const {Service1, wait_for_full_reload} = require('./service_l1');
  const mcc2 = await wait_for_full_reload(compiled_js_filename);
  console.log('loaded');
  const s1 = new Service1(mcc2);
  console.log('about:');
  s1.about();
  console.log('ok');
}

run2();

const example_objects = require('../../js_iteration_2/examples/js/example_objects.js');
