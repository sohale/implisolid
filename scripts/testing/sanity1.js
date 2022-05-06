
// NodeJS code

console.log(process.argv);
const [, , compiled_js_file] = process.argv;
const mcc = require(compiled_js_file);
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

_main(); // native function `main` called before runtime initialization
