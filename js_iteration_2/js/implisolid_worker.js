
var w_impli3 = {};  // not worker side, but uses worker

w_impli3.custom_mc_settings = {};

// Copied from implisolid_main.js
var SET_ME_TO_TRUE = false;
/*
w_impli3.getLiveGeometry  = function(shape_properties, bbox, ignore_root_matrix, geom_callback) {
    var mc_properties = IMPLICIT.make_polygonization_settings(bbox, ignore_root_matrix);

    var shape_json_str = JSON.stringify(shape_properties);
    var polygonization_properties_json = JSON.stringify(mc_properties);
    this.getLiveGeometry_from_json(shape_json_str, polygonization_properties_json, geom_callback);
}
*/
w_impli3.getLiveGeometry_from_json  = function(shape_json_str, polygonization_setttings_json_str, geom_callback) {
    // w_impli3.custom_mc_settings is a non-explicit argument changes the default settings

    assert(typeof shape_json_str === "string");
    assert(typeof polygonization_setttings_json_str === "string");
    assert(typeof geom_callback === "function");

    // service2.make_geometry
    //var geom =   // no geom returned here anymore
    wapi_make_geometry(shape_json_str, polygonization_setttings_json_str,
        function (vf_dict) {
            var verts = vf_dict.verts;
            var faces = vf_dict.faces;
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
            service3.make_normals_into_geometry(geom, shape_json_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

            //this_.aaaaaaaaaA(verts);

            // no geom anymore
            // return geom;
            var shape_id = 99; // vf_dict.shape_id
            geom_callback(geom, shape_id);
        }
    );

    // return geom;
};

