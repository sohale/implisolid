
var w_impli3 = {};  // not worker side, but uses worker

w_impli3.custom_mc_settings = {};

// Copied from implisolid_main.js
var SET_ME_TO_TRUE = false;
w_impli3.getLiveGeometry  = function(dict, bbox, ignore_root_matrix, geom_callback) {
    // w_impli3.custom_mc_settings is a non-explicit argument changes the default settings
    var shape_properties = dict;
    var s = 1;

    my_assert(bbox, "You need to specify the bounding box");
    my_assert("min" in bbox, "You need to specify *.min.x");
    my_assert("max" in bbox, "You need to specify *.max.x");
    my_assert("x" in bbox["min"], "You need to specify *.max.x");

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

    // todo: Dependency-inject CONFIG.
    var CONFIG;
    var CONFIG_implisolid = CONFIG ? CONFIG.implisolid : {default_mc_resolution: 28, use_II: 1, repeats: 1,  use_III: 0,  use_II: 0};

    // Designer-specific
    var getResolution  = function(bb) {
        return CONFIG_implisolid.default_mc_resolution;
        const max_value = 40;
        const min_value = 14;
        const factor = CONFIG.implisolid.default_mc_resolution;
        var max_length = Math.max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"], bb["zmax"] - bb["zmin"]);
        var tmp =  Math.min(max_value,max_length*factor);
        return  Math.floor(Math.max(tmp,min_value));

        // 1 -> 28
        // 2 -> 48
    };

    var mc_res = CONFIG_implisolid.default_mc_resolution;
    var mc_properties = {
        resolution: getResolution(bb),
        box: bb,
        ignore_root_matrix: ignore_root_matrix,

        vresampl: CONFIG_implisolid.use_II? {iters: 1, c: 0.4} : {iters: 0, c: 1.0},
        projection: {enabled: CONFIG_implisolid.use_II? 1 : 0},
        qem: {enabled: CONFIG_implisolid.use_II_qem? 1 : 0},
        //subdiv: {enabled: 1},
        //overall_repeats: 2,
        subdiv: {enabled: CONFIG_implisolid.use_III? 1 : 0},
        overall_repeats: CONFIG_implisolid.repeats,

        debug: {
            enabled_pointsets: 0,
            post_subdiv_noise: CONFIG_implisolid.use_noise? 0.01 : 0.0,
        },
    };

    // w_impli3.custom_mc_settings is a non-explicit argument changes the default settings
    if (w_impli3.custom_mc_settings) {
        // nonrecursive is enough, because we want to replace everything except for "box".
        // note: this is not tested.
        // Example usage:
        //    IMPLICIT.custom_mc_settings = {vresampl: {iters: 1, c: 1} };
        mc_properties = merge_dicts_nonrecursive(mc_properties, w_impli3.custom_mc_settings);
    }

    var mp5_json_str = JSON.stringify(shape_properties);
    var polygonization_properties_json = JSON.stringify(mc_properties);
// service2.make_geometry
    //var geom =   // no geom returned here anymore
    wapi_make_geometry(mp5_json_str, polygonization_properties_json,
        function (vf_dict) {
            var verts = vf_dict.verts;
            var faces = vf_dict.faces;
            // ThreeJS-specific code

            // var ignore_root_matrix: Does not need other (MC-related) arguments.

            var allocate_buffers = true;
            var geom = new LiveBufferGeometry79(verts, faces, allocate_buffers);

            // Set the normals
            // var ignore_root_matrix = mc_params.ignore_root_matrix;  // Does not need other (MC-related) arguments.
            //geom.update_normals(this_, verts, mp5_json_str, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.
            //
            if (!mp5_json_str)
                console.error(mp5_json_str);
            if (SET_ME_TO_TRUE)
            service3.make_normals_into_geometry(geom, mp5_json_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

            //this_.aaaaaaaaaA(verts);

            // no geom anymore
            // return geom;
            var shape_id = 99; // vf_dict.shape_id
            geom_callback(geom, shape_id);
        }
    );

    // return geom;
};

