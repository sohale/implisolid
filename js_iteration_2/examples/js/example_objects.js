'use strict';

function setShapeMatrix2Eye(shape_matrix, sz) {
    if (typeof sz == "undefined")
        sz = 1.0;

    for(var ii=0; ii<16; ii++) {
        shape_matrix[0] = 0;
    }

    shape_matrix[0] = 1 * sz;
    shape_matrix[5] = 1 * sz;
    shape_matrix[10] = 1 * sz;
    shape_matrix[15] = 1;

    shape_matrix[3] = 0;
    shape_matrix[7] = 0;
    shape_matrix[11] = 0;
}


var MP5_GENERIC_EXAMPLE_MOON = {
   "printerSettings":{
        "name":"test",  "layerThickness":0.2, "emptyLayer":0, "infillSpace":4, "topThickness":0.6, "paramSamples":75,
        "speedRate":1000, "circleSpeedRate":1000, "temperature":220, "inAirSpeed":7000, "flowRate":0.035,
        "criticalLength":35, "retractionSpeed":2400, "retractionLength":5, "shellNumber":3, "material":"PLA 2.85mm",
        "autoZScar":true, "zScarGap":0.5, "critLayerTime":6, "filamentDiameter":2.85 },
   "mp5-version":"0.3",
   "root":{
      "type":"root",
      "children":[
         {
            "type":"Difference",
            "protected":false,
            "children":[
               {
                  "type":"cylinder",
                  "displayColor":{ "x":0.77, "y":0.04, "z":0.18 },
                  "matrix":[ 42.5, 0, 0, 0, 0, 49.6, 0, 0, 0, 0, 10, 0, 0, 0, 0, 1 ],
                  "index":652818
               },
               {
                  "type":"Difference",
                  "protected":false,
                  "children":[
                     {
                        "type":"cylinder",
                        "displayColor":{"x":0.82, "y":0.66, "z":0.74 },
                        "matrix":[10, 0, 0, 0, 0, 10, 0, 0, 0, 0, 10, 0, 0, 0, 0, 1],
                        "index":1272174
                     },
                     {
                        "type":"cylinder",
                        "displayColor":{"x":0.11, "y":0.076,"z":0.63 },
                        "matrix":[10, 0, 0, 0.66, 0, 10, 0, 6.22, 0, 0, 10, 0.00001, 0, 0, 0, 1 ],
                        "index":2463576
                     }
                  ],
                  "displayColor":{ "x":0.66, "y":0.45, "z":0.72 },
                  "matrix":[2.38, 0, 0, 0.36, 0, 2.38, 0, 0.56, 0, 0, 2.38, 6.9, 0, 0, 0, 1],
                  "index":413872
               }
            ],
            "displayColor":{"x":0.55, "y":0.07, "z":0.12 },
            "matrix":[ 1, 0, 0, 0.33, 0, 1, 0, 0.156, 0, 0, 1, 0.000000000000014, 0, 0, 0, 1],
            "index":6565922
         }
      ]
   }
};
var sz0 = 100.0;
var HEART = '{"printerSettings":{"name":"test","layerThickness":0.2,"emptyLayer":0,"infillSpace":4,"topThickness":0.6,"paramSamples":75,"speedRate":1000,"circleSpeedRate":1000,"temperature":220,"inAirSpeed":7000,"flowRate":0.035,"criticalLength":35,"retractionSpeed":2400,"retractionLength":5,"shellNumber":3,"material":"PLA 2.85mm","autoZScar":true,"zScarGap":0.5,"critLayerTime":6,"filamentDiameter":2.85},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"Union","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.9661750055896976,"y":0.8085857086202395,"z":0.41578037212168595},"matrix":[10,0,0,-3.7368777450135866,0,10,0,-1.9559832356144682,0,0,10,1.7323194345664206e-7,0,0,0,1],"index":7575510},{"type":"cube","displayColor":{"x":0.23399071378141634,"y":0.31584816496653323,"z":0.35457351563365425},"matrix":[10,0,0,1.867397834493545,0,10,0,-1.7325527119763677,0,0,10,-9.734106853898084e-10,0,0,0,1],"index":5587759},{"type":"cylinder","displayColor":{"x":0.43814645627496795,"y":0.39556472441055845,"z":0.3415798414286939},"matrix":[10,0,0,1.8694799105200275,0,10,0,3.688535947590836,0,0,10,-1.7225853365943067e-7,0,0,0,1],"index":6657333}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.6470588235294118,"y":0.2784313725490196,"z":0.5882352941176471},"matrix":[10.0,0,0,82.63768850593796,0,10.0,0,126.37324151118989,0,0,10.0,5.000000079354265,0,0,0,1],"index":4872526}]}}';
HEART = JSON.parse(HEART);
setShapeMatrix2Eye(HEART.root.children[0].matrix, sz0);
HEART = JSON.stringify(HEART);

var MOON = '{"printerSettings":{"name":"test","layerThickness":0.2,"emptyLayer":0,"infillSpace":4,"topThickness":0.6,"paramSamples":75,"speedRate":1000,"circleSpeedRate":1000,"temperature":220,"inAirSpeed":7000,"flowRate":0.035,"criticalLength":35,"retractionSpeed":2400,"retractionLength":5,"shellNumber":3,"material":"PLA 2.85mm","autoZScar":true,"zScarGap":0.5,"critLayerTime":6,"filamentDiameter":2.85},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"Difference","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.7675997200783986,"y":0.03892568708507049,"z":0.1754374135888661},"matrix":[35,0,0,0,0,35,0,0,0,0,9,0,0,0,0,1],"index":652818},{"type":"Difference","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.8122645344236872,"y":0.657334404743416,"z":0.7357336310755096},"matrix":[10,0,0,0,0,10,0,0,0,0,10,0,0,0,0,1],"index":1272174},{"type":"cylinder","displayColor":{"x":0.11421729990684737,"y":0.07562705374348999,"z":0.6324600862122098},"matrix":[10,0,0,0.658889604636343,0,10,0,6.215549332615993,0,0,10,1.3327027659215673e-7,0,0,0,1],"index":2463576}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.6627450980392157,"y":0.4549019607843137,"z":0.7215686274509804},"matrix":[2.381193509886417,0,0,0.3600215429489424,0,2.381193509886417,0,0.5604901669421452,0,0,2.381193509886417,6.9059681360437395,0,0,0,1],"index":413872}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.5529411764705883,"y":0.06666666666666667,"z":0.11764705882352941},"matrix":[1,0,0,0.32938436512727,0,1,0,0.15604124684634,0,0,1,0.000000000000014,0,0,0,1],"index":6565922}]}}';

var SIMPLE_CONE = '{"printerSettings":{},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"icone","displayColor":{"x":0.7,"y":0.7,"z":0.7},"matrix":[8,0,0,0,0,8,0,0,0,0,8,0,0,0,0,1],"index":9185154}]}}';






function provide_input (subjective_time, is_update_mode, globals) {
    // is_update_mode: 0 => init
    // is_update_mode: 1 => update
    // is_update_mode: 2 => deprecated (used in make_geometry_old1() )

    // asked for:
    var shape_json, polygonization_json;
    const DONT_CHANGE = false;  // Whether use the mode=0 for mode 1. Should be false.

    if (DONT_CHANGE || is_update_mode == 0) {
        // var subjective_time = subjective_time;

        // var ellipsoid_radius = 8.0 * 0.1; // subjective_time;
        var ellipsoid_radius = 8.0; // subjective_time;
        var sphere_dict1 = {
            type: "ellipsoid",
            //radius: subjective_time,
            matrix:[
                ellipsoid_radius, 0,0, 0,
                0,ellipsoid_radius, 0, 0,
                0,0,ellipsoid_radius * 0.5 ,  0,
                0,0,0,   1]
        };
        var sphere_dict = {
            type: "Union",
            matrix:[
                1, 0,0, 0,
                0,1, 0, 0,
                0,0,1,  0,
                0,0,0,   1],
            children: [sphere_dict1, sphere_dict1]
        };
        //globals.global_mp5 = mp5_json_0;
        var mp5_json_sphere_dict= JSON.parse(JSON.stringify(MP5_GENERIC_EXAMPLE_MOON));
        mp5_json_sphere_dict.root.children=[sphere_dict];
        var SPHERE = JSON.stringify(mp5_json_sphere_dict);
        // console.log(SPHERE);

        var tetrahedron_dict = {
            type: "tetrahedron",
            corners: [[0, 0, 5],[0, 10, 0],[10, 0, 0],[0, 0, -5]],
            matrix:[
                1, 0,0, 0,
                0,1, 0, 0,
                0,0,1,  0,
                0,0,0,   1]
        };
        var mp5_json_tetrahedron_dict= JSON.parse(JSON.stringify(MP5_GENERIC_EXAMPLE_MOON));
        mp5_json_tetrahedron_dict.root.children=[tetrahedron_dict];
        var TETRAHEDRON = JSON.stringify(mp5_json_tetrahedron_dict);

        // Select the object to display

        // var mp5_json = HEART;     const BB_SIZE = 9;
        // var mp5_json = MOON;         const BB_SIZE = 9+9;
        // var mp5_json = TETRAHEDRON;  const BB_SIZE = 9;
        // var mp5_json = SPHERE;  const BB_SIZE = 9;
        var mp5_json = SIMPLE_CONE; const BB_SIZE = 12.0/10.0;

        var shape_dict = JSON.parse(mp5_json).root.children[0];

        if(is_update_mode == 0) {
            var sz = 1.;
        } else {
            var sz = subjective_time;
        }
        /*
        HEART.root.children[0].matrix[0] = sz;
        HEART.root.children[0].matrix[5] = sz;
        HEART.root.children[0].matrix[10] = sz;
        */

        // Scale the object
        shape_dict.matrix[0] = sz;
        shape_dict.matrix[5] = sz;
        shape_dict.matrix[10] = sz;

        /*
            // Move the object
            shape_dict.matrix[3] = 0;
            shape_dict.matrix[7] = 0;
            shape_dict.matrix[11] = 0;
        */

        if(false) {
            setShapeMatrix2Eye(shape_dict.matrix);
        }

        //var bb = getBoundingBoxForTree();

        mp5_json = JSON.stringify(shape_dict);

        globals.global_mp5 = mp5_json;
        console.log(mp5_json);



        var x0=0, y0=0, z0=0;
        //var BB_SIZE = 18;
        //const BB_SIZE = 9;
        var bbox = {xmin: x0-BB_SIZE, xmax: x0+BB_SIZE, ymin: y0-BB_SIZE , ymax: y0+BB_SIZE, zmin: z0-BB_SIZE, zmax: z0+BB_SIZE};

        // create new geometry
        var mc_properties_json = JSON.stringify({
            resolution: Math.floor(14),
            box: bbox,
            ignore_root_matrix: false,

            vresampl: {iters: 3, c: 1.0},
            projection: {enabled: 1},
            qem: {enabled: 1},
            subdiv: {enabled: 1},
            overall_repeats: 3,
            debug: {
                enabled_pointsets: 0,
                // only_rank
                post_subdiv_noise: 0.0,
            },
            // bug: When vresampl.iters is large, the centroid projection projects into 0 (0,0,0)
        });

        return {shape_json: mp5_json, polygonization_json: mc_properties_json};
    }

    //Alternatively, use a growing sphere only.
    if(is_update_mode == 1 || is_update_mode == 2) {

        var x0=0, y0=0, z0=0;
        const BB_SIZE = 9 /10.0;
        var bbox = {xmin: x0-BB_SIZE, xmax: x0+BB_SIZE, ymin: y0-BB_SIZE , ymax: y0+BB_SIZE, zmin: z0-BB_SIZE, zmax: z0+BB_SIZE};
        // console.error(bbox);

        var subjective_time = subjective_time * 5;

        // update
        var mc_properties_json__update = {
            resolution: 28,
            box: bbox,
            ignore_root_matrix: false,

            // makes things very slow
            vresampl: {iters: 1, c: 1.0},
            projection: {enabled: 1},
            qem: {enabled: 1},

            debug: {
                post_subdiv_noise: 0.001 / 3.0,
            },
        };
        //implisolid_.build_geometry( 28, mc_properties_json__update, "sphere", subjective_time);
        /*
        var mp5_json = JSON.stringify({
            type: "ellipsoid",
            //radius: subjective_time, //boudbl_mushroom only
            matrix:[
                subjective_time, 0,0, x0,
                0,subjective_time, 0, y0,
                0,0, subjective_time, z0,
                0,0,0,   1]
        });

        shape_json = mp5_json;
        */

        var d = JSON.parse(globals.global_mp5);
        var updated_size = subjective_time / 10.;
        d.matrix[0] = updated_size;
        d.matrix[5] = updated_size;
        d.matrix[10] = updated_size;
        shape_json = d;
        // shape_json = globals.global_mp5;
        // polygonization_json = mc_properties_json__update;

        // update
        var simple_MC_only_polygonization_json = {
            resolution: Math.floor( 14 ),
            box: bbox,
            ignore_root_matrix: false,
            vresampl: {iters: 0, c: 1.0},
            projection: {enabled: 0},
            qem: {enabled: 0},
            subdiv: {enabled: 0},
            overall_repeats: 0,

            debug: {
                enabled_pointsets: 0,
                // only_rank
                post_subdiv_noise: 0.0,
            },
        };
        //polygonization_json = simple_MC_only_polygonization_json;
        //globals.global_mp5 is not updated here. It is the original object before the update.
        //var tuple = [shape_json, polygonization_json];

        //return {shape_json: shape_json, polygonization_json: simple_MC_only_polygonization_json};
        return {shape_json: shape_json, polygonization_json: mc_properties_json__update};
    }

}
