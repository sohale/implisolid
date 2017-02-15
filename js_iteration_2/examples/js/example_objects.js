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

var SIMPLE_SCREW = '{"printerSettings":{},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"screw","matrix":[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],"v": [0,2,0],"pitch": 0.5,"profile":"sin","delta_ratio":1.5,"end_type": "0","index":9185154}]}}';
// union of ball and cone
// var SIMPLE_SCREW = '{"printerSettings":{"PRINTER":"Ultimaker Origin","FILAMENT":"PLA","DEFAULT":0},"mp5-version":"0.4","root":{"type":"root","children":[{"type":"Union","protected":false,"children":[{"type":"icone","displayColor":[0,0,1],"matrix":[2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -1, 0, 0, 1],"index":1769767},{"type":"iellipsoid","displayColor":[0,0,1],"matrix":[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],"index":117621}],"displayColor":[0.8705882352941177,0.4196078431372549,0.8705882352941177],"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],"index":6389872}]},"createdAt":"2017-01-05T17:47:01.921Z","title":"Unnamed (2017-01-05T17-52-32)","contributors":[],"unique_id":"197bba0a-1e1e-45f0-8aba-2796307ba05c","licence":{}}'
// cone subtracts two cones
// var SIMPLE_SCREW = '{"printerSettings":{"PRINTER":"Ultimaker Origin","FILAMENT":"PLA","DEFAULT":0},"mp5-version":"0.4","root":{"type":"root","children":[{"type":"Difference","protected":false,"children":[{"type":"Difference","protected":false,"children":[{"type":"icone","displayColor":[1,0,0],"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],"index":4423185},{"type":"half_plane","plane_vector":[0,0,1],"plane_point":[0,0,0.2],"displayColor":[0.0057981841651877,0.76608476471992,0.025145201935641],"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],"index":4532637}],"displayColor":[0.50980392156863,0.031372549019608,0.37647058823529],"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],"index":2937721},{"type":"half_plane","plane_vector":[0,0,-1],"plane_point":[0,0,-0.2],"displayColor":[0.0057981841651877,0.76608476471992,0.025145201935641],"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],"index":5610308}],"displayColor":[0.66274509803922,0.45490196078431,0.72156862745098],"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],"index":6578438}]},"createdAt":"2017-01-05T17:47:01.921Z","title":"Unnamed (2017-01-05T17-52-32)","contributors":[],"unique_id":"197bba0a-1e1e-45f0-8aba-2796307ba05c","licence":{}}'
// cone top_bottom_lid
// var SIMPLE_SCREW = '{"printerSettings":{},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"top_bottom_lid","matrix":[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],"v": [0,2,0],"pitch": 0.5,"profile":"sin","delta_ratio":1.5,"end_type": "0","index":9185154}]}}';

var asmjs = '{"printerSettings":{"PRINTER":"Ultimaker Origin","FILAMENT":"PLA","DEFAULT":0},"mp5-version":"0.4","root":{"type":"root","children":[{"type":"sdf_3d","displayColor":[0.114,0.075,0.63],"param1":0.5,"matrix":[10,0,0,0,0,10,0,0,0,0,10,0,0,0,0,1],"index":2463577,"implicit":"_f =  (+param1)*(+param1) - ((_x/0.7) *(_x/0.7) + _y * _y + _z * _z);","gradient":"_gx = (-2) * _x /(0.7); _gy = -2 * _y; _gz = -2 * _z;"}]},"createdAt":"2016-12-27","title":"asmjs example 2","description":"Multiple colours","contributors":["sohail"],"unique_id":"3ef4a9f2-6175-4a9d-b421-720eac674d89","licence":{}}';


//test for extrusion
var SIMPLE_EXTRUSION = '{"printerSettings":{},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"extrusion","matrix":[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 1],"size": 6,"end_type": "0","index":9185154}]}}';


 // {
 //                "type": "screw",
 //                "axis": [
 //                    [0,  0,  0.3165],    // A
 //                    [0.5,  0,  -0.5],     // A+w*len
 //                ],
 //                "start_orientation": [0,1,0],   // the u vector

 //                "pitch": 2,   //  2mm distance between two dents
 //                "profile" : "sin+",        // "sin+" -> max{sin(t),0},   "sin"  -> sin(t)   profile of the indentation (thread profile)

 //                "diameter_inner":  12,   // 12mm : all sizes in millimeters
 //                "diameter_outer":  8,    // 8mm

 //                "end_type": "0",   // Now only one type "0", but there are seveeal types: "C", "R", "F"

 //                 // common in all mp5 objects:
 //                "index": 6125140,
 //                "displayColor": [ 0.71, 0.164, 0    ],
 //                 // note that there is no "matrix" field.
 //            }
 //        ]


/**
    @return {
       shape_json: ,  // JSONified (sting)
       polygonization_json:  // not JSONified (dict)
   }
*/
function provide_input (subjective_time, is_update_mode, globals) {
    // is_update_mode: 0 => init
    // is_update_mode: 1 => update
    // is_update_mode: 2 => deprecated (used in make_geometry_old1() )

    // asked for:
    var shape_json, polygonization_json;
    const DONT_CHANGE = false;  // Whether use the mode=0 for mode 1. Should be false.


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



        /*var screw_dict = {
            type:"screw",
            axis:[[0,0,0.3165],[0.5, 0, -0.5]],
            start_orientation: [0,1,0],
            pitch: 2,
            profile:"sin",
            diameter_inner:12.0,
            diameter_outer:8.0,
            end_type: "0"
        };*/
        var screw_dict = JSON.parse(SIMPLE_SCREW).root.children[0];


        var mp5_json_screw_dict= JSON.parse(JSON.stringify(MP5_GENERIC_EXAMPLE_MOON));
        mp5_json_screw_dict.root.children=[screw_dict];
        var SCREW = JSON.stringify(mp5_json_screw_dict);

        var extrusion_dict = JSON.parse(SIMPLE_EXTRUSION).root.children[0];


        var mp5_json_extrusion_dict= JSON.parse(JSON.stringify(MP5_GENERIC_EXAMPLE_MOON));
        mp5_json_extrusion_dict.root.children=[extrusion_dict];
        var EXTRUSION = JSON.stringify(mp5_json_extrusion_dict);
        // Select the object to display

        var used_matrix = true;
        // var mp5_json = HEART;     const BB_SIZE = 9;
        // var mp5_json = MOON;         const BB_SIZE = 9+9;
        // var mp5_json = TETRAHEDRON;  const BB_SIZE = 9;
        // var mp5_json = SPHERE;  const BB_SIZE = 9;
        // var mp5_json = SIMPLE_SCREW; var BB_SIZE = 3.0; used_matrix = false;
        // var shape_dict = JSON.parse(mp5_json).root.children[0];
        var mp5_json = null;

        // var obj_selector = "screw";
        //var obj_selector = "cone";
        // var obj_selector = "metaballs";
        //var obj_selector = "extrusion";
        var obj_selector = "asmjs";
        //var obj_selector = "metaballs";
        //var obj_selector = "extrusion";

        var resize_mp5 = function(){console.error("dont know how to resize.");}

        switch (obj_selector) {
            case "screw":
                console.log("OK SCREW")
                var mp5_json = SIMPLE_SCREW; var BB_SIZE = 0.5;var used_matrix = false;
                // var mp5_json = SCREW; var BB_SIZE = 12.0/10.0; used_matrix = false;
                var shape_dict = JSON.parse(mp5_json).root.children[0];

                if(is_update_mode == 1 || is_update_mode == 2) {
                    BB_SIZE = 9 /10.0;
                }

                resize_mp5 = function (d, sz) {
                    d.diameter_outer = sz;
                    d.diameter_inner = sz * 0.7;
                    var df = [0,0,0]; var dd=0;
                    for (var j = 0 ; j < 3 ; ++j) {
                        //d.axis[1][j] = sz + d.axis[0][j];
                        df[j] = d.axis[1][j] - d.axis[0][j];
                        dd += df[j] ** 2;
                    }
                    for (var j = 0 ; j < 3 ; ++j) {
                        d.axis[1][j] = sz * df[j] / Math.sqrt(dd) + d.axis[0][j];
                    }
                }

            break;
            case "cone":
                console.log("OK CONE", BB_SIZE)
                var mp5_json = SIMPLE_CONE; var BB_SIZE = 12.0/10.0; used_matrix = true;
                var shape_dict = JSON.parse(mp5_json).root.children[0];

                if(is_update_mode == 1 || is_update_mode == 2) {
                    BB_SIZE = 9 /10.0;
                }

                resize_mp5 = function (d, sz) {
                    d.matrix[0] = sz;
                    d.matrix[5] = sz;
                    d.matrix[10] = sz;
                    console.error("resizing", sz);
                }

            break;
            case "metaballs":

                var mp5_json = SIMPLE_CONE;
                var shape_dict = JSON.parse(mp5_json).root.children[0];
                shape_dict.type = "meta_balls";
                // shape_dict.time = subjective_time;  // 1 is not good.

                var BB_SIZE = 0.55; used_matrix = false;

                console.error(shape_dict);
                if(is_update_mode == 0) {
                    //shape_dict.time = 0.1;
                }

                if(is_update_mode == 1 || is_update_mode == 2) {
                }

                resize_mp5 = function (d, sz) {
                    for( var i=0; i<16; ++i) d.matrix[i] = 0;
                    d.matrix[0] = 1;
                    d.matrix[5] = 1;
                    d.matrix[10] = 1;
                    d.matrix[15] = 1;
                    d.time = sz;
                    console.error("no resizing. time=", sz);
                }

                resize_mp5(shape_dict, 0);  //is called later

            break;
            case "extrusion":
                console.log("OK EXTRUSION", BB_SIZE)
                var mp5_json = EXTRUSION; var BB_SIZE = 1.2; used_matrix = true;
                var shape_dict = JSON.parse(mp5_json).root.children[0];

                if(is_update_mode == 1 || is_update_mode == 2) {
                    BB_SIZE = 0.9;
                }

                resize_mp5 = function (d, sz) {
                    d.matrix[0] = sz;
                    d.matrix[5] = sz;
                    d.matrix[10] = sz;
                    console.error("resizing", sz);
                }
            break;
            case "asmjs":
                console.log("OK asmjs", BB_SIZE)
                var mp5_json = asmjs; var BB_SIZE = 3; used_matrix = true;
                var shape_dict = JSON.parse(mp5_json).root.children[0];


                // resize_mp5 = function (d, sz) {
                //     d.matrix[0] = sz;
                //     d.matrix[5] = sz;
                //     d.matrix[10] = sz;
                //     console.error("resizing", sz);
                // }


            break;
            default:
                console.error("error");
        }


        // is_update_mode == 0 always (here)
        if(is_update_mode == 0) {
            var sz = 1.;
        } else {
            var sz = subjective_time;  // this will be always zero
        }
        /*
        HEART.root.children[0].matrix[0] = sz;
        HEART.root.children[0].matrix[5] = sz;
        HEART.root.children[0].matrix[10] = sz;
        */

        // Scale the object
        if (used_matrix) {
            shape_dict.matrix[0] = sz;
            shape_dict.matrix[5] = sz;
            shape_dict.matrix[10] = sz;
        }


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

        //globals.global_mp5 = mp5_json;
        console.log(mp5_json);

        const simple_MC_only_polygonization_json = function(BSH) { return {
            resolution: Math.floor( 14 ),
            box: //bbox,
                {xmin: -BSH, xmax: +BSH, ymin: -BSH, ymax: +BSH, zmin: -BSH, zmax: +BSH},
            ignore_root_matrix: false,

            vresampl: {iters: 0, c: 1.0},

            projection: {enabled: 0},
            qem: {enabled: 0},
            subdiv: {enabled: 0},
            overall_repeats: 1,

            debug: {
                enabled_pointsets: 0,
                // only_rank
                post_subdiv_noise: 0.0,
            },
        }};


    ;
    if (DONT_CHANGE || is_update_mode == 0) {

        var x0=0, y0=0, z0=0;
        //var BB_SIZE = 18;
        const BB_SIZE = 15;
        var bbox = {xmin: x0-BB_SIZE, xmax: x0+BB_SIZE, ymin: y0-BB_SIZE , ymax: y0+BB_SIZE, zmin: z0-BB_SIZE, zmax: z0+BB_SIZE};

        var REPEATS = 0; // has no effect! ()?!)
        // create new geometry
        // tiger
        var mc_properties_json = {
            resolution: Math.floor(80),
            box: bbox,
            ignore_root_matrix: false,

            vresampl: {iters: 1, c: 1.0},
            projection: {enabled: 1},
            qem: {enabled: 1},
            subdiv: {enabled: 1},
            overall_repeats: REPEATS,
            debug: {
                enabled_pointsets: 0,
                // only_rank
                post_subdiv_noise: 0.0,
            },
            // bug: When vresampl.iters is large, the centroid projection projects into 0 (0,0,0)
        };

        // mc_properties_json = simple_MC_only_polygonization_json(0.55);
        console.log("mc_properties_json__update", mc_properties_json);
        // // watch: String.fromCharCode.apply(null, Module.HEAPU8.subarray(i2,i2+100) )
        //function ss(i2,n){return String.fromCharCode.apply(null, Module.HEAPU8.subarray(i2,i2+n));}
        return {shape_json: mp5_json, polygonization_json: mc_properties_json};
    }

    //Alternatively, use a growing sphere only.
    if(is_update_mode == 1 || is_update_mode == 2) {

        var x0=0, y0=0, z0=0;
        // const BB_SIZE = 9 /10.0;
        var bbox = {xmin: x0-BB_SIZE, xmax: x0+BB_SIZE, ymin: y0-BB_SIZE , ymax: y0+BB_SIZE, zmin: z0-BB_SIZE, zmax: z0+BB_SIZE};
        // console.error(bbox);

        //var subjective_time = subjective_time * 5;

        var REPEATS = 0;  // 3-> slow
        // update
        var mc_properties_json__update = {
            resolution: 28,
            box: bbox,
            ignore_root_matrix: false,

            // makes things very slow
            vresampl: {iters: 1, c: 1.0},
            projection: {enabled: 1},
            qem: {enabled: 1},

            overall_repeats: REPEATS,

            debug: {
                post_subdiv_noise: 0.001 / 3.0 * 0,
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

        var d = JSON.parse(mp5_json); //(globals.global_mp5);
        var updated_size = subjective_time / 10.;
        resize_mp5(d, updated_size);
        shape_json = d;
        // shape_json = globals.global_mp5;
        // polygonization_json = mc_properties_json__update;

        // update

        //globals.global_mp5 is not updated here. It is the original object before the update.

        //return {shape_json: shape_json, polygonization_json: simple_MC_only_polygonization_json(0.55)};
        return {shape_json: shape_json, polygonization_json: mc_properties_json__update};
    }

}
