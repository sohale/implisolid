'use strict';

/*
See:
  - docs/implisolid-build/demo1/js/example_objects.js
*/

// import checkMP5Object

// based on MP5_GENERIC_EXAMPLE_MOON
var _EXAMPLE_0 = {
  "printerSettings":{
       "name":"test",  "layerThickness":0.2, "emptyLayer":0, "infillSpace":4, "topThickness":0.6, "paramSamples":75,
       "speedRate":1000, "circleSpeedRate":1000, "temperature":220, "inAirSpeed":7000, "flowRate":0.035,
       "criticalLength":35, "retractionSpeed":2400, "retractionLength":5, "shellNumber":3, "material":"PLA 2.85mm",
       "autoZScar":true, "zScarGap":0.5, "critLayerTime":6, "filamentDiameter":2.85 },
  "mp5-version":"0.3",
  "root":{
     "type":"root",
     "children":[
        //
     ]
  }
};

const _make_sphere = function(ellipsoid_radius) {
  return {
    type: "ellipsoid",
    //radius: subjective_time,
    matrix:[
        ellipsoid_radius, 0,0, 0,
        0,ellipsoid_radius, 0, 0,
        0,0,ellipsoid_radius ,  0,
        0,0,0,   1]
  };
}

const PRIMITIVES = {
  sphere: _make_sphere,
};

class Mp5builder {
    constructor() {
      this.mp5 = _EXAMPLE_0;
      // this.cursor = 0;
      this.root_cursor = this.mp5.root.children;
      this.cursor = this.root_cursor; // make it a function?
    }

    primitives = PRIMITIVES

    sphere() {
      this.cursor.push(this.primitives.sphere(1.0));
      checkMP5Object(this.mp5);
      return this;
    }
}

var mp5builder = new Mp5builder();


/*
(()=> {
    mp5builder.sphere = function () {
      ;
    }
})();
*/

// move to a differnt file:
// polygonization_json
function create_polygonization_json() {

  // see docs/implisolid-build/demo1/js/example_objects.js
      //  `if (DONT_CHANGE || is_update_mode == 0) {`
  const x0 = 0, y0 = 0, z0 = 0;
  const BB_SIZE = 15; // 18
  const bbox = {xmin: x0-BB_SIZE, xmax: x0+BB_SIZE, ymin: y0-BB_SIZE , ymax: y0+BB_SIZE, zmin: z0-BB_SIZE, zmax: z0+BB_SIZE};

  const REPEATS = 0; // has no effect?

  const mc_properties_json = {
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
  return mc_properties_json;
}

if (typeof module !== 'undefined') {
  module.exports = {mp5builder, create_polygonization_json};
}
