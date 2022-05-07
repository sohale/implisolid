'use strict';
/*
    todo: move
        from js_iteration_2/examples/simple_viewer1/siviewer1.js
        to js_iteration_2/examples/js-lib/siviewer1.js
        or js_iteration_2/examples/js-lib/sv/siviewer1.js

    todo: move: js-lib, lib-js
        js_iteration_2/examples/js
        to js_iteration_2/examples/js-lib
    ✓ (done)

    movnig:
        js_iteration_2/examples
        to ./examples

    todo: refactor with elements in js_iteration_1/mcc2_3js_r79-sohailver.html
        into common files in `js-lib`
*/

const CONF = {
  NAV: true,  // ENABLE_NAVIGATION
};
freeze(CONF);

function getCanvasDOM(document, id) {
  return document.getElementById(id);
}


class SimpleImplicitViewer {
    constructor() {
        // requirements:
        assert (THREE.REVISION >= 79, 'You need to use threejs r79+');
        CONF.NAV && assert (THREE.OrbitControls, 'You need to use ornbtconfitl');
    }
    bindCanvas(canvasDOMId) {
      this.canvasDOM = getCanvasDOM(canvasDOMId);
    }
    setObject(mp5Obj) {
      console.log(mp5Obj);
    }
    setPolygonizingSettings(tesselSettings) {
      console.log(tesselSettings);
    }
}
