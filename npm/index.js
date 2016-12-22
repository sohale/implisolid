'use strict';
module.exports = function() {
    console.log("Not ready. Please wait while we finalise the npm package.");
    loadScript('./js_iteration_2/js/worker_api.js');
    //importScripts('../../build/mcc2.compiled.js');
    loadScript('./js_iteration_2/js/implisolid_main.js');
    console.log("ImpliSolid");
    return IMPLICIT;
};

