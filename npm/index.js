'use strict';
console.log("IINDEXX000.JS");
/*
module.exports = function() {
    console.log("IINDEXX.JS22");

    console.log("Not ready. Please wait while we finalise the npm package.");
    loadScript('./js_iteration_2/js/worker_api.js');
    //importScripts('../../build/mcc2.compiled.js');
    loadScript('./js_iteration_2/js/implisolid_main.js');
    console.log("ImpliSolid");
    return IMPLICIT;
};
*/

/*
var load_imp = function() {
    console.log("IINDEXX.JS22");

    console.log("Not ready. Please wait while we finalise the npm package.");
    //loadScript('../js_iteration_2/js/worker_api.js');
    var jj;
    jj = require('../js_iteration_2/js/worker_api.js');
    //importScripts('../../build/mcc2.compiled.js');
    loadScript('./js_iteration_2/js/implisolid_main.js');
    console.log("ImpliSolid");
    return ;
};
module.exports = load_imp();
//console.log("jj", jj);   // global
*/

//var jj = require('../js_iteration_2/js/worker_api.js');
var Module={};
Module.locateFile = function(){return '../build/mcc2.compiled.js.mem';};

var module=require('../build/mcc2.compiled.js');
var qq=require('./js_iteration_2/js/implisolid_main.js');
module.exports = [jj,module,qq];

