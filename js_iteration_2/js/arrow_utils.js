'use strict';

// todo: refactor imlisolid function calls

//*************** Arrow Functions ****************************
/*
function make_arrow(){
    // Create an arrow and add it to the scene.
    var dir  = new THREE.Vector3( 10, 4, 0);
    var origin = new THREE.Vector3(0, 0, 0);
    var length = 10;
    var hex = 0xffff00;

    var arrowHelper = new THREE.ArrowHelper(dir, origin, length, hex);
    scene.add(arrowHelper);
    console.log("Arrow");
    return arrowHelper;
}
*/

function make_multiple_arrows(n) {
    /**
     * Create an arrow and add it to the scene.
     */
    var arrowHelpers = [];
    for(var i=0;i<n;i++) {
        var dir  = new THREE.Vector3( 10, 4, 0);
        var origin = new THREE.Vector3(0, 0, 0);
        var length = 10;
        var hex = 0xffff00;

        var arrowHelper = new THREE.ArrowHelper(dir, origin, length, hex);
        scene.add(arrowHelper);
        arrowHelpers.push(arrowHelper);
        console.log("Arrow");
    }
    return arrowHelpers;
}

function update_arrows(arrows, implisolid_, mp5_str, n, ignore_root_matrix) {
    /**
     * Update the arrow parameters
    */

    var vlist = new Float32Array(n*3);
    for(var i=0;i<n;i++) {
        var x = MESH_SCALE*10 * (Math.random()*2 - 1.);
        var y = MESH_SCALE*10 * (Math.random()*2 - 1.);
        var z = MESH_SCALE*10 * (Math.random()*2 - 1.);
        /*vlist.push(x);
        vlist.push(y);
        vlist.push(z);*/
        vlist[i*3+0] = x;
        vlist[i*3+1] = y;
        vlist[i*3+2] = z;
    }
    //var _verts = new Float32Array([x, y, z]);
    //var _verts = new Float32Array(vlist);
    var _verts = vlist;
    my_assert(_verts instanceof Float32Array);

var result1 =
    implisolid_.set_object(mp5_str, ignore_root_matrix);
console.log("SET_OBJECT() RESULTS:", result1);
    implisolid_.set_vect(_verts);
    implisolid_.calculate_implicit_gradients();
    const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
    var ptr = implisolid_.get_gradients_ptr();
    var ptr_len = implisolid_.get_gradients_size();
    var g = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
    //console.log("grad len = " +  ptr_len+ "  grad = " + g);  // x 4

    var dir = null;
    for(var i=0;i<n;i++) {

        //if not bush
        var x = _verts[i*3 + 0];
        var y = _verts[i*3 + 1];
        var z = _verts[i*3 + 2];

        var gx = g[i*3 + 0];
        var gy = g[i*3 + 1];
        var gz = g[i*3 + 2];
        //implisolid_.unset_x();
        //implisolid_.unset_object(1);


        //console.log("*****************************************")
        //console.log(arrow.position);
        //console.log("*****************************************")
        //arrow.setLength(200*Math.random());
        var l = Math.sqrt(gx*gx + gy*gy + gz*gz);
        //l = l * 100/10000;
        //console.log("l: "+l);
        var ls = MESH_SCALE*1.5;
        arrows[i].position.set(x, y, z);
        arrows[i].setLength(ls);
        // This can easily get broken. i.e. by multiplying ls.
        if (!dir) {
            dir = new THREE.Vector3( 0, 0, 0);
        }
        dir.set(-gx/l, -gy/l, -gz/l);
        arrows[i].setDirection(dir);
        //console.log("change_arrow");
    }

    implisolid_.unset_x();
var result1 =
    implisolid_.unset_object(1);
console.log("UNSET_OBJ() RESULTS:", result1);
}
