'use strict';


// DRY 6 folds!
function clone_array_float(verts) {
    var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * verts.length);
    var g = new Float32Array(raw_buffer);
    g.subarray(0, verts.length).set(new Float32Array(verts));
    return g;
}

function concatenate_array_float(verts1, verts2) {
    var n1 = verts1.length;
    var n2 = verts2.length;
    var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * (n1 + n2) );
    var verts12 = new Float32Array(raw_buffer);
    verts12.subarray(0, n1).set(new Float32Array(verts1));
    verts12.subarray(n1, (n1+n2)).set(new Float32Array(verts2));  // deliberate bug for TDD
    return verts12;
}

function make_indices_for_n1_n2(n) {
    var n1 = n;
    var raw_buffer = new ArrayBuffer( Uint32Array.BYTES_PER_ELEMENT * (n1) * 3 );
    var indices12 = new Uint32Array(raw_buffer);
    for (var i=0; i < n1; i++) {
        indices12[i * 3] = i;
        indices12[i * 3 + 1] = i + n1; //i + n1;
        indices12[i * 3 + 2] = i; //i + n1;
    }
    my_assert(indices12.length == n1 * 3);
    return indices12;
}

/*
function clone_array_int(indices) {
    var raw_buffer = new ArrayBuffer( Uint32Array.BYTES_PER_ELEMENT * indices.length);
    var g = new Uint32Array(raw_buffer);
    g.subarray(0, indices.length).set(new Uint32Array(indices));
    return g;
}
*/


