
function grid_bed(step) {
    my_assert(step > 0.0 && "Grid step has to be non zero");
    // var step = 25;
    var y_floor =  step * -3;  // -75;
    var ngrids = 40;
    var span = step * ngrids;  // 1000;
    //note: span/step == ngrids
    var line_material = new THREE.LineBasicMaterial( { color: 0x303030 } );
    var _geometry = new THREE.Geometry();
    for ( var i = 0; i <= ngrids; i++ ) {

        var v = i*step-span/2;

        _geometry.vertices.push( new THREE.Vector3( -span/2,   y_floor,    v ) );
        _geometry.vertices.push( new THREE.Vector3( +span/2,    y_floor,    v ) );

        _geometry.vertices.push( new THREE.Vector3( v, y_floor, -span/2 ) );
        _geometry.vertices.push( new THREE.Vector3( v, y_floor, +span/2 ) );
    }
    if(threejs_rev == THREEJS_R71) {
        throw new Error("Grid: to be implemented on r71- where THREE.LineSegments() does not exist");
    }
    var grid_lines = new THREE.LineSegments( _geometry, line_material );
    return grid_lines;
}


