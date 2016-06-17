'use strict';

//Core contains the mesh (or changes to the mesh)
function make_geometry_core( verts, faces) {


    var vertexCount = verts.length/3;
    var facecount = faces.length/3;
    var indexCount = facecount*3;

    if(VERBOSE){
        console.log("vertexCount="+vertexCount+ "   ,  facecount=" + facecount+ "   facecount*3="+(facecount*3));
    }

    // buffers
    var indices = new ( indexCount > 65535 ? Uint32Array : Uint16Array )( indexCount );
    var vertices = new Float32Array( vertexCount * 3 );
    var normals = new Float32Array( vertexCount * 3 );
    var uvs = new Float32Array( vertexCount * 2 );

    // offset variables
    //var vertexBufferOffset = 0;
    //var uvBufferOffset = 0;
    //var indexBufferOffset = 0;
    //var numberOfVertices = 0;

    // group variables
    //var groupStart = 0;

    // build each side of the box geometry
    //output: vertices, normals, uvs, indices

    var nans_warnings = 0;
    //console.log("vertexCount " + vertexCount)
    var SCALE = 1.
    var x,y,z,d;
    for(var i=0; i < vertexCount; i++)
    {
        for(var di = 0; di < 3; di++)
            vertices[i*3+di] = verts[i*3+di] * SCALE;
        x = verts[i*3+0];
        y = verts[i*3+1];
        z = verts[i*3+2];
        //console.log(x+" , " + y + " , " + z);

        d = Math.sqrt(x*x+y*y+z*z);
        d=d+0.;
        //if(d==0) d=1.;
        if(isNaN(d)) {
            d=1.;x=1.;y=1.;z=1.;
            console.log("Warning");
        }
        //console.log("x y z"+x+" "+y+" "+z+"   / "+ d)

        var sgn = +1;
        normals[i*3 + 0] = x/d*sgn;
        normals[i*3 + 1] = y/d*sgn;
        normals[i*3 + 2] = z/d*sgn;

        if(isNaN(x/d))
            nans_warnings ++;

        var d2 = Math.sqrt(x*x+y*y);
        uvs[i*2+0] = x/d2;
        uvs[i*2+1] = y/d2;
    }

    for(var i=0; i < facecount; i++)
        for(var si=0; si<3; si++)
            indices[i*3+si] = faces[i*3+si];

    //console.log(verts.subarray(0,3*3*3));
    //console.log(indices.subarray(0,3*3*3));


    if(nans_warnings > 0)
        console.error("WARNING: NaN in vertices. "+nans_warnings/978+" out of "+vertexCount/978+ "  subtract:"+(vertexCount-nans_warnings)/978);

    return {
        indices: indices,
        vertices: vertices,
        normals: normals,
        uvs: uvs
    };
}

/** Simply creates a geometry . This is static and cannot be modified when displayed. Instantiate a new one and make a new THREE.Mesh() */
var
MyBufferGeometry77 = function ( verts, faces ) {

    //MyBufferGeometry77 = function ( verts, faces, width, height, depth, widthSegments, heightSegments, depthSegments ) {

    //console.log(faces);
    if(VERBOSE)
        console.log("MyBufferGeometry77");

    THREE.BufferGeometry.call( this );

    //console.log("raw v size " + verts.length);
    //console.log("vl " + verts.length/3);
    //console.log("fl " + faces.length/3);

    this.type = 'MyBufferGeometry77';

    this.parameters = {
        //width: width,
        //height: height,
        //depth: depth,
        //widthSegments: widthSegments,
        //heightSegments: heightSegments,
        //depthSegments: depthSegments
    };

    var mesh_core = make_geometry_core(verts, faces);

    var materialIndex = 0;
    this.addGroup( 0, mesh_core.indices.length, materialIndex ); //not sure about *3 . Why??

    //this.addGroup( groupStart, groupCount, materialIndex );  //groupCount is same as indices' index.

    //modified, but not output: indexBufferOffset, vertexBufferOffset, uvBufferOffset, numberOfVertices, groupStart

    // build geometry
    this.setIndex( new THREE.BufferAttribute( mesh_core.indices, 3 ) );
    this.addAttribute( 'position', new THREE.BufferAttribute( mesh_core.vertices, 3 ) );
    this.addAttribute( 'normal', new THREE.BufferAttribute( mesh_core.normals, 3 ) );
    this.addAttribute( 'uv', new THREE.BufferAttribute( mesh_core.uvs, 2 ) );

};

MyBufferGeometry77.prototype = Object.create( THREE.BufferGeometry.prototype );
MyBufferGeometry77.prototype.constructor = MyBufferGeometry77;



function WGeometry77(verts, faces) {
    //vects, faces

    THREE.Geometry.call( this );

    this.type = 'ImplicitGeometry'; //?

    this.parameters = {
        /*width: width,
        height: height,
        depth: depth,
        widthSegments: widthSegments,
        heightSegments: heightSegments,
        depthSegments: depthSegments*/
    };


    //faces = faces.subarray(0, 300);
    this.fromBufferGeometry( new MyBufferGeometry77( verts, faces ) );
    //this.mergeVertices();

};

WGeometry77.prototype = Object.create( THREE.Geometry.prototype );
WGeometry77.prototype.constructor = WGeometry77;

