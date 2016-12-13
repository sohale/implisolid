'use strict';

/*
Three ways:
1- Using WGeometry77  using  faster=true
2- " " " using faster=false which uses  make_geometry_core()
3- " " " using faster=false which uses  make_geometry_core_slower()
4- UsingUpdatable geometry; WGeometry77.prototype.update()
*/

function make_geometry_core( verts, faces) {

    var indices = faces;
    var vertices = verts;
    var normals = null, uvs = null;
    //vertices = verts;
    //indices = faces;

    return {
        indices: indices,
        vertices: vertices,
        normals: normals,
        uvs: uvs
    };
}

//Core contains the mesh (or changes to the mesh)
function make_geometry_core_slower( verts, faces) {

    const ENABLE_NORMALS = false;

    var vertexCount = verts.length/3;
    var facecount = faces.length/3;
    var indexCount = facecount*3;

    if(VERBOSE){
        console.log("vertexCount="+vertexCount+ "   ,  facecount=" + facecount+ "   facecount*3="+(facecount*3));
    }

    // buffers
    var indices = new ( indexCount > 65535 ? Uint32Array : Uint16Array )( indexCount );
    var vertices = new Float32Array( vertexCount * 3 );

    var normals, uvs;  // plan
    if(ENABLE_NORMALS){
        normals = new Float32Array( vertexCount * 3 );
        uvs = new Float32Array( vertexCount * 2 );
    }
    else{
        normals = null;
        uvs = null;
    }

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

        if(ENABLE_NORMALS){
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


       function checktype(src, _type){
            //if(typeof src !== type_name)
            if(!(src instanceof _type))
             {
                console.error();
                throw "Developer's error; " + (typeof src) + (Object.prototype.toString.call(src.buffer));
            }
        }

        // See  https://developer.mozilla.org/en-US/docs/Web/JavaScript/Typed_arrays

        function copy_Float32Array(src)  {
            checktype(src, Float32Array);
            return new Float32Array(src);
        }
        function copy_Uint32Array(src)  {
            checktype(src, Uint32Array);
            return new Uint32Array(src);
        }
        function copy_Uint16Array(src)  {
            checktype(src, Uint16Array);
            return new Uint16Array(src);
        }

        function copy_Float32Array_preallocated(src, prealloc_size, randomize)  {
            checktype(src, Float32Array);
            var TYPE_SIZE = 4;
            var len_bytes = Math.max(prealloc_size*TYPE_SIZE, src.byteLength);
            var dst = new ArrayBuffer(len_bytes);
            var r = new Float32Array(dst);
            r.set(new Float32Array(src));
            if(randomize){
                for(var i=0*src.byteLength/TYPE_SIZE;i<r.length;i++){
                    r[i] = (Math.random()-0.5);
                }
            }
            return r;
        }
        function copy_Uint32Array_preallocated(src, prealloc_size)  {
            checktype(src, Uint32Array);
            var TYPE_SIZE = 4;
            var len_bytes = Math.max(prealloc_size*TYPE_SIZE, src.byteLength);
            var dst = new ArrayBuffer(len_bytes);
            var r = new Uint32Array(dst);
            r.set(new Uint32Array(src));
            return r;
        }
        function copy_Uint16Array_preallocated(src, prealloc_size)  {
            checktype(src, Uint16Array);
            var TYPE_SIZE = 2;
            var len_bytes = Math.max(prealloc_size*TYPE_SIZE, src.byteLength);;
            var dst = new ArrayBuffer(len_bytes);
            var r = new Uint16Array(dst);
            r.set(new Uint16Array(src));
            return r;
        }

/** Simply creates a geometry . This is static and cannot be modified when displayed. Instantiate a new one and make a new THREE.Mesh() */
function MyBufferGeometry77( verts, faces,  re_allocate) {

    THREE.BufferGeometry.call( this );
    this.type = 'MyBufferGeometry77';

    this.parameters = { };

    this.initbg = function(){
        const faster = true;
        if(faster){

            if(1){
                var materialIndex = 0;
                this.addGroup( 0, faces.length*1-10, materialIndex ); //not sure about *3 . Why??
                //console.log("ok2");
            }

            if(re_allocate){
                //Float32Array, Uint32Array
                /*
                faces = copy_Uint32Array(faces, prealloc_size);
                verts = copy_Float32Array(verts, prealloc_size);
                */
                console.log("Aloocating separate space for verts,faces.");
                faces = copy_Uint32Array_preallocated(faces, 30000*3);
                verts = copy_Float32Array_preallocated(verts, 30000*3);
            }
            //WRONG! WHEN   re_allocate is false

            // build geometry
            this.addAttribute( 'position', new THREE.BufferAttribute( verts, 3 ) );
            this.setIndex( new THREE.BufferAttribute( faces, 3 ) ); //new Uint32Array(faces) ??
            for(var j=0;j<faces.length;j++)
                if(this.index.array[j] !== faces[j]){
                    console.error(j);break;
                }
            //my_assert(this.index.array === faces);

            if(re_allocate){
                console.log("Aloocating separate space for norms, colors.");
                var normals = copy_Float32Array_preallocated(verts, 30000*3/10000, true);
                var colors = copy_Float32Array_preallocated(verts, 30000*3/10000, true);
                //var uvs = copy_Float32Array_preallocated(new Float32Array([]), 30000*3 * 0, true);

            }

            this.addAttribute( 'normal', new THREE.BufferAttribute( normals, 3, true ) );
            this.addAttribute( 'color', new THREE.BufferAttribute( colors, 3, true ) );
            //this.addAttribute( 'uv', new THREE.BufferAttribute( uvs, 2 ) );

        }else{
            var mesh_core = make_geometry_core(verts, faces);
            //var mesh_core = make_geometry_core_slower(verts, faces);

            if(0){
                var materialIndex = 0;
                this.addGroup( 0, mesh_core.indices.length, materialIndex ); //not sure about *3 . Why??
            }

            // build geometry
            this.setIndex( new THREE.BufferAttribute( mesh_core.indices, 3 ) );
            this.addAttribute( 'position', new THREE.BufferAttribute( mesh_core.vertices, 3 ) );
            if(mesh_core.normals) {
                this.addAttribute( 'normal', new THREE.BufferAttribute( mesh_core.normals, 3 ) );
                this.addAttribute( 'uv', new THREE.BufferAttribute( mesh_core.uvs, 2 ) );
            }
        }
    }

    // ThreeJS does not use prototype-based OOP.

    this.update_geometry = function(newPolygonizer) {
        var geometry = this;
        var nverts = newPolygonizer.get_v_size();
        var nfaces = newPolygonizer.get_f_size();
        var verts_address = newPolygonizer.get_v_ptr();
        var faces_address = newPolygonizer.get_f_ptr();
        var verts = Module.HEAPF32.subarray(
            verts_address/_FLOAT_SIZE,
            verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(
            faces_address/_INT_SIZE,
            faces_address/_INT_SIZE + 3*nfaces);

        // REFACTORED. NOTE TESTED.
        geometry.update_geometry1 = function(verts, faces);
    };

    this.update_geometry1 = function(verts, faces) {
        //NOT TESTED.


        //var g_nverts = geometry.attributes.position.count/3;  // Displayed size
        //check allocated space
        var g_nverts = geometry.attributes.position.array.length/3;  // Physical space size.
        var g_nfaces = geometry.index.array.length/3;

        //console.log(nverts); console.log(g_nverts);
        //console.log(nfaces); console.log(g_nfaces);

        //Bug revealed by commenting the following! Changing the geometry broke the C++ code !
        var nv3 = Math.min(nverts, g_nverts) * 3;
        var nf3 = Math.min(nfaces, g_nfaces) * 3;
        //console.log(nv3/3. +" === "+ verts.length/3.)
        my_assert(nv3 === verts.length);
        my_assert(nf3 === faces.length);
        /*
        for(var i=0;i<nv3;i++){
            geometry.attributes.position.array[i] = verts[i];  //or use copy
        }
        for(var i=0;i<nf3;i++){
            geometry.index.array[ i ] = faces[i];
        }
        */
        /*
            var geometryAttributes = geometry.attributes;
            var geometryAttribute = geometryAttributes[ name ];
            geometryAttribute instanceof THREE.InterleavedBufferAttribute
            if ( geometryAttribute instanceof THREE.InterleavedBufferAttribute ) {

            var data = geometryAttribute.data;
            if ( data instanceof THREE.InstancedInterleavedBuffer ) {
        */
        if(0){
            name = 'position';
            var geometryAttributes = geometry.attributes;
            var geometryAttribute = geometryAttributes[ name ];
            var a_IL = geometryAttribute instanceof THREE.InterleavedBufferAttribute
            if(a_IL){
                var data = geometryAttribute.data;
                var i_IL = data instanceof THREE.InstancedInterleavedBuffer;
                console.log( a_IL + " " + i_IL );
            }
            else{
                var instanced = geometryAttribute instanceof THREE.InstancedBufferAttribute;
                console.log( a_IL + " not interleaved " + instanced );  //false, false
            }
        }
        if(1){
            // *************************************
            // * The following code works only when we use a single instance of the geometry, i.e. used in one Mesh.
            // * It will not work well if the same BufferGeometry is used in more than one Mesh
            // * So the object will look fine if we disable the wireframe mesh.
            // *************************************
            geometry.attributes.position.array.set(verts);
            geometry.index.array.set(faces);
        }

        geometry.setDrawRange( 0, nf3 );
        /*geometry.clearGroups();
        geometry.addGroup( 0, nf3, 0 );
        */

        geometry.attributes.position.needsUpdate = true;
        geometry.index.needsUpdate = true;

        //geometry.index.array[i] === faces[i]


        //geometry.setBou

        //geometry.computeFaceNormals();
        //geometry.computeBoundingSphere();

        if(0){
        for(var j=0;j<faces.length;j++)
            if(!(geometry.index.array[j] === faces[j])){
                console.error(j);break;
            }
        for(var j=0;j<verts.length;j++){
            if(!(geometry.attributes.position.array[j] === verts[j]))
                console.error(j);break;
        }
        }
    };

    this.initbg()
};

MyBufferGeometry77.prototype = Object.create( THREE.BufferGeometry.prototype );
MyBufferGeometry77.prototype.constructor = MyBufferGeometry77;


// This class is not actually necessary.
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


/* The following code will not work, but is in progress.
WGeometry77.prototype.update = function(verts, faces) {
    var geom = this;
    //for ( var vi = 0, l = geom.vertices.length; vi < l; vi ++ ) {
    //    var time = Math.random();
    //    geom.vertices[ vi ].y += 0.01*(Math.random()-0.5) ; // 0.1* 35 * Math.sin( vi / 5 + ( time + vi ) / 7 );
    //}
    var l = Math.min(geom.vertices.length, verts.length/3);
    for ( var vi = 0; vi < l; vi ++ ) {
        //geom.vertices[ vi ].x += verts[vi*3+0];
        //geom.vertices[ vi ].y += verts[vi*3+1];
        //geom.vertices[ vi ].z += verts[vi*3+2];
        geom.vertices[ vi ].x = verts[vi*3+0];
        geom.vertices[ vi ].y = verts[vi*3+1];
        geom.vertices[ vi ].z = verts[vi*3+2];
    }
    var lv = Math.min(geom.faces.length, faces.length/3);
    for ( var fi = 0; fi < lv; fi ++ ) {
        geom.faces[ fi ].a = faces[fi*3+0];
        geom.faces[ fi ].b = faces[fi*3+1];
        geom.faces[ fi ].c = faces[fi*3+2];
    }

    geom.verticesNeedUpdate = true;
}
*/

/*
function UpdatableGeometry77(verts, faces)
{
    WGeometry77.call( this, verts, faces );
    this.type = 'UpdatableGeometry77';

}

//this.geo = new WGeometry77(verts, faces);
UpdatableGeometry77.prototype.update = function(verts, faces) {
    geom = this;
    var time = math.random();
    for ( var vi = 0, l = geom.vertices.length; vi < l; vi ++ ) {
            geom.vertices[ vi ].y += 35 * Math.sin( vi / 5 + ( time + vi ) / 7 );
        }
    geom.verticesNeedUpdate = true;
}

UpdatableGeometry77.prototype = Object.create( WGeometry77.prototype );
UpdatableGeometry77.prototype.constructor = UpdatableGeometry77;
*/
