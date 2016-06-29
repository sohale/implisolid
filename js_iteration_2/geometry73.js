'use strict';




function __check_TypedArray_type(src, _type){
    //if(typeof src !== type_name)
    if(!(src instanceof _type))
     {
        console.error();
        throw "Developer's error; " + (typeof src) + (Object.prototype.toString.call(src.buffer));
    }
}

/**
* @param {Number | String} initial_value: if "random" then randomize color else if it's a number color set to this value
*/
function copy_Float32Array_preallocated(src, prealloc_size, initial_value)  {
    __check_TypedArray_type(src, Float32Array);
    var TYPE_SIZE = 4;
    var len_bytes = Math.max(prealloc_size*TYPE_SIZE, src.byteLength);
    var dst = new ArrayBuffer(len_bytes);
    var r = new Float32Array(dst);
    r.set(new Float32Array(src));
    if(initial_value === "random"){
        for(var i=0*src.byteLength/TYPE_SIZE;i<r.length;i++){
            r[i] = (Math.random()-0.5);
        }
    }else if (Number(initial_value) === initial_value){
        console.log("color number :" + initial_value + " " + len_bytes);
        for(var i=0*src.byteLength/TYPE_SIZE;i<r.length;i++){
            r[i] = initial_value;
        }

    }
    return r;
}
function copy_Uint32Array_preallocated(src, prealloc_size)  {
    __check_TypedArray_type(src, Uint32Array);
    var TYPE_SIZE = 4;
    var len_bytes = Math.max(prealloc_size*TYPE_SIZE, src.byteLength);
    var dst = new ArrayBuffer(len_bytes);
    var r = new Uint32Array(dst);
    r.set(new Uint32Array(src));
    return r;
}


/** Simply creates a geometry . This is static and cannot be modified when displayed. Instantiate a new one and make a new THREE.Mesh() */
function LiveBufferGeometry71( verts, faces,  pre_allocate) {

    THREE.BufferGeometry.call( this );
    this.type = 'LiveBufferGeometry71';

    this.parameters = { };

    this.allocate_buffers = function(){

        if(pre_allocate){
            console.log("Allocating separate space for verts,faces.");
            faces = copy_Uint32Array_preallocated(faces, 30000*3);
            verts = copy_Float32Array_preallocated(verts, 30000*3);
        }
        assert(pre_allocate);  // Will be wrong if pre_allocate is false

        // build geometry

        this.addAttribute( 'index', new THREE.BufferAttribute( faces, 3 ) );
        this.addAttribute( 'position', new THREE.BufferAttribute( verts, 3 ) );
        //this.setIndex( new THREE.BufferAttribute( faces, 3 ) ); //new Uint32Array(faces) ??
        for(var j=0;j<faces.length;j++)
            if(this.attributes.index.array[j] !== faces[j]){
                console.error(j);break;
            }
        //assert(this.attributes.index.array === faces);

        if(pre_allocate){
            console.log("Allocating separate space for norms, colors.");
            var normals = copy_Float32Array_preallocated(verts, 30000*3/10000, "random");
            //var colors = copy_Float32Array_preallocated(verts, 30000*3/10000, "random");
            //var uvs = copy_Float32Array_preallocated(new Float32Array([]), 30000*3 * 0, "random");
            normals.set(verts);
        }

        this.addAttribute( 'normal', new THREE.BufferAttribute( normals, 3, true ) );
        //this.addAttribute( 'color', new THREE.BufferAttribute( colors, 3, true ) ); color is overidden
        //this.addAttribute( 'uv', new THREE.BufferAttribute( uvs, 2 ) );

        //var materialIndex = 0;
        //this.addGroup( 0, faces.length*1-10, materialIndex );
    }

    // ThreeJS does not use prototype-based OOP.

    this.update_geometry = function(implicit_service) {

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;

        var geometry = this;
        var nverts = implicit_service.get_v_size();
        var nfaces = implicit_service.get_f_size();
        var verts_address = implicit_service.get_v_ptr();
        var faces_address = implicit_service.get_f_ptr();
        var verts = Module.HEAPF32.subarray(
            verts_address/_FLOAT_SIZE,
            verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(
            faces_address/_INT_SIZE,
            faces_address/_INT_SIZE + 3*nfaces);

        var g_nverts = geometry.attributes.position.array.length/3;  // Physical space size.
        var g_nfaces = geometry.attributes.index.array.length/3;

        var nv3 = Math.min(nverts, g_nverts) * 3;
        var nf3 = Math.min(nfaces, g_nfaces) * 3;
        assert(nv3 === verts.length);
        assert(nf3 === faces.length);

        // *************************************
        // * The following code works only when we use a single instance of the geometry, i.e. used in one Mesh.
        // * It will not work well if the same BufferGeometry is used in more than one Mesh
        // * So the object will look fine if we disable the wireframe mesh.
        // *************************************
        geometry.attributes.position.array.set(verts);
        geometry.attributes.index.array.set(faces);


        //geometry.setDrawRange( 0, nf3 );
        //geometry.attributes.index.length = nf3;
        var gl_chunkSize=21845;

        var ii = 0;
        geometry.offsets.push({start:ii, index: ii , count: Math.min(nf3 - ii, gl_chunkSize*3)});
        /*geometry.clearGroups();
        geometry.addGroup( 0, nf3, 0 );
        */

        //geometry.computeBoundingSphere();

        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.index.needsUpdate = true;

        /*
        g =currentMeshes[0].geometry

        IMPLICIT.finish_geometry();
        IMPLICIT.build_geometry(28,0.9);
        g.update_geometry(IMPLICIT)


        for(var i=0;i<1000;i++){ IMPLICIT.finish_geometry();IMPLICIT.build_geometry(28,i*0.1); g.update_geometry(IMPLICIT);}

        */

    };

    this.allocate_buffers()
};

LiveBufferGeometry71.prototype = Object.create( THREE.BufferGeometry.prototype );
LiveBufferGeometry71.prototype.constructor = LiveBufferGeometry71;

