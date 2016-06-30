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
function copy_Uint16Array_preallocated(src, prealloc_size)  {
    __check_TypedArray_type(src, Uint16Array);
    var TYPE_SIZE = 2;
    var len_bytes = Math.max(prealloc_size*TYPE_SIZE, src.byteLength);;
    var dst = new ArrayBuffer(len_bytes);
    var r = new Uint16Array(dst);
    r.set(new Uint16Array(src));
    return r;
}

/** Simply creates a geometry . This is static and cannot be modified when displayed. Instantiate a new one and make a new THREE.Mesh() */
function LiveBufferGeometry71( verts, faces,  pre_allocate, faces_capacity, verts_capacity) {

    THREE.BufferGeometry.call( this );
    this.type = 'LiveBufferGeometry71';

    this.parameters = { };

    if(faces_capacity === undefined) faces_capacity = 9000*0;
    if(verts_capacity === undefined) verts_capacity = 8000*0;
    //_expect(faces_capacity);
    //_expect(verts_capacity);
    console.log("verts : "+ verts_capacity + " faces : "+ faces_capacity);

    this.allocate_buffers = function(){
        var padded_faces, padded_verts;
        if(pre_allocate){
            console.log("Allocating separate space for verts,faces.");
            padded_faces = copy_Uint32Array_preallocated(faces, faces_capacity*3);
            padded_verts = copy_Float32Array_preallocated(verts, verts_capacity*3);
        }
        else
        {
            // Warning: not cloned
            //padded_faces = faces;
            //padded_verts = verts;
            console.error("padded_faces = ?");
        }
        assert(pre_allocate);  // Will be wrong if pre_allocate is false

        // build geometry

        this.addAttribute( 'index', new THREE.BufferAttribute( padded_faces, 3 ) );
        this.addAttribute( 'position', new THREE.BufferAttribute( padded_verts, 3 ) );
        //this.setIndex( new THREE.BufferAttribute( padded_faces, 3 ) ); //new Uint32Array(padded_faces) ??
        for(var j=0;j<faces.length;j++)
            if(this.attributes.index.array[j] !== faces[j]){
                console.error(j);break;
            }
        //assert(this.attributes.index.array === padded_faces);

        var padded_normals;
        if(pre_allocate){
            console.log("Allocating separate space for norms, colors.");
            padded_normals = copy_Float32Array_preallocated(verts, verts_capacity*3, "random");
            padded_normals.set(padded_verts);
            //var padded_colors = copy_Float32Array_preallocated(verts, verts_capacity*3, "random");
            //var uvs = copy_Float32Array_preallocated(new Float32Array([]), verts_capacity*3 * 0, "random");
            //padded_colors.set(padded_verts);
        }
        else{
            console.error("padded_normals = ?");
        }

        this.addAttribute( 'normal', new THREE.BufferAttribute( padded_normals, 3, true ) );
        //this.addAttribute( 'color', new THREE.BufferAttribute( padded_colors, 3, true ) ); //color is overidden
        //this.addAttribute( 'uv', new THREE.BufferAttribute( uvs, 2 ) );

        //var materialIndex = 0;
        //this.addGroup( 0, faces.length*1-10, materialIndex );

        var nf3 = faces.length;
        var gl_chunkSize=21845;

        /*
        var ii = 0;
        this.offsets.push({start:ii, index: ii , count: Math.min(nf3 - ii, gl_chunkSize*3)});
        */
        var triangles = nf3/3;
        this.offsets = [];
        var offsets = triangles / gl_chunkSize;
        for ( var i = 0; i < offsets; i ++ ) {

            var offset = {
                start: i * gl_chunkSize * 3,
                index: i * gl_chunkSize * 3,
                count: Math.min( triangles - ( i * gl_chunkSize ), gl_chunkSize ) * 3
            };

            this.offsets.push( offset );

        }
    }

    // ThreeJS does not use prototype-based OOP.
    /*
    This function may return new LiveBufferGeometry71 object, in this case the old geometry has to be replaced (by the new returned object) in the caller function.
    */
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
        //assert(nv3 === verts.length);
        //assert(nf3 === faces.length);

        // *************************************
        // * The following code works only when we use a single instance of the geometry, i.e. used in one Mesh.
        // * It will not work well if the same BufferGeometry is used in more than one Mesh
        // * So the object will look fine if we disable the wireframe mesh.
        // *************************************
        var grow_needed = false;
        var availableVertsSize = geometry.attributes.position.array.length;
        if(verts.length > availableVertsSize){

            grow_needed = true;
            
            //console.log("needed "+verts.length );
/*
            console.log(geometry.attributes.position.array.length);  */

            //geometry.attributes.position.array.set(verts.subarray(0, availableVertsSize));


/*
            var reserved_size = (verts.length/3) * 2;
            console.log("reserved "+reserved_size);
            var new_pos_array = copy_Float32Array_preallocated(verts, reserved_size*3);
            geometry.attributes.position.array = new_pos_array;
            //this.addAttribute( 'position', new THREE.BufferAttribute( new_pos_array, 3 ) );
            geometry.attributes.position.array.set(verts);


            var new_pos_array = copy_Float32Array_preallocated(verts, reserved_size*3);
            geometry.attributes.normal.array = new_pos_array;
            geometry.attributes.normal.array.set(verts);

            //var new_pos_array = copy_Float32Array_preallocated(verts, reserved_size*3);
            //geometry.attributes.color = new_pos_array;
            //geometry.attributes.color.set(verts);

            //var positionAttr = geometry.attributes.position;
            //positionAttr.updateRange.offset = 0; // where to start updating
            //positionAttr.updateRange.count = 1;

            //The problem is geometry.attributes.position.buffer is not updated.


            console.log(geometry.attributes.position.array.length);
            */
        }else{
            geometry.attributes.position.array.set(verts);
        }

        var availableFacesSize = geometry.attributes.index.array.length;
        if(faces.length > availableFacesSize){
            //geometry.attributes.index.array.set(faces.subarray(0, availableFacesSize));
            grow_needed = true;
        }else{
            geometry.attributes.index.array.set(faces);
        }

        if(grow_needed){
            console.log("increasing capacity");
            //this.dispose();

            return new LiveBufferGeometry71( verts, faces,  true, Math.max(availableFacesSize/3, faces.length/3) * 1.5 + 1, Math.max(availableVertsSize/3, verts.length/3) * 1.5 +1);
        }
        var copied_faces = Math.min(faces.length, availableFacesSize);

        //copied_faces = nf3;

        //geometry.setDrawRange( 0, nf3 );
        //geometry.attributes.index.length = nf3;

        /*
        var ii = 0;
        var gl_chunkSize=21845;
        geometry.offsets = [];
        geometry.offsets.push({start:ii, index: ii , count: Math.min( copied_faces - ii, gl_chunkSize*3)});
        //FIXME: when many vertices
        _expect(copied_faces - ii <= gl_chunkSize*3);
        */

        var triangles = copied_faces/3;
        var gl_chunkSize=21845;
        geometry.offsets = [];
        var offsets = triangles / gl_chunkSize;
        for ( var i = 0; i < offsets; i ++ ) {

            var offset = {
                start: i * gl_chunkSize * 3,
                index: i * gl_chunkSize * 3,
                count: Math.min( triangles - ( i * gl_chunkSize ), gl_chunkSize ) * 3
            };

            geometry.offsets.push( offset );

        }
        //geometry.computeBoundingSphere();


        /*geometry.clearGroups();
        geometry.addGroup( 0, nf3, 0 );
        */



        //geometry.computeOffsets();
        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.index.needsUpdate = true;
        geometry.attributes.normal.needsUpdate = true;
        //geometry.attributes.color.needsUpdate = true;

        /*
        g =currentMeshes[0].geometry

        IMPLICIT.finish_geometry();
        IMPLICIT.build_geometry(28,0.9);
        g.update_geometry(IMPLICIT)


        for(var i=0;i<1000;i++){ IMPLICIT.finish_geometry();IMPLICIT.build_geometry(28,i*0.1); g.update_geometry(IMPLICIT);}

        */
        return false;




    };

    this.allocate_buffers()
};

LiveBufferGeometry71.prototype = Object.create( THREE.BufferGeometry.prototype );
LiveBufferGeometry71.prototype.constructor = LiveBufferGeometry71;




function test_update(t){
    var g = currentMeshes[0].geometry;

    IMPLICIT.finish_geometry();
    IMPLICIT.needsFinish = false;

    var new_geometry = IMPLICIT.build_geometry(28, t);
    IMPLICIT.needsFinish = true;

    if(new_geometry)
        currentMeshes[0].geometry = new_geometry
    
    g.update_geometry(IMPLICIT);
}

/*
var t=0;var iid=setInterval(function(){test_update(t+=0.02);},6);
*/
