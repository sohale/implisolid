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
function copy_Float32Array_preallocated(src, min_prealloc_size, initial_value)  {
    __check_TypedArray_type(src, Float32Array);
    var TYPE_SIZE = 4;
    if(min_prealloc_size % 1 !== 0) console.error("min_prealloc_size must be integer: " + min_prealloc_size);
    //if(min_prealloc_size !== int(min_prealloc_size)) console.error("min_prealloc_size must be integer: " + min_prealloc_size);
    var len_bytes = Math.max(min_prealloc_size*TYPE_SIZE, src.byteLength);
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
function copy_Uint32Array_preallocated(src, min_prealloc_size)  {
    __check_TypedArray_type(src, Uint32Array);
    var TYPE_SIZE = 4;
    //if(min_prealloc_size !== int(min_prealloc_size)) console.error("min_prealloc_size must be integer: " + min_prealloc_size);
    if(min_prealloc_size % 1 !== 0) console.error("min_prealloc_size must be integer: " + min_prealloc_size);
    var len_bytes = Math.max(min_prealloc_size*TYPE_SIZE, src.byteLength);
    var dst = new ArrayBuffer(len_bytes);
    // console.log("dst[0] : " + dst[0]);  // output:  undefined
    var r = new Uint32Array(dst);
    r.set(new Uint32Array(src));
    return r;
}
function copy_Uint16Array_preallocated(src, min_prealloc_size)  {
    __check_TypedArray_type(src, Uint16Array);
    var TYPE_SIZE = 2;
    if(min_prealloc_size !== int(min_prealloc_size)) console.error("min_prealloc_size must be integer: " + min_prealloc_size);
    var len_bytes = Math.max(min_prealloc_size*TYPE_SIZE, src.byteLength);;
    var dst = new ArrayBuffer(len_bytes);
    var r = new Uint16Array(dst);
    r.set(new Uint16Array(src));
    return r;
}

const threejs_r71 = false;
const threejs_r79 = true;
assert( threejs_r71 || threejs_r79 );

/** Simply creates a geometry . This is static and cannot be modified when displayed. Instantiate a new one and make a new THREE.Mesh() */
function LiveBufferGeometry79( verts_, faces_,  pre_allocate_, min_faces_capacity_, min_verts_capacity_) {
    // Alternative design: another argument: given_normals

    THREE.BufferGeometry.call( this );
    this.type = 'LiveBufferGeometry79';
    // Use THREE.InstancedBufferGeometry for geometry with multiple instances.

    this.parameters = { };

    if(min_faces_capacity_ === undefined) min_faces_capacity_ = 9000*0;
    if(min_verts_capacity_ === undefined) min_verts_capacity_ = 8000*0;

    const GROWTH_FACTOR = 1.5;
    const GROWTH_ADDITIONAL = 1;

    if (pre_allocate_ == "Deliberately Empty") {
        //  todo => set some flags, dont create this box, etc.
        // for now, it creates the box, which will be replaced with something else => it's fine.
        // But if I dont create anything, things might go wrong.
        // Also we should not allow update() to zero faces. i.e.
        this.deliberately_empty = true;
    } else {
        this.deliberately_empty = false;
    }
    pre_allocate_ = !!pre_allocate_;

    //_expect(min_faces_capacity_);
    //_expect(min_verts_capacity_);
    //console.log("verts : "+ min_verts_capacity_ + " faces : "+ min_faces_capacity_);

    this.allocate_buffers = function(verts, faces,  pre_allocate, faces_capacity, min_verts_capacity) {
        /*
            Allocates buffers, also copies. Can be used to update the geometry object.
        */

        if(faces.length == 0) {
            if (!this.deliberately_empty) {
                console.error("emptyimplicit");
            }
            var verts = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 ]);
            var faces = new Uint32Array([0,1,2, 0,2,3, 0,4,5, 0,5,1, 1,5,6, 1,6,2, 2,6,3, 3,6,7, 4,5,6, 5,6,7]);
        }

        /*if (faces_capacity == 0 || faces.length == 0 || verts.length == 0 ) {
            console.warn("faces_capacity == 0");
        }
        */
        if (faces.length == 0 || verts.length == 0 ) {
            if (this.deliberately_empty) {
                // ok
            } else {
                console.error("faces.length, verts.length == 0");
            }
        }
        if (this.deliberately_empty) {
            // assert certain things
        }

        var POINTS_PER_FACE = 3;

        // console.log("min capacity for verts : "+ min_verts_capacity + " faces : "+ faces_capacity);
        var padded_faces, padded_verts;
        //if(pre_allocate){
            // console.log("Allocating separate space for verts,faces.");
        padded_faces = copy_Uint32Array_preallocated(faces, faces_capacity * POINTS_PER_FACE);
        padded_verts = copy_Float32Array_preallocated(verts, min_verts_capacity*3);
        //}
        /*
        else
        {
            // Warning: not cloned
            //padded_faces = faces;
            //padded_verts = verts;
            console.error("padded_faces = ?");
        }
        */
        assert(pre_allocate);  // Will be wrong if pre_allocate is false

        // build geometry

        /***********************************
         * Adding the attributes: position
         ***********************************/

         // todo:
         //  if exists (i.e. it's an update) then: geometry.removeAttribute('position')
        assert(padded_verts instanceof Float32Array);
        this.addAttribute( 'position', new THREE.BufferAttribute( padded_verts, 3 ) );  // does this replace the old one?

        /***********************************
         * Adding the attributes: index (i.e. faces)
         ***********************************/

        assert(padded_faces instanceof Uint32Array);
        _expect(padded_faces.length < 65535);
        if(threejs_r71) {
            if (POINTS_PER_FACE !== 3) {
                console.error("In ThreeJS r71 POINTS_PER_FACE should be (& can be) 3.");
            }
            this.addAttribute( 'index', new THREE.BufferAttribute( padded_faces, POINTS_PER_FACE ) );
        }
        if(threejs_r79) {
            // If you choose 3, it will not choose them in .setDrawRange() correctly.
            this.setIndex( new THREE.BufferAttribute( padded_faces, 1 ) );  // Why 1??
        }

        //this.setIndex( new THREE.BufferAttribute( padded_faces, POINTS_PER_FACE ) ); //new Uint32Array(padded_faces) ??
        for(var j=0;j<faces.length;j++)
            //if(this.attributes.index.array[j] !== faces[j]){  //threejs_r71
            if(this.index.array[j] !== faces[j]){  //r77
                console.error(j);break;
            }
        //assert(this.attributes.index.array === padded_faces);

        /***********************************
         * Adding the attributes: normal
         ***********************************/

        //if(!ignoreNormals) {
        var padded_normals;
        if(pre_allocate){
            console.log("Allocating separate space for norms, colors.");
            padded_normals = copy_Float32Array_preallocated(verts, min_verts_capacity*3, "random");
            // if !ignoreNormals ...
            //is this correct? : padded_verts
            padded_normals.set(padded_verts);

            /*
            if (given_normals) {
                console.log("given_normals = ", given_normals);
                console.log(padded_normals.length);
                console.log(given_normals.length);
                padded_normals.set(given_normals);
            } else if (given_normals == null ) {
                ; //it's fine.
            }
            else {
                console.error(" no given_normals", given_normals);
            }
            */

            //var padded_colors = copy_Float32Array_preallocated(verts, min_verts_capacity*3, "random");
            //var uvs = copy_Float32Array_preallocated(new Float32Array([]), min_verts_capacity*3 * 0, "random");
            //padded_colors.set(padded_verts);
        }
        else{
            console.error("padded_normals = ?");
        }

        // We always need to allocate normals
        this.addAttribute( 'normal', new THREE.BufferAttribute( padded_normals, 3, true ) );

        //}
        /***********************************
         * Color
         ***********************************/

        var padded_color = copy_Float32Array_preallocated(verts, min_verts_capacity*3, "random");
        padded_color.set(padded_verts);
        this.addAttribute( 'color', new THREE.BufferAttribute( padded_color, 3, true ) );

        /***********************************
         * Other attributes
         ***********************************/

        //Other attributes:
        //    Specific to InstancedBufferGeometry
        //      values added to each coordinate:
        //       var offsets = new THREE.InstancedBufferAttribute( new Float32Array( instances * 3 ), 3, 1 );
        //         or     offsets = new THREE.InstancedBufferAttribute( DELTA_XYZ, 3,1 )??
        //     * this.addAttribute( 'offset', offsets );

        //     * this.addAttribute( 'orientation', new THREE.InstancedBufferAttribute( DELTA_XYZW, 4,1 ) );

        //this.addAttribute( 'uv', new THREE.BufferAttribute( new Float32Array( [..]), 2 );
        //geometry.faces

        //this.addAttribute( 'color', new THREE.BufferAttribute( padded_colors, 3, true ) ); //color is overidden
        //this.addAttribute( 'uv', new THREE.BufferAttribute( uvs, 2 ) );


        /***********************************
         * Setting the sizes (number of faces)
         ***********************************/

        // var nf3 = faces.length;
        var ntriangles = faces.length / POINTS_PER_FACE;

        this.__set_range_of_used_faces(ntriangles, POINTS_PER_FACE);

        //var materialIndex = 0;
        //this.addGroup( 0, faces.length*1-10, materialIndex );
        //scope.addGroup( groupStart, groupCount, materialIndex );

        /***********************************
         * Other settings of the attributes
         ***********************************/

        //this.index has a .count attribute  this.index.count
        this.index.dynamic=true;

        this.attributes.position.dynamic=true;
        this.attributes.normal.dynamic=true;
        // this.index.needsUpdate = true;

    };

    // ThreeJS does not use prototype-based OOP.
    /*
    This function may return new LiveBufferGeometry79 object, in this case the old geometry has to be replaced (by the new returned object) in the caller function.
    */

    /* deprecaed. Use ImpliSolid.update_geometry() instead.
     i.e. swap geometry & implicit_service:
     geom.update_geometry(IMPLISOLID, true)  ->  IMPLISOLID.update_geometry(geom, true)
     */
    this.update_geometry_ = function(implicit_service, ignoreDefaultNormals) {
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;
        const POINTS_PER_FACE = 3;

        const geometry = this;
        var nverts = implicit_service.get_v_size();
        var nfaces = implicit_service.get_f_size();

        if(nfaces > 0){
            var verts_address = implicit_service.get_v_ptr();
            var faces_address = implicit_service.get_f_ptr();

            var verts = Module.HEAPF32.subarray(
                verts_address/_FLOAT_SIZE,
                verts_address/_FLOAT_SIZE + 3*nverts);

            var faces = Module.HEAPU32.subarray(
                faces_address/_INT_SIZE,
                faces_address/_INT_SIZE + nfaces * POINTS_PER_FACE);
        }
        else{
            console.log("empty implicit");
            var verts = new Float32Array([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, 0,1,1 ]);
            var faces = new Uint32Array([0,1,2, 0,2,3, 0,4,5, 0,5,1, 1,5,6, 1,6,2, 2,6,3, 3,6,7, 4,5,6, 5,6,7]);
        }

        return this.update_geometry1(verts, faces, ignoreDefaultNormals, false);
    };

    this.__set_range_of_used_faces = function (ntriangles, POINTS_PER_FACE)
    {
        /*
            Sets which part of the faces array is actually used. This is necessary when we preallocate an extra capacity for updatable geometryies to avoid need for dynamically growing/changing the size at each update.
        */

        if(threejs_r71) {
            var gl_chunkSize=21845;
            // var ii = 0;this.offsets.push({start:ii, index: ii , count: Math.min(faces.length - ii, gl_chunkSize*3)});
            // if(threejs_r71) {...}
            this.offsets = [];
            var offsets = ntriangles / gl_chunkSize;
            for ( var i = 0; i < offsets; i ++ ) {

                var offset = {
                    start: i * gl_chunkSize * POINTS_PER_FACE,
                    index: i * gl_chunkSize * POINTS_PER_FACE,
                    count: Math.min( ntriangles - ( i * gl_chunkSize ), gl_chunkSize ) * POINTS_PER_FACE
                };

                this.offsets.push( offset );

            }
            this.used_faces_range = ntriangles * POINTS_PER_FACE;
            //this.computeBoundingSphere();
            //this.clearGroups();this.addGroup( 0, ntriangles * POINTS_PER_FACE, 0 );

        } else if(threejs_r79) {
            // assert(ntriangles);
            if (!ntriangles) {
                console.error("Warning: ntriangles==0 : " + ntriangles);
            }
            //vertices or faces? faces.
            this.setDrawRange( 0, ntriangles * POINTS_PER_FACE );
            this.used_faces_range = ntriangles * POINTS_PER_FACE;  // array size
        }
    };

    this.__set_needsUpdate_flag = function (ignoreNormals) {
        /*
            Sets the flags so that WebGL is notified of change/update in the geometry's faces or verts.
        */

        var geometry = this;
        //geometry.computeOffsets();
        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.position.dynamic = true;
        if(threejs_r71) {
            geometry.attributes.index.needsUpdate = true;
        } else if (threejs_r79) {
            geometry.index.needsUpdate = true;
            geometry.index.dynamic = true;
        }
        if(!ignoreNormals) {
            geometry.attributes.normal.needsUpdate = true;
            geometry.attributes.normal.dynamic = true;
        }
        //geometry.attributes.color.needsUpdate = true;
    };

    /**
        Known bug: The "Wireframe update" bug: When the BufferGrometry is shown using a wirefram material, it does not update the lines correctly. But the non-wirefram solid works perfectly.
        @param ignoreDefaultNormals=true; if false, uses the vertex's vector as the "default normal" (which is incorrect), i.e. the vector from O=(0,0,0) to each vector's position.
        So it is always suggested to be `true`. Note that we use a vertex-based normal here (i.e. a normal vector for each vertex, not for each face).
        @param use_wireframe=false should always be false. This is an attempt to fix the "wireframe update" bug.
    */
    this.update_geometry1 = function(verts, faces, ignoreDefaultNormals, use_wireframe) {
        use_wireframe = false;  // experimental. should be always false;

        if (faces.length == 0) {
            console.warn("Warning: faces.length == 0");
        }
        /*
            A wireframe BufferGeometry may have "2" vertex "indices" per line.
            A non-wireframe one uses "3" vertex "indices" per face.
        */
        if (use_wireframe)
            var POINTS_PER_FACE = 2;
        else
            var POINTS_PER_FACE = 3;

        // nverts actually not used
        // var nverts = verts.length / 3;  // todo: double check
        // var nfaces = faces.length / POINTS_PER_FACE;  // todo: double check
        // A BUG DETECTED HERE! IT SHOULD BE LIKE THIS INSTEAD:
        // var nfaces = faces.length / POINTS_PER_FACE;  // todo: double check

        const geometry = this;

        if(threejs_r71) {
            var old_nverts = geometry.attributes.position.array.length/3;  // Physical space size.
            var old_nfaces = geometry.attributes.index.array.length/POINTS_PER_FACE;
        } else if (threejs_r79) {
            var old_nverts = geometry.attributes.position.array.length/3;
            var old_nfaces = geometry.index.array.length/POINTS_PER_FACE;
        }

        // var nv3 = Math.min(nverts, old_nverts) * 3;
        // var nf3 = Math.min(nfaces, old_nfaces) * POINTS_PER_FACE;
        //assert(nv3 === verts.length);
        //assert(nf3 === faces.length);

        // *************************************
        // * The following code works only when we use a single instance of the geometry, i.e. used in one Mesh.
        // * It will not work well if the same BufferGeometry is used in more than one Mesh
        // * So the object will look fine if we disable the wireframe mesh.
        // *************************************
        var growth_needed = false;

        if(threejs_r71) {
            var availableVertsSize = geometry.attributes.position.array.length;
        } else if (threejs_r79) {
            var availableVertsSize = geometry.attributes.position.array.length;
        }

        if(verts.length > availableVertsSize){

            growth_needed = true;

        }else{
            // copy now
            geometry.attributes.position.array.set(verts);  // .set() can copy if the number of soure (verts) is smaller than target
            if(!ignoreDefaultNormals) {
                geometry.attributes.normal.array.set(verts);
            }
        }

        if(threejs_r71) {
            var availableFacesSize = geometry.attributes.index.array.length;
        } else if (threejs_r79) {
            var availableFacesSize = geometry.index.array.length;
        }
        if(faces.length > availableFacesSize){
            //geometry.attributes.index.array.set(faces.subarray(0, availableFacesSize));
            growth_needed = true;
        }else{
            // copy now

            if(threejs_r71) {
                geometry.attributes.index.array.set(faces);
            } else if (threejs_r79) {
                  // .set() can copy if the number of soure (verts) is smaller than target
                geometry.index.array.set(faces);
            }
        }

        if(growth_needed){
            // Why are we using faces_ here ??
            console.log("increasing capacity : availableFacesSize : " + availableFacesSize + " facesLength : " + faces.length);
            //this.dispose();
            var min_faces_capacity = Math.floor(Math.max(availableFacesSize/POINTS_PER_FACE, faces.length/POINTS_PER_FACE) * GROWTH_FACTOR + GROWTH_ADDITIONAL);
            var min_verts_capacity = Math.floor(Math.max(availableVertsSize/3, verts.length/3) * GROWTH_FACTOR + GROWTH_ADDITIONAL);

            var min_faces_capacity = faces.length/POINTS_PER_FACE;
            var min_verts_capacity = verts.length/3;

            console.error("capacity grow: "+min_faces_capacity+"  "+min_verts_capacity);

            //also copies
            this.allocate_buffers(verts, faces,  true, min_faces_capacity, min_verts_capacity);

            // It seems this is needed:
            // faces_ = faces;
            // verts_ = verts;

            //old solution was: create a new object
            //  var new_geometry= new LiveBufferGeometry79( verts, faces,  true, Math.max(availableFacesSize/POINTS_PER_FACE, faces.length/POINTS_PER_FACE) * 1.5 + 1, Math.max(availableVertsSize/3, verts.length/3) * 1.5 +1);

            geometry.__set_needsUpdate_flag(ignoreDefaultNormals);

            return false; //new_geometry;

        } else {

            assert(!growth_needed);
            var copied_faces = Math.min(faces.length, availableFacesSize);
            var ntriangles = copied_faces/POINTS_PER_FACE;
            this.__set_range_of_used_faces(ntriangles, POINTS_PER_FACE);
            assert(!growth_needed);

            //if(growth_needed){

            // Notify the changes

            geometry.__set_needsUpdate_flag(ignoreDefaultNormals);

            /*
            g =currentMeshes[0].geometry

            IMPLICIT.finish_geometry();
            var mc_properties_json = JSON.stringify({resolution: CONFIG.implisolid.default_mc_resolution, box: {xmin: -1, xmax: 1, ymin: -1 , ymax: 1, zmin: -1, zmax: 1}});
            IMPLICIT.build_geometry(......., mc_properties_json,);
            //g.update_geometry(IMPLICIT, false)
            IMPLICIT.update_geometry(g, false)


            for(var i=0;i<1000;i++){ IMPLICIT.finish_geometry();IMPLICIT.build_geometry(..........28, mc_properties_json, "sphere", i*0.1); g.update_geometry(IMPLICIT, false);}

            */
            return false;
        }
    };

    this.update_normals_from_array = function(normals) {
        this.computeVertexNormals();
        this.__set_needsUpdate_flag(false);
        return;

        //i.e. setNormals()
        // Only called after updateGeometry(), or after LiveGeometry()
        var geometry = this;
        // assert ensure length
        geometry.attributes.normal.array.set(normals);

        geometry.computeVertexNormals();
        geometry.attributes.normal.needsUpdate = true;
    };

    // Separation of concern: Should not use implicit_service. Only provide .update_normals_from_array()
    this.update_normals__ = function(implicit_service, mp5_str, x, ignore_root_matrix) {
        // Only called after updateGeometry(), or after LiveGeometry()

        //this.aaaaaaaaaA(x, mp5_str, ignore_root_matrix) {

        // var x = new Float32Array(nverts);
        console.error("geometry.update_normals() deprecated");

        // throw new Error("Deprecated method");

        var geom = this;
        //var gradients = implicit_service.make_normals_into_geometry(geom, mp5_str, x, ignore_root_matrix);
        implicit_service.make_normals_into_geometry(geom, mp5_str, x, ignore_root_matrix);

        /*

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;

        implicit_service.set_object(mp5_str, ignore_root_matrix);
        implicit_service.set_vect(x);  // overhead
        implicit_service.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
        var ptr = implicit_service.get_gradients_ptr();
        var ptr_len = implicit_service.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        //console.log("grad len = " +  ptr_len+ "  grad = " + gradients);  // x 4

        //for( var i = 0 ; i < ptr_len; i++) {
        //    gradients[i] += Math.random() * 0.2;
        //}


        geom.update_normals_from_array(gradients);


        implicit_service.unset_x();
        implicit_service.unset_object();
        */
    };

    this.use_default_normals_from_vertices = function () {
        // use default normals
        // candidate name: remove_normals()
        this.removeAttribute('normal');
        this.computeVertexNormals();
        this.__set_needsUpdate_flag(false);
    }


    // rename: get_min_z_from_mesh()
    this.get_minz = function(matrix, shape) {
        // shape is not needed really.
        return +1;
    }

    this.getVerticesArray = function (){
        return this.attributes.position.array;
    };
    this.getFacesIndexArray = function (){
        return this.index.array;
    };
    this.getNumVerts = function (){
        return this.attributes.position.count;
    };
    this.getNumFaces = function (){
        return this.index.count;
    };


    // this.computeBoundingBox = function...; // not needed.
    // this.computeBoundingSphere = function...; // not needed.

    this.allocate_buffers(verts_, faces_,  pre_allocate_, min_faces_capacity_, min_verts_capacity_);
};

LiveBufferGeometry79.prototype = Object.create( THREE.BufferGeometry.prototype );
LiveBufferGeometry79.prototype.constructor = LiveBufferGeometry79;
