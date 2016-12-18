// requires js_utils.js

//utils3d
var PS_UTILS_CLASS = function(IMPLICIT) {
    
    this.IMPLI1 = IMPLICIT.service2.service1;
    
    // todo: should work with point_count==0
    this.get_emc_array = function (verts_address, point_count) {
        // var verts_address = this.IMPLI1.get_pointset_ptr(name);
        // var point_count = this.IMPLI1.get_pointset_size(name);

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var float_array = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 1 * point_count);
        // console.error(float_array);

        if (point_count == 0) {
            console.warn("Error: zero points", point_count);
            //return // null; // this._make_empty_lines_mesh();
        }

        return float_array;
    }

    this.get_pointset_ = function (pointset_label) {
        assert(typeof pointset_label === "string");
        var n = this.IMPLI1.get_pointset_size(pointset_label);
        if (n == 0) {
            console.error("zero pointset ", pointset_label);
            return null;
        }
        var verts = this.get_emc_array(
            this.IMPLI1.get_pointset_ptr(pointset_label),
            n * 3 );
        return verts;
    }

    this.make_pointset = function (name, point_size, color_r, color_g, color_b) {

        //console.error("get_pointset: skipping");
        //return new THREE.Object3D();

        // r79 only

        // var point_size = 15.0; // millimeters
        var scale_up = 1.0;

        point_size = point_size !== undefined? point_size : 15.0;

        color_r = color_r !== undefined? color_r : 0.9;
        color_g = color_g !== undefined? color_g : 0.9;
        color_b = color_b !== undefined? color_b : 0.9;

        // var dx = Math.random() * 1.0 * scale_up * 1;
        // var dy = Math.random() * 1.0 * scale_up * 1;
        // var dz = Math.random() * 1.0 * scale_up * 1;

        // var dx = 1.0 * scale_up * 1;
        // var dy = 1.0 * scale_up * 1;
        // var dz = 1.0 * scale_up * 1;
        var dx = 0, dy = 0, dz = 0;

        /*
        var verts_address = this.IMPLI1.get_pointset_ptr(name);
        var point_count = this.IMPLI1.get_pointset_size(name);

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3 * point_count);
        */
        var point_count = this.IMPLI1.get_pointset_size(name);
        if (point_count == 0) {
            console.error("0 points from point-set ", name);
            return this._make_empty_lines_mesh();
        }

        var verts = this.get_emc_array(this.IMPLI1.get_pointset_ptr(name), 3 * point_count );
        // console.error(verts);


        var geometry = new THREE.BufferGeometry();
        var positions = verts; //new Float32Array( point_count * 3 );

        var colors = new Float32Array( point_count * 3 );
        for ( var i = 0; i < positions.length; i += 3 ) {
            positions[ i     ] = positions[ i     ] * scale_up + dx;
            positions[ i + 1 ] = positions[ i + 1 ] * scale_up + dy;
            positions[ i + 2 ] = positions[ i + 2 ] * scale_up + dz;
        }

        for ( var i = 0; i < positions.length; i += 3 ) {
            colors[ i ]     = color_r;
            colors[ i + 1 ] = color_g;
            colors[ i + 2 ] = color_b;
        }

        geometry.addAttribute( 'position', new THREE.BufferAttribute( positions, 3 ) );
        geometry.addAttribute( 'color', new THREE.BufferAttribute( colors, 3 ) );
        geometry.computeBoundingSphere();

        //geometry is ready

        var material = new THREE.PointsMaterial( { size: point_size, vertexColors: THREE.VertexColors } );
        var points = new THREE.Points( geometry, material );

        return points;
    }



    this.make_lines_mesh_v1 = function (verts12, indices12, color) {

        if (indices12.length == 0) {
            console.error("Warning: indices12.length =", indices12.length);
        }
        var re_allocate_buffers;
        //var geom = new LiveBufferGeometry79( verts12, indices12, re_allocate_buffers=true, 0, 0);
        var geom = new LiveBufferGeometry79( new Float32Array(new ArrayBuffer(0)), new Uint32Array(new ArrayBuffer(0)), re_allocate_buffers="Deliberately Empty", 0, 0);
        geom.update_geometry1(verts12, indices12, true, true);

        var wireframe_material = new THREE.MeshBasicMaterial( {
            color: color, wireframe: true, opacity:0.9,  transparent: true,} );
        var wire_msh = new THREE.Mesh(geom, wireframe_material);
        return wire_msh;
    }

    this._make_empty_lines_mesh = function () {
        return new THREE.LineSegments(new THREE.BoxGeometry(0.33333,0.3333333,0.33333333));
    }

    this.make_vectorfield_flow_lines_v1 = function (name1, name2, point_size) {
        //console.error("skipping   make_vectorfield_flow_lines_v1()");
        //return THREE.Object3D();

        // works!
        // vector field flow
        /*
        var color_r, color_g, color_b;
        color_r = color_r !== undefined? color_r : 0.9;
        color_g = color_g !== undefined? color_g : 0.9;
        color_b = color_b !== undefined? color_b : 0.9;
        */

        if (this.IMPLI1.get_pointset_ptr(name1) == 0) {

        }
        var verts1 = this.get_emc_array(this.IMPLI1.get_pointset_ptr(name1)+120*3*4*0, 3 * this.IMPLI1.get_pointset_size(name1) + 0  );
        var verts2 = this.get_emc_array(this.IMPLI1.get_pointset_ptr(name2)+120*3*4*0, 3 * this.IMPLI1.get_pointset_size(name2) + 0 );

        var n1 = verts1.length / 3;
        var n2 = verts2.length / 3;

        if (n1 == 0) {
            return this._make_empty_lines_mesh();
        }
        /*
        n1 = 4*10;
        n2 = 4*10;
        var debug_mode = true;
        verts1 = verts1.subarray(0, n1 * 3);
        verts2 = verts2.subarray(0, n2 * 3);
        */

        var verts12 = concatenate_array_float(verts1, verts2);
        /*
        var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * (n1 + n2) * 3);
        var verts12 = new Float32Array(raw_buffer);
        // var verts12 = new Float32Array(n1 + n2);
        verts12.subarray(0, n1 * 3).set(new Float32Array(verts1));
        verts12.subarray(n1 * 3, (n1+n2) * 3).set(new Float32Array(verts2));  // deliberate bug for TDD
        */

        /*
        for( var i=0; i<n2; i++) {
            var b = (i+n1)*3;
            verts12[b + 0] += (Math.random()*2-1)*10;
            verts12[b + 1] += (Math.random()*2-1)*10;
            verts12[b + 2] += (Math.random()*2-1)*10;
        }
        */


        if (n1 !== n2) {
            console.error("Need to have the same number of elements.");
        }

        /*
        // I don't use a geometry mesh with wiredframe. I use a LineSegments that uses a BufferGeometry

        var raw_buffer = new ArrayBuffer( Uint32Array.BYTES_PER_ELEMENT * (n1) * 2 );
        var indices12 = new Uint32Array(raw_buffer);
        for (var i=0; i < n1; i++) {
            indices12[i * 2] = i;
            indices12[i * 2 + 1] = i ; //+ n1;  //i+1; //i + n1;
            // indices12[i * 2] = 0;
            // indices12[i * 2 + 1] = i;
        }
        */


        var indices12 = make_indices_for_n1_n2(n1);

        /*
        var raw_buffer = new ArrayBuffer( Uint32Array.BYTES_PER_ELEMENT * (n1) * 3 );
        var indices12 = new Uint32Array(raw_buffer);
        for (var i=0; i < n1; i++) {
            indices12[i * 3] = i;
            indices12[i * 3 + 1] = i + n1; //i + n1;
            indices12[i * 3 + 2] = i; //i + n1;
        }
        my_assert(indices12.length == n1 * 3);
        */


        var wire_msh = this.make_lines_mesh_v1(verts12, indices12, 0xffffff);

        /*
        var line_material = new THREE.LineBasicMaterial({ vertexColors: THREE.VertexColors });
        var wire_msh = new THREE.LineSegments(geom, line_material);
        */
        return wire_msh;
    }


    this.make_vectorfield_flow_lines_v2 = function (name1, name2, point_size) {
        // vector field flow. Without using LiveGEometry (and using LineSegments)

        // does not work
        /*
        var color_r, color_g, color_b;
        color_r = color_r !== undefined? color_r : 0.9;
        color_g = color_g !== undefined? color_g : 0.9;
        color_b = color_b !== undefined? color_b : 0.9;
        */



        /*

                    var positions = [];
                    var next_positions_index = 0;
                    var colors = [];
                    var indices_array = [];
        */

        var verts1 = this.get_emc_array(this.IMPLI1.get_pointset_ptr(name1), 3 * this.IMPLI1.get_pointset_size(name1) );
        var verts2 = this.get_emc_array(this.IMPLI1.get_pointset_ptr(name2), 3 * this.IMPLI1.get_pointset_size(name2) );

        var n1 = verts1.length / 3;
        var n2 = verts2.length / 3;

        my_assert(n1 > 0);

        var verts12 = concatenate_array_float(verts1, verts2);
        /*
        var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * (n1 + n2) * 3);
        var verts12 = new Float32Array(raw_buffer);
        // var verts12 = new Float32Array(n1 + n2);
        verts12.subarray(0, n1 * 3).set(new Float32Array(verts1));
        verts12.subarray(n1 * 3, (n1+n2) * 3).set(new Float32Array(verts2));  // deliberate bug for TDD
        */


        if (n1 !== n2) {
            console.error("Need to have the same number of elements.");
        }

        // I don't use a geometry mesh with wiredframe. I use a LineSegments that uses a BufferGeometry

        var raw_buffer = new ArrayBuffer( Uint16Array.BYTES_PER_ELEMENT * (n1) * 2 );
        var indices12 = new Uint16Array(raw_buffer);
        for (var i=0; i < n1; i++) {
            indices12[i * 2] = i;
            indices12[i * 2 + 1] = i + n1;  //i+1; //i + n1;
        }


        var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * (n1) * 3 );
        var colors_ta = new Float32Array(raw_buffer);
        for (var i=0; i < n1; i++) {
            colors_ta[i * 3] = (i/10.0) % 1.0;
            colors_ta[i + 3 + 1] = (i/10.0) % 1.0;
            colors_ta[i + 3 + 2] = (i/10.0) % 1.0;
        }

        // var indices12 = new Uint16Array( indices_array );
        // var verts12 = new Float32Array( positions );
        // var colors_ta = new Float32Array( colors );

        var geom = new THREE.BufferGeometry();
        geom.setIndex( new THREE.BufferAttribute( indices12, 1 ) );
        geom.addAttribute( 'position', new THREE.BufferAttribute( verts12, 3 ) );
        geom.addAttribute( 'color', new THREE.BufferAttribute( colors_ta, 3 ) );
        geom.computeBoundingSphere();

        // var line_material = new THREE.LineBasicMaterial({ vertexColors: THREE.VertexColors });
        var line_material = new THREE.MeshBasicMaterial( {
            color: 0xffffff, wireframe: true, opacity:0.8,  transparent: true,} );
        var wire_msh = new THREE.LineSegments(geom, line_material);

        return wire_msh;
    }


    this.hairy_object_based_on_vectors = function (vects, gradients, alpha1, alpha2, color) {

        var n = gradients.length / 3;
        var n1 = vects.length / 3;
        my_assert(n == n1);

        var arrow_tips = new Float32Array(new ArrayBuffer( n * 3 * Float32Array.BYTES_PER_ELEMENT));
        var arrow_bases = new Float32Array(new ArrayBuffer( n * 3 * Float32Array.BYTES_PER_ELEMENT));
        for (var i = 0; i < n; ++i) {
            for (var d = 0; d < 3; d++) {
                var j = (i * 3) + d;
                arrow_tips[j] = vects[j] + alpha2 * gradients[j];
                arrow_bases[j] = vects[j] + alpha1 * gradients[j];
            }
        }

        var verts12 = concatenate_array_float(arrow_bases, arrow_tips);
        /*
        // repeated code
        var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * (n + n) * 3);
        var verts12 = new Float32Array(raw_buffer);
        verts12.subarray(0, n * 3).set(new Float32Array(arrow_bases));
        verts12.subarray(n * 3, (n+n) * 3).set(new Float32Array(arrow_tips));  //
        */

        var indices12 = make_indices_for_n1_n2(n);
        /*
        //repeated code
        var raw_buffer = new ArrayBuffer( Uint32Array.BYTES_PER_ELEMENT * (n) * 3 );
        var indices12 = new Uint32Array(raw_buffer);
        for (var i=0; i < n; i++) {
            indices12[i * 3] = i;
            indices12[i * 3 + 1] = i + n; //i + n;
            indices12[i * 3 + 2] = i + n; //i + n;
        }
        my_assert(indices12.length == n * 3);
        */

        var wire_msh = this.make_lines_mesh_v1(verts12, indices12, color)
        return wire_msh;
    }



    this.visualise_normals = function (shape_json, pointset_label, alpha1, alpha2, color) {
        /*
            Draws lines from alpha1*grad to alpha2*grad.
            alpha1,2 are the start and end of the lines parallel to the gradients. Most likely negative alpha values are used.
            alpha=0 means points specified from the pointset specified by label pointset_label. Using color.
            Example: alpha1=0.0 and alpha2=-1.0 will draw normals (outwards). color = 0x88ff00;
        */

        /*
        // todo: refactor
        var n = this.IMPLI1.get_pointset_size(pointset_label);
        if (n == 0) {
            return this._make_empty_lines_mesh();
        }
        var verts = this.get_emc_array(this.IMPLI1.get_pointset_ptr(pointset_label),
            n * 3 );
        */
        var ps = this.get_pointset_ (pointset_label);

        if (ps === null) {
            console.error("zero points from pointset", pointset_label);
            return this._make_empty_lines_mesh();
        }
        var verts = ps;
        var n = verts.length/3;
        assert( n > 0);


        // verts = arrow_bases

        // see update_arrows()

        var ignore_root_matrix = false;
        var mp5_str = shape_json;
        // console.error(mp5_str);
    var result1 =
        this.IMPLI1.set_object(mp5_str, ignore_root_matrix);
    console.log("SET_OBJECT() RESULTS:", result1);
        this.IMPLI1.set_vect(verts);
        this.IMPLI1.calculate_implicit_gradients();

        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        var ptr = this.IMPLI1.get_gradients_ptr();
        var ptr_len = this.IMPLI1.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(
            ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        my_assert(ptr_len == n*3);

        for (var i = 0; i < n; ++i) {
            var r = Math.sqrt(gradients[i*3 + 0] * gradients[i*3 + 0] + gradients[i*3 + 1] * gradients[i*3 + 1] + gradients[i*3 + 2] * gradients[i*3 + 2]);
            gradients[i*3 + 0] /= r;
            gradients[i*3 + 1] /= r;
            gradients[i*3 + 2] /= r;
        }


        var wire_msh = this.hairy_object_based_on_vectors(verts, gradients, alpha1, alpha2, color);

        /*
        var raw_buffer = new ArrayBuffer( n * 3 * Float32Array.BYTES_PER_ELEMENT);
        var arrow_tips = new Float32Array(raw_buffer);
        var raw_buffer = new ArrayBuffer( n * 3 * Float32Array.BYTES_PER_ELEMENT);
        var arrow_bases = new Float32Array(raw_buffer);
        for (var i = 0; i < n; ++i) {
            for (var d = 0; d < 3; d++) {
                var j = (i * 3) + d;
                arrow_tips[j] = verts[j] + alpha2 * gradients[j];
                arrow_bases[j] = verts[j] + alpha1 * gradients[j];
            }
        }

        // repeated code
        var raw_buffer = new ArrayBuffer( Float32Array.BYTES_PER_ELEMENT * (n + n) * 3);
        var verts12 = new Float32Array(raw_buffer);
        verts12.subarray(0, n * 3).set(new Float32Array(arrow_bases));
        verts12.subarray(n * 3, (n+n) * 3).set(new Float32Array(arrow_tips));  //

        //repeated code
        var raw_buffer = new ArrayBuffer( Uint32Array.BYTES_PER_ELEMENT * (n) * 3 );
        var indices12 = new Uint32Array(raw_buffer);
        for (var i=0; i < n; i++) {
            indices12[i * 3] = i;
            indices12[i * 3 + 1] = i + n; //i + n;
            indices12[i * 3 + 2] = i + n; //i + n;
        }
        my_assert(indices12.length == n * 3);

        var wire_msh = this.make_lines_mesh_v1(verts12, indices12, color)
        */

        this.IMPLI1.unset_x();
    var result1 =
        this.IMPLI1.unset_object(1);
    console.log("UNSET_OBJ() RESULTS:", result1);
        return wire_msh;
    }
};
