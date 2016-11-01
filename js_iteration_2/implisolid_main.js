'use strict';


function init(service) {
    'use strict';
    //main = Module.cwrap('main', 'number', []);
    //var service={}; //= newProducer //is an interface
    service.build_geometry = Module.cwrap('build_geometry', null, [ 'string', 'string']);
    service.get_v_size = Module.cwrap('get_v_size', 'number', []);
    service.get_f_size = Module.cwrap('get_f_size', 'number', []);
    service.get_v = Module.cwrap('get_v', null, ['number']);
    service.get_f = Module.cwrap('get_f', null, ['number']);
    service.get_v_ptr = Module.cwrap('get_v_ptr', 'number', []);
    service.get_f_ptr = Module.cwrap('get_f_ptr', 'number', []);
    service.finish_geometry = Module.cwrap('finish_geometry', null, []);

    service.set_object = Module.cwrap('set_object', 'number', ['string', 'number']);
    service.unset_object = Module.cwrap('unset_object', 'number', ['number']);
    service.set_x = Module.cwrap('set_x', 'number', ['number', 'number']);  // sends the x points for evaluation of implicit or gradient
    service.unset_x = Module.cwrap('unset_x', null, []);
    service.calculate_implicit_values = Module.cwrap('calculate_implicit_values', null, []);
    service.get_values_ptr = Module.cwrap('get_values_ptr', 'number', []);
    service.get_values_size = Module.cwrap('get_values_size', 'number', []);
    service.calculate_implicit_gradients = Module.cwrap('calculate_implicit_gradients', null, ['number']);  // boolean argument to normalize and reverse the vevtors, suitable for rendering.
    service.get_gradients_ptr = Module.cwrap('get_gradients_ptr', 'number', []);
    service.get_gradients_size = Module.cwrap('get_gradients_size', 'number', []);

    service.get_pointset_ptr = Module.cwrap('get_pointset_ptr', 'number', ['string']);
    service.get_pointset_size = Module.cwrap('get_pointset_size', 'number', ['string']);

    if (Module["_about"]) {
        service.about = Module.cwrap('about', null, []);
    } else {
        service.about = "C++ method problem: no about()";
        console.error("C++ method problem: a correct version not found");
    }
    /*
    try {
        service.about = Module.cwrap('about', null, []);
    } catch(err) {
        service.about = "C++ method problem";
        console.log("C++ method problem: not found");
    }
    */

    service.init = function(){
        service.needs_deallocation = false;

        this.use_II = true;
        this.use_II_qem = true;
    }

    service.finish_with = function (){
        //after the last round.
        if(!this.needs_deallocation){
            console.error("cannot `finish_geometry()`. Geometry not produced.");
        }
        service.finish_geometry();
        service.needs_deallocation = false;
    }

    service.set_vect = function (float32Array) {
        // Accesses module
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        if (float32Array.length % 3 != 0) {console.error("bad input array");};
        var nverts = float32Array.length / 3;
        var verts_space = Module._malloc(_FLOAT_SIZE*3*nverts);
        Module.HEAPF32.subarray(verts_space/_FLOAT_SIZE, verts_space/_FLOAT_SIZE + 3*nverts).set(float32Array);
        var result = this.set_x(verts_space, nverts);
        console.log("result: "+result);
        if (!result) {
            console.error("Something went wrong: ", result);
        }
        Module._free( verts_space );

    };

    service.update_geometry = function(geometry, ignoreNormals) {

        var implicit_service = this;
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT;
        const POINTS_PER_FACE = 3;

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

        return geometry.update_geometry1(verts, faces, ignoreNormals, false);
    };


    //was: geom. update_normals (service ..)
    service.make_normals_into_geometry = function(geom, mp5_str, x, ignore_root_matrix) {
        geom.removeAttribute('normal');
        geom.computeVertexNormals();
        geom.__set_needsUpdate_flag(false);
        return;


        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        // console.error(x);
        /*
        for( var i = 0 ; i < x.length; i++) {
            x[i] += Math.random() * 0.9;
        }
        */

        service.set_object(mp5_str, ignore_root_matrix);
        service.set_vect(x);  // overhead
        service.calculate_implicit_gradients(true);  // Why TRUE doe snot have any effect?
        var ptr = service.get_gradients_ptr();
        var ptr_len = service.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        //console.log("grad len = " +  ptr_len+ "  grad = " + gradients);  // x 4

        //for( var i = 0 ; i < ptr_len; i++) {
        //    gradients[i] += Math.random() * 0.2;
        //}

        geom.update_normals_from_array(gradients);

        service.unset_x();
        service.unset_object();

    };


    service.init();

    return service;
}




/*
Function: query_implicit_values(shape_str, points, reduce_callback)
Evaluates the points **points** in the implicit function. The function is specified as a mp5-fragment (json) string **shape_str**.
The values are not directly returned. Instead, a callback is used to prepare a value.
The reason a callback is used is that the resources need to be freed.

Returned value:
The value returned us the value returned by the callback reduce_callback.
In examples below, the return value is:
- Example1: returns Float32Array of implicit values in the following usage.
- Example2: returns   true for collision, false for no collision,

Example 1:
    var my_shape = mainModel.root.sons[#nr]
    var verts_array = new Float32Array([x1,y1,z1, x2,y2,z2, ..., xn,yn,zn])
    var implicit_vals = query_implicit_values(my_shape, points,  // false
        reduce_callback = function (values_tarray) {
            // here we create a view on the HEAP from position
            // ( ptr / _FLOAT_SIZE ) to position (ptr / _FLOAT_SIZE + ptr_len)
            var result = new Float32Array(values_tarray);  // clones (makes a copy of) data
            return result;
        });
Example 2:
    var my_shape = mainModel.root.sons[#nr]
    var verts_array = new Float32Array([x1,y1,z1, x2,y2,z2, ..., xn,yn,zn])
    var has_positive_vals = query_implicit_values(my_shape, verts_array,  // true
        function (values_tarray) {
            var found_positive_val = false;
            for (var i=0; i< values_tarray.length; i++)
            {
                if (values_tarray[i] >= - CONFIG.collision_detection.epsilon)
                {
                    found_positive_val = true;
                    break;
                }
            }
            var result =  found_positive_val;
            return result;
        });


Implementation Details:

1. We do not have direct access to C++ objects from Javascript, so we call the next functions
    set_object:     sets the object specified by mp5_str as the current object for any calculation to follow until it gets freed by unset_object
    set_x:          does something similar but for the input x .

2. There is a repeating pattern taking place:
    JS asks for C++ pointers, and the pointers' data is accessed using HEAPF32(start, end)

3. Don't forget to Module._free() andter Module._alloc().

*/

function _query_implicit_values(mp5_str, points, reduce_callback)
{
    assert(points instanceof Float32Array);                   /* simple check */
    assert(typeof reduce_callback === 'function');             /* the callback that produces the return value of the function */

    /* Declarations */

    var _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;                    /* We ll need the _FLOAT_SIZE in bytes, when we deal with allocations and HEAPF32*/
    var nverts = points.length / 3;                                       /* the number of vertices for which we want to calculate the implicit value */
    var verts_space_address = Module._malloc(_FLOAT_SIZE * 3 * nverts);  /* This allocates space in the C++ side, in order to create the array of vertices. */

    var ignore_root_matrix = false;

    /* body */
    Module.HEAPF32.subarray(verts_space_address / _FLOAT_SIZE, verts_space_address / _FLOAT_SIZE + 3 * nverts
        ).set(points);

    var obj_id = IMPLICIT.set_object(mp5_str, ignore_root_matrix);                    /* Allocations on the C++ side */
    var success = IMPLICIT.set_x(verts_space_address, nverts);
    //todo: rename "success"
    if (!success){
        console.log("set_x returned false . Probably an error in memory allocation. nverts was: " + nverts);
        IMPLICIT.unset_object(obj_id);
        return false;
    }
    // IMPLICIT.set_x_with_matrix(verts_space_address, nverts, matrix);     /* can create a function like this so we dont have the matrix issue */
    IMPLICIT.calculate_implicit_values();                                /* The actual implicit value calculation*/

    var ptr = IMPLICIT.get_values_ptr();                                 /* retrieve a pointer to the position in memory of the calculated array */
    var ptr_len = IMPLICIT.get_values_size()
    var values_tarray = Module.HEAPF32.subarray(ptr / _FLOAT_SIZE , ptr / _FLOAT_SIZE + ptr_len );

    /*
    var found_positive_val = false;
    if (typeof return_boolean === 'boolean' && return_boolean === true)
    {

       for (var i=0; i< values_tarray.length; i++)
       {
           if (values_tarray[i] >= - CONFIG.collision_detection.epsilon)
           {
               found_positive_val = true;
               break;
           }
       }
       result =  found_positive_val;
    } else if (typeof return_boolean === 'boolean' && return_boolean === false)
    {
        var result = new Float32Array(values_tarray);    // here we create a view on the HEAP from position
                                             // ( ptr / _FLOAT_SIZE ) to position (ptr / _FLOAT_SIZE + ptr_len)
    } else if (typeof return_boolean === 'function') {
        var reduce_callback = return_boolean;
        var result = reduce_callback(values_tarray);
    }
    */
    var result = reduce_callback(values_tarray);

    //Bug! Forgot to FREE!!!
    Module._free( verts_space_address );
    IMPLICIT.unset_object(obj_id);    /* Free allocated C++ memory */
    IMPLICIT.unset_x();
    return  result;
}

var ImplicitService = function(){
    init(this);

    this.make_geometry = function (shape_params, mc_params) {
        var startTime = new Date();
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT

        if(this.needs_deallocation) {
            this.finish_geometry();
            this.needs_deallocation = false;
        }

        //console.log("mc_params.resolution " + mc_params.resolution);
        //mc_params.resolution = 40;

        var mp5_str = JSON.stringify(shape_params);
        this.build_geometry(mp5_str, JSON.stringify(mc_params));
        this.needs_deallocation = true;

        var nverts = this.get_v_size();
        var nfaces = this.get_f_size();

        var verts_address = this.get_v_ptr();
        var faces_address = this.get_f_ptr();

        var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);

        var allocate_buffers = true;
        var geom = new LiveBufferGeometry79(verts, faces, allocate_buffers);

        // Set the normals
        var ignore_root_matrix = mc_params.ignore_root_matrix;  // Does not need other (MC-related) arguments.
        //geom.update_normals(this, verts, mp5_str, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

        this.make_normals_into_geometry(geom, mp5_str, verts, ignore_root_matrix);  // Evaluates the implicit function and sets the goemetry's normals based on it.

        //this.aaaaaaaaa(verts);


        var endTime = new Date();
        var timeDiff = endTime - startTime;

        //report_time(timeDiff, function(){hist();});

        return geom;
    };
    /*
    this.aaaaaaaaa(x, mp5_str, ignore_root_matrix) {

        // var x = new Float32Array(nverts);

        this.set_object(mp5_str, ignore_root_matrix);
        this.set_vect(x);  // overhead
        this.calculate_implicit_gradients();
        var ptr = this.get_gradients_ptr();
        var ptr_len = this.get_gradients_size();
        var gradients = Module.HEAPF32.subarray(ptr/_FLOAT_SIZE, ptr/_FLOAT_SIZE + ptr_len);
        //console.log("grad len = " +  ptr_len+ "  grad = " + gradients);  // x 4

        geom.update_normals(gradients);

        this.unset_x();
        this.unset_object();

    }
    */
    //This method is called by the designer to obtain the geometry from the ImplicitService
    this.getLiveGeometry = function(dict, bbox, ignore_root_matrix) {
        //var mc_properties = {resolution: 28, box: {xmin: -1, xmax: 1, ymin: -1, ymax: 1, zmin: -1, zmax: 1}};

        //var shape_properties = {type:"sphere",displayColor:{x:0.38015037447759337,y:0.6015094592616681,z:0.9774198226067741},matrix:[10,0,0,92.9405888205127,0,10,0,101.93969389296757,0,0,10,8.59828143220919,0,0,0,1],index:7935813}
        //var shape_properties = {type:"sphere",displayColor:{x:0.38015037447759337,y:0.6015094592616681,z:0.9774198226067741},matrix:[10,0,0,92.9405888205127,0,10,0,101.93969389296757,0,0,10,8.59828143220919,0,0,0,1],index:7935813}
        /*{subjective_time: 0.0, implicit_obj_name: "sphere"*/

        //var shape_properties = {type:"meta_balls",time: 0.0};
        var shape_properties = dict;
        //var radius = 3.0;
        //var mc_properties = {resolution: 28, box: {xmin: -radius+1, xmax: radius, ymin: -radius, ymax: radius, zmin: -radius, zmax: radius}};
        //var shape_properties = {type:"simple_sphere", radius: radius};

        // var shape_properties = {type:"egg",displayColor:{x:0.38015037447759337,y:0.6015094592616681,z:0.9774198226067741},matrix:[10,0,0,92.9405888205127,0,10,0,101.93969389296757,0,0,10,8.59828143220919,0,0,0,1],index:7935813}
        // var s = 10;
        // var mc_properties = {resolution: 28, box: {xmin: 100-s, xmax: 100+s, ymin: 100-s, ymax: 100+s, zmin: 5-s, zmax: 5+s}};

        //var shape_properties = {type:"egg",displayColor:{x:0.38015037447759337,y:0.6015094592616681,z:0.9774198226067741},matrix:[1,0,0,92.9405888205127-100,0,1,0,101.93969389296757-100,0,0,1,8.59828143220919-5,0,0,0,1],index:7935813}
        var s = 1;
        //var mc_properties = {resolution: 28, box: {xmin: -s, xmax: s, ymin: -s, ymax: s, zmin: -s, zmax: s}};
        //var mc_properties = {resolution: 28, box: {xmin: 92.9405888205127-100-s, xmax: 92.9405888205127-100+s, ymin: 101.93969389296757-100-s, ymax: 101.93969389296757-100+s, zmin: 8.59828143220919-5-s, zmax: 8.59828143220919-5+s}};

        //implicit_double_mushroom center will be zero.

        /*
        //shape_properties.type = "egg";
        var m = shape_properties.matrix;
        var bb ={};
        var dd = 0.0;
        var wx = Math.sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        var wy = Math.sqrt(m[4]*m[4] + m[5]*m[5] + m[6]*m[6]);
        var wz = Math.sqrt(m[8]*m[8] + m[9]*m[9] + m[10]*m[10]);
        bb["xmin"] = m[3] - wx/2 +dd;
        bb["xmax"] = m[3] + wx/2-dd;

        bb["ymin"] = m[7] - wy/2 +dd;
        bb["ymax"] = m[7] + wy/2-dd;

        bb["zmin"] = m[11] - wz/2 +dd;
        bb["zmax"] = m[11] + wz/2-dd;
        */
        assert(bbox, "yOU need to specify the bounding box");


        var bb ={};
        var sc = 1.0;
        var test = 0.;
        bb["xmin"] = bbox.min.x * sc + test;
        bb["xmax"] = bbox.max.x * sc - test;

        bb["ymin"] = bbox.min.y * sc + test;
        bb["ymax"] = bbox.max.y * sc - test;

        bb["zmin"] = bbox.min.z * sc + test;
        bb["zmax"] = bbox.max.z * sc - test;

        _expect(bb["xmin"], "boundingbox has null");
        _expect(bb["xmax"], "boundingbox has null");
        _expect(bb["ymin"], "boundingbox has null");
        _expect(bb["ymax"], "boundingbox has null");
        _expect(bb["zmin"], "boundingbox has null");
        _expect(bb["zmax"], "boundingbox has null");


        var mc_res = CONFIG.implisolid.default_mc_resolution;
        var mc_properties = {
            resolution: getResolution(bb),
            box: bb,
            ignore_root_matrix: ignore_root_matrix,

            vresampl: {iters: this.use_II?1:0, c: 1.0},
            projection: {enabled: this.use_II?1:0},
            qem: {enabled: (this.use_II && this.use_II_qem)?1:0},
            subdiv: {enabled: 0},
            overall_repeats : 1,
            debug: {
                enabled_pointsets: 0,
                post_subdiv_noise: 0.0
            }
        };


        console.log (" mc properties : " + JSON.stringify(mc_properties));
        var geom = this.make_geometry(shape_properties, mc_properties);
        return geom;
    };

    this.query_implicit_values = _query_implicit_values;

};

var IMPLICIT = null;
function _on_cpp_loaded() {
    console.log("C++ ready.");
    //IMPLISOLID.
    IMPLICIT = new ImplicitService();

    assert = _assert_000;
};

function getResolution(bb){
    return CONFIG.implisolid.default_mc_resolution;
    const max_value = 40;
    const min_value = 14;
    const factor = CONFIG.implisolid.default_mc_resolution;
    var max_length = Math.max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"], bb["zmax"] - bb["zmin"]);
    var tmp =  Math.min(max_value,max_length*factor);
    return  Math.floor(Math.max(tmp,min_value));

    // 1 -> 28
    // 2 -> 48


}

/* Put the following in the HTML
<script>
        Module={preRun:[],
        onRuntimeInitialized: _on_cpp_loaded,
    };
</script>
<script type="text/javascript" src="mcc2.cpp.js"></script>
*/


//function test_update1(t, mesh){
function test_update1(t, mesh, dict){
    var g = mesh.geometry;

    IMPLICIT.finish_geometry();
    IMPLICIT.needs_deallocation = false;

    var s = Math.sin(t)*3+3;
    console.log("s="+s);

    var mc_properties = {resolution: 28, box: {xmin: (-1-s)*0, xmax: 1+s, ymin: (-1-s)*0 , ymax: 1+s, zmin: (-1-s)*0, zmax: 1+s}};

    // Aliasing test:
    //var res = Math.floor(28*(s+1)/3)-4; var sz = (res) * 0.14/2;
    //var mc_properties = {resolution: res, box: {xmin: 0, xmax: sz, ymin: 0 , ymax: sz, zmin: 0, zmax: sz}};

    //var mc_properties = {resolution: 28, box: {xmin: -1, xmax: 1, ymin: -1, ymax: 1, zmin: -1, zmax: 1}};
    //var new_geometry = IMPLICIT.build_geometry(28, mc_properties, "sphere", 0);

    //var shape_properties = {type: "sphere",matrix:[10,0,0,92.9405888205127,0,10,0,101.93969389296757,0,0,10,8.59828143220919,0,0,0,1]};

    //var shape_properties = mesh.parentShape.getDict1();
    console.log()
    if(!dict){
        //var shape_properties = {type:"meta_balls",time: t };
        var shape_properties = {type:"simple_sphere", radius: 3.0};
    }
    else
    {
        var shape_properties = dict;
    }

    var new_geometry = IMPLICIT.build_geometry(
        JSON.stringify(shape_properties) ,
        JSON.stringify(mc_properties));


    IMPLICIT.needs_deallocation = true;

    if(new_geometry){
        mesh.geometry = new_geometry;
        g = new_geometry;
    }
    // g.update_geometry(IMPLICIT, true);
    IMPLICIT.update_geometry(g, true);
    //??????? g.update_normals(IMPLICIT);
    var x = "?";
    IMPLICIT.make_normals_into_geometry(g, JSON.stringify(shape_properties), x, false);

}


function test_update2(t){
    var g = currentMeshes[0].geometry;

    // var new_geometry = g.update_geometry(IMPLICIT, false);
    var new_geometry = IMPLICIT.update_geometry(g, false);
    // g.update_normals(IMPLICIT);
    var x = "?";
    var mp5_string = "?";
    IMPLICIT.make_normals_into_geometry(g, mp5_String, x, false);

    if(new_geometry){
        currentMeshes[0].geometry = new_geometry;
    }

}

/*
var t=0;m=currentMeshes[0];test_update1(t, m);var iid=setInterval(function(){test_update1(t+=0.02, m);},6+9);

var t=0;m=currentMeshes[0]; d=m.parentShape.getDict1();test_update1(t, m);var iid=setInterval(function(){test_update1(t+=0.02, m,d);}, 66);
*/
