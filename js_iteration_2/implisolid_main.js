'use strict';


function init(service) {
    'use strict';
    //main = Module.cwrap('main', 'number', []);
    //var service={}; //= newProducer //is an interface
    service.build_geometry = Module.cwrap('build_geometry', null, ['number', 'number']);
    service.get_v_size = Module.cwrap('get_v_size', 'number', []);
    service.get_f_size = Module.cwrap('get_f_size', 'number', []);
    service.get_v = Module.cwrap('get_v', null, ['number']);
    service.get_f = Module.cwrap('get_f', null, ['number']);
    service.get_v_ptr = Module.cwrap('get_v_ptr', 'number', []);
    service.get_f_ptr = Module.cwrap('get_f_ptr', 'number', []);
    service.finish_geometry = Module.cwrap('finish_geometry', null, []);

    service.init = function(){ service.needsFinish = false; }
    service.finish_with = function (){
        //after the last round.
        if(!this.needsFinish){
            console.error("cannot `finish_geometry()`. Geometry not produced.");
        }
        service.finish_geometry();
        service.needsFinish = false;
    }

    service.init();
    return service;
}


var ImplicitService = function(){
    init(this);
    this.make_geometry = function (params) {
        var startTime = new Date();
        const _FLOAT_SIZE = Float32Array.BYTES_PER_ELEMENT;
        const _INT_SIZE = Uint32Array.BYTES_PER_ELEMENT

        if(this.needsFinish) {
            this.finish_geometry();
            this.needsFinish = false;
        }
        this.build_geometry( 28, params["subjective_time"]);
        this.needsFinish = true;

        var nverts = this.get_v_size();
        var nfaces = this.get_f_size();

        var verts_address = this.get_v_ptr();
        var faces_address = this.get_f_ptr();

        var verts = Module.HEAPF32.subarray(verts_address/_FLOAT_SIZE, verts_address/_FLOAT_SIZE + 3*nverts);
        var faces = Module.HEAPU32.subarray(faces_address/_INT_SIZE, faces_address/_INT_SIZE + 3*nfaces);

        var allocate_buffers = true;
        var geom = new LiveBufferGeometry73(verts, faces, allocate_buffers);

        var endTime = new Date();
        var timeDiff = endTime - startTime;

        //report_time(timeDiff, function(){hist();});

        return geom;
    };
    this.getLiveGeometry = function(){
        var geom = this.make_geometry( {subjective_time: 0.0} );
        return geom;
    }

};

var IMPLICIT = null;
function _on_cpp_loaded() {
    console.log("C++ ready.");
    //IMPLISOLID.
    IMPLICIT = new ImplicitService();
};


/* Put hte following in the HTML
<script>
        Module={preRun:[],
        onRuntimeInitialized: _on_cpp_loaded,
    };
</script>
<script type="text/javascript" src="mcc2.cpp.js"></script>
*/
