/*
    this.getBoundingBox_ = function(ignore_root_matrix) {

        return getBoundingBoxForSingleShapeMatrix(this.matrix, ignore_root_matrix, 0.5, 0.5, 0.5);

    }
    return getBoundingBoxForTree(this, ignore_root_matrix);

*/

const getBoundingBoxForSingleShape = function (matrix, ignore_root_matrix) {
        return getBoundingBoxForSingleShapeMatrix(matrix, ignore_root_matrix, 0.5, 0.5, 0.5);
}

function dictIsSingleShape(dict) {
    var type = dict.type;
    if (type === "Union" || type === "Difference" || type === "Intersection" || type === "root") {
        return false;
    }
    return true;
}
function dictIsTree(dict) {
    return ! dictIsSingleShape(dict);
}

function xxxxxxx(type, children, root_matrix, ignore_root_matrix) {
       getBoundingBoxForTree = getBoundingBoxForDict;
       var result_bbox = { min: {x: Infinity, y: Infinity, z: Infinity}, max: {x: -Infinity, y: -Infinity, z: -Infinity}};
        var n = children.length;

        if(type == "Difference") {
            n = 1;
        }

        for(var i = 0 ; i < n; i++) {
            var sbbox = getBoundingBoxForTree(children[i], false);  // never ignore the matrix if not root
            //updateBoundingBoxForSubBoundingBox(result_bbox, children[i].matrix ,  sbbox, true);
            updateBoundingBoxForSubBoundingBox(result_bbox, null,  sbbox, true);  // just update, dont apply the matrix
            //console.log(JSON.stringify(result_bbox));
        }
        // rotate the final result
        //return result_bbox;
        //console.log("ignor toot m : " + ignore_root_matrix);
        //console.log(JSON.stringify(result_bbox));
        //console.log(root_matrix);
        var ret = updateBoundingBoxForSubBoundingBox(null, root_matrix ,  result_bbox, ignore_root_matrix);
        //console.log(JSON.stringify(ret));
        return ret;
}

const getBoundingBoxForDict = function (dict, ignore_root_matrix) {
    if(dictIsSingleShape(dict)) {
        //return shape.getBoundingBox_(ignore_root_matrix);
        return getBoundingBoxForSingleShape(dict.matrix, ignore_root_matrix);
    } else if (dictIsTree(dict)) {
        return xxxxxxx(dict.type, dict.children, dict.matrix, ignore_root_matrix);
    } else {
        console.error("not a shape neither a tree");
    }
};


var getBoundingBoxForTree = function (shape, ignore_root_matrix) {
    console.error("not here");

    if(isSingleShape(shape)) {

        return shape.getBoundingBox_(ignore_root_matrix);

    } else if (isTree(shape)) {

        return xxxxxxx(shape.type, shape.sons, shape.matrix, ignore_root_matrix);

    } else {
        console.error("not a shape neither a tree");
    }
};

// sizeOfMultiple
getBoundingBoxForGroup = function (shape_list) {
    getBoundingBoxForTree = getBoundingBoxForDict;

    //assert(shape_list[0]);
    assert(!!shape_list.push, "shape_list has to be of type list ", shape_list);
    var result_bbox = { min: {x: Infinity, y: Infinity, z: Infinity}, max: {x: -Infinity, y: -Infinity, z: -Infinity}};
    var n = shape_list.length;

    for(var i = 0 ; i < n; i++) {
        var sbbox = getBoundingBoxForTree(shape_list[i], false);  // never ignore the matrix if not root
        updateBoundingBoxForSubBoundingBox(result_bbox, null,  sbbox, true);  // just update, dont apply the matrix
    }
    var matrix = null;
    //if(matrix) {
    //    var ret = updateBoundingBoxForSubBoundingBox(
    //        null, shape.matrix ,  result_bbox, ignore_root_matrix );
    //    return ret;
    //}
    return result_bbox;
};

getBoundingBoxForSingleShapeMatrix = function (this_matrix, ignore_root_matrix, half_size_x, half_size_y, half_size_z) {

    var result_bbox = { min: {x: Infinity, y: Infinity, z: Infinity}, max: {x: -Infinity, y: -Infinity, z: -Infinity}};
    //var updatedData = decomposeMatrix(this_matrix);
    var matrix4 = new THREE.Matrix4();
    // ignore_root_matrix = true;

    if(!ignore_root_matrix) {
        THREE.Matrix4.prototype.set.apply(matrix4, this_matrix);
    }
    //console.error(this_matrix);
    if(!half_size_x){
        //half_size = 0.5;
        console.error("half_size not specified");
    }
    var start_x = -half_size_x;
    var stop_x =  +half_size_x;
    var start_y = -half_size_y;
    var stop_y =  +half_size_y;
    var start_z = -half_size_z;
    var stop_z =  +half_size_z;
    var v = new THREE.Vector3(0, 0, 0);
    var new_x, new_y, new_z;
    for(var x=start_x; x < stop_x + 0.01; x += stop_x - start_x) {
        for(var y=start_y; y < stop_y + 0.01; y += stop_y - start_y) {
            for(var z=start_z; z < stop_z + 0.01; z += stop_z - start_z) {

                v.x = x;
                v.y = y;
                v.z = z;
                v.applyMatrix4(matrix4);
                new_x = v.x;
                new_y = v.y;
                new_z = v.z;

                if(new_x < result_bbox.min.x) result_bbox.min.x = new_x;
                if(new_y < result_bbox.min.y) result_bbox.min.y = new_y;
                if(new_z < result_bbox.min.z) result_bbox.min.z = new_z;
                if(new_x > result_bbox.max.x) result_bbox.max.x = new_x;
                if(new_y > result_bbox.max.y) result_bbox.max.y = new_y;
                if(new_z > result_bbox.max.z) result_bbox.max.z = new_z;
            }

        }

    }
    return result_bbox;
};

/* Multiplies the this_matrix into input_sbbox and updates result_bbox
*/
updateBoundingBoxForSubBoundingBox = function (result_bbox, this_matrix, input_sbbox, ignore_root_matrix) {
    if(result_bbox == null) {
        result_bbox = { min: {x: Infinity, y: Infinity, z: Infinity}, max: {x: -Infinity, y: -Infinity, z: -Infinity}};
    }

    var matrix4 = new THREE.Matrix4();
    // ignore_root_matrix = true;
    if(!ignore_root_matrix) {
        THREE.Matrix4.prototype.set.apply(matrix4, this_matrix);
    }
    var new_x, new_y, new_z;
    var v = new THREE.Vector3(0, 0, 0);
    var x, y, z;
    for(var xi=0; xi <= 1; xi++) {
        x = (xi==0) ? (input_sbbox.min.x) : (input_sbbox.max.x);

        for(var yi=0; yi < 2; yi++) {
            y = (yi==0) ? (input_sbbox.min.y) : (input_sbbox.max.y);

            for(var zi=0; zi < 2; zi++) {
                z = (zi==0) ? (input_sbbox.min.z) : (input_sbbox.max.z);

                //console.log(v);
                v.x = x;
                v.y = y;
                v.z = z;
                v.applyMatrix4(matrix4);
                new_x = v.x;
                new_y = v.y;
                new_z = v.z;

                if(new_x < result_bbox.min.x) result_bbox.min.x = new_x;
                if(new_y < result_bbox.min.y) result_bbox.min.y = new_y;
                if(new_z < result_bbox.min.z) result_bbox.min.z = new_z;
                if(new_x > result_bbox.max.x) result_bbox.max.x = new_x;
                if(new_y > result_bbox.max.y) result_bbox.max.y = new_y;
                if(new_z > result_bbox.max.z) result_bbox.max.z = new_z;
            }

        }

    }
    //console.log("tmp bbox :" + JSON.stringify(result_bbox));
    return result_bbox;
};



function get_bbox_centre(bbox) {
    // todo: reuse the object (avoid gc)
    return new Vector3D(
            (bbox.min.x + bbox.max.x) / 2,
            (bbox.min.y + bbox.max.y) / 2,
            (bbox.min.z + bbox.max.z) / 2
        );
}

function bbox_to_list(input_sbbox) {
    var list = [];
    var x, y, z;
    for(var xi=0; xi <= 1; xi++) {
        x = (xi==0) ? (input_sbbox.min.x) : (input_sbbox.max.x);

        for(var yi=0; yi < 2; yi++) {
            y = (yi==0) ? (input_sbbox.min.y) : (input_sbbox.max.y);

            for(var zi=0; zi < 2; zi++) {
                z = (zi==0) ? (input_sbbox.min.z) : (input_sbbox.max.z);

                list.push(x);
                list.push(y);
                list.push(z);
            }

        }

    }
    return list;
}

function bb_isequal(bb1, bb2) {
    var l1 = bbox_to_list(bb1);
    var l2 = bbox_to_list(bb2);
    for(var i=0; i<l1.length; i++) {
        if (Math.abs(l1[i] - l2[i]) > 0.0000001)
            return false;
    }
    return true;
}


function test_bb_1() {

    var sz = 5;
    var mp5_json = '{"printerSettings":{"name":"test","layerThickness":0.2,"emptyLayer":0,"infillSpace":4,"topThickness":0.6,"paramSamples":75,"speedRate":1000,"circleSpeedRate":1000,"temperature":220,"inAirSpeed":7000,"flowRate":0.035,"criticalLength":35,"retractionSpeed":2400,"retractionLength":5,"shellNumber":3,"material":"PLA 2.85mm","autoZScar":true,"zScarGap":0.5,"critLayerTime":6,"filamentDiameter":2.85},"mp5-version":"0.3","root":{"type":"root","children":[{"type":"Union","protected":false,"children":[{"type":"cylinder","displayColor":{"x":0.9661750055896976,"y":0.8085857086202395,"z":0.41578037212168595},"matrix":[10,0,0,-3.7368777450135866,0,10,0,-1.9559832356144682,0,0,10,1.7323194345664206e-7,0,0,0,1],"index":7575510},{"type":"cube","displayColor":{"x":0.23399071378141634,"y":0.31584816496653323,"z":0.35457351563365425},"matrix":[10,0,0,1.867397834493545,0,10,0,-1.7325527119763677,0,0,10,-9.734106853898084e-10,0,0,0,1],"index":5587759},{"type":"cylinder","displayColor":{"x":0.43814645627496795,"y":0.39556472441055845,"z":0.3415798414286939},"matrix":[10,0,0,1.8694799105200275,0,10,0,3.688535947590836,0,0,10,-1.7225853365943067e-7,0,0,0,1],"index":6657333}],"initialSize":{"x":1,"y":1,"z":1},"displayColor":{"x":0.6470588235294118,"y":0.2784313725490196,"z":0.5882352941176471},"matrix":['+sz+',0,0,82.63768850593796,0,'+sz+',0,126.37324151118989,0,0,'+sz+',5.000000079354265,0,0,0,1],"index":4872526}]}}';

    var mp5_dict = JSON.parse(mp5_json);
    console.log(mp5_dict);
    //console.log(mp5_dict.root.sons[0]);
    var model = dictToShapesOrTreeOrModel(mp5_dict);
    var shape = model.root.sons[0];
    console.log(shape);
    var bb = getBoundingBoxForTree(shape, false);
    console.log(bb);

    var correct_bb__ignore_root_matrix =
        {min:{"x":38.95329928398132,"y":91.59332883358002,"z":-20.00000086129269},"max":{"x":116.98508715629578,"y":169.8159248828888,"z":30.000000866159695}};
    assert(bb_isequal(bb, correct_bb__ignore_root_matrix));

    assert(bb_isequal(
        getBoundingBoxForTree(shape, true),
        {"min":{"x":-8.736877679824829,"y":-6.955983281135559,"z":-5.000000172258538},"max":{"x":6.8694798946380615,"y":8.688535928726196,"z":5.000000173231939}}
    ));

    console.log("Tests OK.")
}


// see tree.js, shape.js

