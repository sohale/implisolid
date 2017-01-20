////////////////////////////////////////////////////////////
// util
var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;
function assert(condition, message) {
    if (!condition) {
        message = message || "Assertion failed";
        if (typeof Error !== "undefined") {
            throw new Error(message);
        }
        throw message; // Fallback
    }
}

function bbox_center(xmin, xmax, ymin, ymax, zmin, zmax) {
  assert(xmin!=xmax, "wrong usage xmin == xmax");
  assert(ymin!=ymax, "wrong usage xmin == xmax");
  return [(xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2];
}

function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}

// modelling technique

function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}

function subtract(val_0, val_1) {
  return Math.min(val_0, -val_1);
}


function intersect(val_0, val_1) {
  return Math.min(val_0, val_1);
}

function displacement(x, y, z) {
  // usage displacement(x, y, z) + shape(x, y, z) see displacement circle below
  return Math.sin(10*x) + Math.sin(10*y) + Math.sin(10*z);
}

function twist(x, y, z) {
    // usage : change the coordinates x, y, z
    var c = Math.cos(4*y);
    var s = Math.sin(4*y);
    return [c*x - s*z, s*x + c*z, y];
}

// exponential smooth max (k = 32);
var smooth_union = function(val_0, val_1)
{
    // console.log("smooth_union...");
    // console.log(val_0, val_1);
    var k = 8;
    var res = Math.exp( k*val_0 ) + Math.exp( k*val_1 );
    res = -(-Math.log( res )/k);
    // console.log(res);
    return res;
}

// based shape
function circle(x, y, z, radius_square, center_x, center_y, center_z) {
    var x = x - center_x;
    var y = y - center_y;
    var z = z - center_z;
    return radius_square - (x*x + y*y + z*z);
}

///////////////////////////////////////////////////////
// iterative shape
var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;

// utilities
function assert(condition, message) {
    if (!condition) {
        message = message || "Assertion failed";
        if (typeof Error !== "undefined") {
            throw new Error(message);
        }
        throw message; // Fallback
    }
}

function bbox_center(xmin, xmax, ymin, ymax, zmin, zmax) {
  assert(xmin!=xmax, "wrong usage xmin == xmax");
  assert(ymin!=ymax, "wrong usage xmin == xmax");
  return [(xmin + xmax)/2, (ymin + ymax)/2, (zmin + zmax)/2];
}

function length(x, y, z) {
  if (z === undefined) {
    return sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return sqrt(x*x + y*y + z*z);
  }
}

// modelling technique

function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}

function subtract(val_0, val_1) {
  return Math.min(val_0, -val_1);
}


function base_circle(x, y, z, radius) {
  return radius - (x*x + y*y + z*z);
}

// float sdBox( vec3 p, vec3 b )
// {
//   vec3 d = abs(p) - b;
//   return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
// }

function base_cube(x, y, z, base_radius) {

  var cube_x = base_radius;
  var cube_y = base_radius;
  var cube_z = base_radius;

  var d_x = abs(x) - cube_x;
  var d_y = abs(y) - cube_y;
  var d_z = abs(z) - cube_z;

  return  0.01 - (min(max(d_x,max(d_y,d_z)),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0),
                                                  max(d_z,0.0)));
}


<!-- function cross(x, y, z, radius, center_x, center_y, center_z) {

  var d_x = abs(x) - center_x;
  var d_y = abs(y) - center_y;
  var d_z = abs(z) - center_z;

  // var d_x = d_x - radius;
  // var d_y = d_y - radius;
  // var d_z = d_z - radius;

  var x_y = radius - (min(max(d_x,d_y),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0)));

  var x_z = radius - (min(max(d_x,d_z),0.0) + length(max(d_x,0.0),
                                                  max(d_z,0.0)));

  var y_z = radius - (min(max(d_y,d_z),0.0) + length(max(d_y,0.0),
                                                  max(d_z,0.0)));
  return union(union(x_y, x_z), y_z);
} -->

function cross(x, y, z, radius, center_x, center_y, center_z) {

    var new_x = x - center_x;
    var new_y = y - center_y;
    var new_z = z - center_z;

    var cylinder_x_y = radius - (new_x*new_x + new_y*new_y); // your powerful stuff!
    var cylinder_x_z = radius - (new_x*new_x + new_z*new_z);
    var cylinder_y_z = radius - (new_y*new_y + new_z*new_z);

    var res = union(cylinder_x_y, cylinder_x_z);
    return union(res, cylinder_y_z);
    // return union(union(cylinder_x_y, cylinder_x_z), cylinder_y_z);
}

function iterative_center(xmin, xmax, ymin, ymax, zmin, zmax, iter_num, this_cross_radius) {
  if (iter_num===0) {
    return;
  } else {
   iter_num -= 1;
   this_cross_radius -= cross_radius_decrease_per_inter;
  }

  var outer_bbox_cen = bbox_center(xmin, xmax, ymin, ymax, zmin, zmax);
  var x_mid = outer_bbox_cen[0]; // naming for readability
  var y_mid = outer_bbox_cen[1]; // ditto
  var z_mid = outer_bbox_cen[2]; // ditto

  // console.log("num", iter_num, "this_cross_radius", this_cross_radius);
  // console.log(x_mid, y_mid, z_mid);

  assert(this_cross_radius > 0, "cross radius <= 0 ? " + this_cross_radius);
  implicit_value = subtract(implicit_value,
        cross(x, y, z, this_cross_radius, x_mid, y_mid, z_mid));
//   console.log(xmin, xmax, ymin, ymax);

  /*
              xmax/ymax/zmax
          2  3
        x_mid/y_mid/z_mid
          0  1
    xmin/ymin/z_min
  */

  var bbox_0_xmin = xmin;
  var bbox_0_xmax = x_mid;
  var bbox_0_ymin = ymin
  var bbox_0_ymax = y_mid;
  var bbox_0_zmin = zmin;
  var bbox_0_zmax = z_mid;

  var bbox_1_xmin = x_mid;
  var bbox_1_xmax = xmax;
  var bbox_1_ymin = ymin;
  var bbox_1_ymax = y_mid;
  var bbox_1_zmin = zmin;
  var bbox_1_zmax = z_mid;

  var bbox_2_xmin = xmin;
  var bbox_2_xmax = x_mid;
  var bbox_2_ymin = y_mid;
  var bbox_2_ymax = ymax;
  var bbox_2_zmin = zmin;
  var bbox_2_zmax = z_mid;

  var bbox_3_xmin = x_mid;
  var bbox_3_xmax = xmax;
  var bbox_3_ymin = y_mid;
  var bbox_3_ymax = ymax;
  var bbox_3_zmin = zmin;
  var bbox_3_zmax = z_mid;

  var bbox_4_xmin = xmin;
  var bbox_4_xmax = x_mid;
  var bbox_4_ymin = ymin
  var bbox_4_ymax = y_mid;
  var bbox_4_zmin = z_mid;
  var bbox_4_zmax = zmax;

  var bbox_5_xmin = x_mid;
  var bbox_5_xmax = xmax;
  var bbox_5_ymin = ymin;
  var bbox_5_ymax = y_mid;
  var bbox_5_zmin = z_mid;
  var bbox_5_zmax = zmax;

  var bbox_6_xmin = xmin;
  var bbox_6_xmax = x_mid;
  var bbox_6_ymin = y_mid;
  var bbox_6_ymax = ymax;
  var bbox_6_zmin = z_mid;
  var bbox_6_zmax = zmax;

  var bbox_7_xmin = x_mid;
  var bbox_7_xmax = xmax;
  var bbox_7_ymin = y_mid;
  var bbox_7_ymax = ymax;
  var bbox_7_zmin = z_mid;
  var bbox_7_zmax = zmax;

  iterative_center(bbox_0_xmin, bbox_0_xmax, bbox_0_ymin, bbox_0_ymax, bbox_0_zmin, bbox_0_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_1_xmin, bbox_1_xmax, bbox_1_ymin, bbox_1_ymax, bbox_1_zmin, bbox_1_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_2_xmin, bbox_2_xmax, bbox_2_ymin, bbox_2_ymax, bbox_2_zmin, bbox_2_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_3_xmin, bbox_3_xmax, bbox_3_ymin, bbox_3_ymax, bbox_3_zmin, bbox_3_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_4_xmin, bbox_4_xmax, bbox_4_ymin, bbox_4_ymax, bbox_4_zmin, bbox_4_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_5_xmin, bbox_5_xmax, bbox_5_ymin, bbox_5_ymax, bbox_5_zmin, bbox_5_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_6_xmin, bbox_6_xmax, bbox_6_ymin, bbox_6_ymax, bbox_6_zmin, bbox_6_zmax, iter_num, this_cross_radius);
  iterative_center(bbox_7_xmin, bbox_7_xmax, bbox_7_ymin, bbox_7_ymax, bbox_7_zmin, bbox_7_zmax, iter_num, this_cross_radius);

}

var iter_num = 3;
var base_radius = 1;
var cross_radius = 0.2;
var cross_radius_decrease_per_inter = 0.06;

assert(cross_radius - cross_radius_decrease_per_inter*iter_num >0);

// var implicit_value = base_circle(x, y, z, base_radius);
var implicit_value = base_cube(x, y, z, base_radius);


iterative_center(-1,
                  1,
                 -1,
                  1,
                 -1,
                  1,
                 iter_num,
                 cross_radius);

// return cross(x, y, z, 0.1, 0, 0, 0);
return implicit_value;

//////////////////////////////////////////////////////////////

// rugby ball
 // tiger
function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}


var center_x = 1;
var center_y = 1;
var center_z = 1;

var d_x = Math.abs(x) - center_x;
var d_y = Math.abs(y) - center_y;
var d_z = Math.abs(z) - center_z;

return - (Math.min(Math.max(d_x, d_y), 0) + length(Math.abs(x), Math.abs(y), Math.abs(z)));

// }

///////////////////////////////////////////////////////////////
// inflated cube

 // tiger
function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}


var center_x = 1;
var center_y = 1;
var center_z = 1;

var d_x = Math.abs(x) - center_x;
var d_y = Math.abs(y) - center_y;
var d_z = Math.abs(z) - center_z;

return 1 - (Math.min(Math.max(d_x, Math.max(d_y, d_z)), 0) + length(Math.abs(x), Math.abs(y), Math.abs(z)));

// }



///////////////////////////////////


var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;

  function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}


  function length(x, y, z) {
  if (z === undefined) {
    return sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return sqrt(x*x + y*y + z*z);
  }
}

  var base_radius = 1;

  var cube_x = base_radius;
  var cube_y = base_radius;
  var cube_z = base_radius;

  var d_x = abs(x) - cube_x;
  var d_y = abs(y) - cube_y;
  var d_z = abs(z) - cube_z;

  var x_y = 0. - (min(max(d_x,d_y),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0)));


  var x_z = 0.1 - (min(max(d_x,d_z),0.0) + length(max(d_x,0.0),
                                                  max(d_z,0.0)));


return union(x_y, x_z);
//////////////////////////////////////////////////////
// heart
x_0 = x;
x_1 = y;
x_2 = z;
implicit_value = 0.00001 - (Math.pow((x_0*x_0 + (9/4)*x_1*x_1 + x_2*x_2 - 1), 3) - x_0*x_0*x_2*x_2*x_2 - (9/80)*x_1*x_1*x_2*x_2*x_2 ); // your powerful stuff!
return implicit_value;

/////////////////
// displacement circle
// tiger
function displacement(x, y, z) {
    return Math.sin(10*x) + Math.sin(10*y) + Math.sin(10*z);
}

function circle(x, y, z, radius_square, center_x, center_y, center_z) {
    var x = x - center_x;
    var y = y - center_y;
    var z = z - center_z;
    return radius_square - (x*x + y*y + z*z);
}

var cir = circle(x, y, z, 1, 0, 0, 0);
var dis = displacement(x, y, z);
return cir + dis;
///////////

// twist demo

function twist(x, y, z) {
    // changing x, y, z inplace
    var c = Math.cos(4*y);
    var s = Math.sin(4*y);
    return [c*x - s*z, s*x + c*z, y];
}
function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}

new_coor = twist(x, y, z);
x = new_coor[0];
y = new_coor[1];
z = new_coor[2];

function sd_torus(x, y, z, para_outer, para_inner) {
//    vec2 q = vec2(length(p.xz)-t.x,p.y);
//   return length(q)-t.y;
    var q_0 = length(x, z) - para_outer;
    var q_1 = y;
    return 0.5 - length(q_0 - para_inner, q_1 - para_inner);
}

return sd_torus(x, y, z, 0.5, 0.5);

////////////////////////////////////////////////
// smooth union

function circle_left(x, y, z) {
    var x = x - 0.5;
    return 0.5 - (x*x + y*y + z*z);
}


function circle_right(x, y, z) {
    var x = x + 0.5;
    return 0.5 - (x*x + y*y + z*z);
}

var smooth_union = function(val_0, val_1)
{
    // console.log("smooth_union...");
    // console.log(val_0, val_1);
    var k = 8;
    var res = Math.exp( k*val_0 ) + Math.exp( k*val_1 );
    res = Math.log( res )/k;
    // console.log(res);
    return res;
}


return smooth_union(circle_left(x, y, z), circle_right(x, y, z));

/////////////////////////////////////
// random nice shape launching rocker

function displacement(x, y, z) {
    // usage displacement(x, y, z) + shape(x, y, z) see displacement circle below
    if (y < - 1)
        return Math.sin(10*x) + Math.sin(10*y) + Math.sin(10*z);
    else
        return 0;
}

var dis = displacement(x, y, z);

function twist(x, y, z) {
    // changing x, y, z inplace
    var c = Math.cos(2*y);
    var s = Math.sin(2*y);
    return [c*x - s*z, s*x + c*z, y];
}

new_coor = twist(x, y, z);
x = new_coor[0];
y = new_coor[1];
z = new_coor[2];

function circle_left(x, y, z) {
    var x = x - 0.5;
    return 0.5 - (x*x + y*y + z*z);
}


function circle_right(x, y, z) {
    var x = x + 0.5;
    return 0.5 - (x*x + y*y + z*z);
}

function circle_y_left(x, y, z) {
    var y = y - 0.5;
    return 0.5 - (x*x + y*y + z*z);
}


function circle_y_right(x, y, z) {
    var y = y + 0.5;
    return 0.5 - (x*x + y*y + z*z);
}

function circle_z_left(x, y, z) {
    var z = z - 0.5;
    return 0.5 - (x*x + y*y + z*z);
}


function circle_z_right(x, y, z) {
    var z = z + 0.5;
    return 0.5 - (x*x + y*y + z*z);
}


function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}

res = union(circle_left(x, y, z), circle_right(x, y, z));
res = union(circle_y_left(x, y, z), res);
res = union(circle_y_right(x, y, z), res);

res = union(circle_z_left(x, y, z), res);
res = union(circle_z_right(x, y, z), res);

return res + dis;
//////////////////////////////


function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}


var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;

function twist(x, y, z) {
    // changing x, y, z inplace
    var c = Math.cos(2*y);
    var s = Math.sin(2*y);
    return [c*x - s*z, s*x + c*z, y];
}

new_coor = twist(x, y, z);
x = new_coor[0];
y = new_coor[1];
z = new_coor[2];

function base_cube(x, y, z, base_radius) {

  var cube_x = base_radius;
  var cube_y = base_radius;
  var cube_z = base_radius;

  var d_x = abs(x) - cube_x;
  var d_y = abs(y) - cube_y;
  var d_z = abs(z) - cube_z;

  return  0.5 - (min(max(d_x,max(d_y,d_z)),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0),
                                                  max(d_z,0.0)));
}


return base_cube(x, y, z, 1);

/////////////////////////////////////////////////
/// Beyblade


var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;

function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}

// float sdCone( vec3 p, vec2 c )
// {
//    // c must be normalized
//    float q = length(p.xy);
//    return dot(c,vec2(q,p.z));
// }

function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}



function intersect(val_0, val_1) {
  return Math.min(val_0, val_1);
}

function sdCone(x, y, z, c_0, c_1) {
    var q = length(x, z);
    return 0.00001 - (c_0 * q + c_1 * -y);
}

function base_cube(x, y, z, base_radius) {
  var cube_x = base_radius;
  var cube_y = base_radius;
  var cube_z = base_radius;

  var d_x = abs(x) - cube_x;
  var d_y = abs(y) - cube_y;
  var d_z = abs(z) - cube_z;

  return  0.001 - (min(max(d_x,max(d_y,d_z)),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0),
                                                  max(d_z,0.0)));
}



var cone = intersect(sdCone(x, y + 2.7, z, 0.5, 0.5),  base_cube(x, y + 2.7, z, 1));
var cone_down = intersect(sdCone(x, - (y) + 2.7, z, 0.5, 0.5),  base_cube(x, - (y) + 2.7, z, 1));

var cone = union(cone, cone_down);


function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}


var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;

function twist(x, y, z) {
    // changing x, y, z inplace
    var c = Math.cos(2*y);
    var s = Math.sin(2*y);
    return [c*x - s*z, s*x + c*z, y];
}

new_coor = twist(x, y, z);
x = new_coor[0];
y = new_coor[1];
z = new_coor[2];

function base_cube(x, y, z, base_radius) {
  var cube_x = base_radius;
  var cube_y = base_radius;
  var cube_z = base_radius;

  var d_x = abs(x) - cube_x;
  var d_y = abs(y) - cube_y;
  var d_z = abs(z) - cube_z;

  return  0.5 - (min(max(d_x,max(d_y,d_z)),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0),
                                                  max(d_z,0.0)));
}


var twisted_cube = base_cube(x, y, z, 1);

var smooth_union = function(val_0, val_1)
{
    // console.log("smooth_union...");
    // console.log(val_0, val_1);
    var k = 8;
    var res = Math.exp( k*val_0 ) + Math.exp( k*val_1 );
    res = -(-Math.log( res )/k);
    // console.log(res);
    return res;
}

return smooth_union(twisted_cube, cone)
<!-- return cone; -->
//////////////////////////////////////////
// pillar

function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}

// float opRep( vec3 p, vec3 c )
// {
//     vec3 q = mod(p,c)-0.5*c;
//     return primitve( q );
// }

function repeat(x, y, z, c_x, c_y, c_z) {
    var q_x = x % c_x - 0.5 * c_x;
    var q_y = y % c_y - 0.5 * c_y;
    var q_z = z % c_z - 0.5 * c_z;
    return [q_x, q_y, q_z];
}

var new_coor = repeat(x, y, z, 2, 2, 2);

x = new_coor[0];
y = new_coor[1];
z = new_coor[2];

x = x*x;
y = y*x;
function sd_torus(x, y, z, para_outer, para_inner) {
    var q_0 = length(x, z) - para_outer;
    var q_1 = y;
    return 0.5 - length(q_0 - para_inner, q_1 - para_inner);
}


function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}
return sd_torus(x, y, z, 0.2, 0.1);

///////////////////////////////
// space ship


var min = Math.min;
var max = Math.max;
var sqrt = Math.sqrt;
var abs = Math.abs;

function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}

function length(x, y, z) {
  if (z === undefined) {
    return Math.sqrt(x*x + y*y);
  } else {
    if (x === undefined | y === undefined) {
      throw SyntaxError('x or y should not be undefined');
    }
    return Math.sqrt(x*x + y*y + z*z);
  }
}

var smooth_union = function(val_0, val_1)
{
    // console.log("smooth_union...");
    // console.log(val_0, val_1);
    var k = 8;
    var res = Math.exp( k*val_0 ) + Math.exp( k*val_1 );
    res = -(-Math.log( res )/k);
    // console.log(res);
    return res;
}


function sd_ellipsoid(x, y, z, r_x, r_y, r_z, size, center_x, center_y, center_z) {
    return size - length(x/r_x + center_x, y/r_y + center_y, z/r_z + center_z) * min(min(r_x, r_y), r_z);
}

function base_cube(x, y, z, base_radius, center_x, center_y, center_z) {

  var x = x + center_x;
  var y = y + center_y;
  var z = z + center_z;

  var cube_x = base_radius;
  var cube_y = base_radius;
  var cube_z = base_radius;

  var d_x = abs(x) - cube_x;
  var d_y = abs(y) - cube_y;
  var d_z = abs(z) - cube_z;

  return  0.5 - (min(max(d_x,max(d_y,d_z)),0.0) + length(max(d_x,0.0),
                                                  max(d_y,0.0),
                                                  max(d_z,0.0)));
}

fist_l = base_cube(x, y, z, 0.1, 2.7, -2.5, 0);
fist_r = base_cube(x, y, z, 0.1, -2.7, -2.5, 0);

ellipsoid_head = sd_ellipsoid(x, y, z,
                              0.7, 0.7, 0.7,
                              1.5,
                              0, 3.5, 0);

ellipsoid_body = sd_ellipsoid(x, y, z,
                              3, 4, 3,
                              2,
                              0, 0, 0);

ellipsoid_shoulder = sd_ellipsoid(x, y, z,
                              3.5, 1, 1,
                              1,
                              0, 1, 0);

ellipsoid_hand_l = sd_ellipsoid(x, y, z,
                              0.5, 2, 0.5,
                              0.4,
                              5.5, -0.5, 0);

ellipsoid_hand_r = sd_ellipsoid(x, y, z,
                              0.5, 2, 0.5,
                              0.4,
                              -5.5, -0.5, 0);

ellipsoid_leg_l = sd_ellipsoid(x, y, z,
                              0.5, 1, 0.5,
                              1,
                              5.5, -4.5, 0);

ellipsoid_leg_r = sd_ellipsoid(x, y, z,
                              0.5, 1, 0.5,
                              1,
                              -5.5, -4.5, 0);

res = ellipsoid_head;
res = smooth_union(res, ellipsoid_body);
res = smooth_union(res, ellipsoid_shoulder);
res = smooth_union(res, ellipsoid_hand_l);
res = smooth_union(res, ellipsoid_hand_r);
res = smooth_union(res, fist_l);
res = smooth_union(res, fist_r);
res = smooth_union(res, ellipsoid_leg_l);
res = smooth_union(res, ellipsoid_leg_r);


function twist(x, y, z) {
    // changing x, y, z inplace
    var c = Math.cos(2*y);
    var s = Math.sin(1*y);
    return [c*x - s*z, s*x + c*z, y];
}

new_coor = twist(x, y, z);
twisted_x = new_coor[0];
twisted_y = new_coor[1];
twisted_z = new_coor[2];

function sd_torus(x, y, z, para_outer, para_inner) {
//    vec2 q = vec2(length(p.xz)-t.x,p.y);
//   return length(q)-t.y;
    var q_0 = length(x, z) - para_outer;
    var q_1 = y;
    return 0.5 - length(q_0 - para_inner, q_1 - para_inner);
}

res = union(res, sd_torus(twisted_x, twisted_y, twisted_z, 1, 1.5));




function displacement(x, y, z) {
    // usage displacement(x, y, z) + shape(x, y, z) see displacement circle below
    if (y < - 1)
        return Math.sin(10*x) + Math.sin(10*y) + Math.sin(10*z);
    else
        return 0;
}

var dis = displacement(x, -y + 2, z);

function twist(x, y, z) {
    // changing x, y, z inplace
    var c = Math.cos(2*y);
    var s = Math.sin(2*y);
    return [c*x - s*z, s*x + c*z, y];
}

var r = 0.2;
function circle_left(x, y, z) {
    var x = x - 0.5;
    return r - (x*x + y*y + z*z);
}


function circle_right(x, y, z) {
    var x = x + 0.5;
    return r - (x*x + y*y + z*z);
}

function circle_y_left(x, y, z) {
    var y = y - 0.5;
    return r - (x*x + y*y + z*z);
}


function circle_y_right(x, y, z) {
    var y = y + 0.5;
    return r - (x*x + y*y + z*z);
}

function circle_z_left(x, y, z) {
    var z = z - 0.5;
    return r - (x*x + y*y + z*z);
}


function circle_z_right(x, y, z) {
    var z = z + 0.5;
    return 0.5 - (x*x + y*y + z*z);
}


function union(val_0, val_1) {
  return Math.max(val_0, val_1);
}

res = union(circle_left(x, y, z), res);
res = union(circle_right(x, y, z), res);
res = union(circle_y_left(x, y, z), res);
res = union(circle_y_right(x, y, z), res);

res = union(circle_z_left(x, y, z), res);
res = union(circle_z_right(x, y, z), res);

res =  res + dis;


function intersect(val_0, val_1) {
  return Math.min(val_0, val_1);
}

function sdCone(x, y, z, c_0, c_1) {
    var q = length(x, z);
    return 0.0001 - (c_0 * q + c_1 * -y);
}


var cone = intersect(sdCone(x, y + 4.5, z, 5, 5),  base_cube(x, y + 4.5, z, 0.5, 0, 0, 0));

return smooth_union(res, cone);

