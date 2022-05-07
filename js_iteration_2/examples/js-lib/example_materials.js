'use strict';

function look_nice() {
    var materials = [];

    // Materials
    var texture = new THREE.Texture( generateTexture() );
    texture.needsUpdate = true;
    materials.push( new THREE.MeshLambertMaterial( { map: texture, transparent: true } ) );
    materials.push( new THREE.MeshLambertMaterial( { color: 0xdddddd, shading: THREE.FlatShading } ) );
    materials.push( new THREE.MeshPhongMaterial( { color: 0xdddddd, specular: 0x009900, shininess: 30,
        shading: THREE.FlatShading,
        side:THREE.DoubleSide,
        opacity: 0.99, transparent: true,
        } ) );

    materials.push( new THREE.MeshNormalMaterial( ) );
    materials.push( new THREE.MeshBasicMaterial( { color: 0xffaa00, transparent: true, blending: THREE.AdditiveBlending } ) );
    //materials.push( new THREE.MeshBasicMaterial( { color: 0xff0000, blending: THREE.SubtractiveBlending } ) );
    materials.push( new THREE.MeshLambertMaterial( { color: 0xdddddd, shading: THREE.SmoothShading } ) );
    materials.push( new THREE.MeshPhongMaterial( { color: 0xdddddd, specular: 0x009900, shininess: 30, shading: THREE.SmoothShading, map: texture, transparent: true } ) );
    materials.push( new THREE.MeshNormalMaterial( { color: 0x00dd00,  shading: THREE.SmoothShading,
        //opacity:0.5,  transparent: true, side:THREE.DoubleSide
    } ));
    materials.push( new THREE.MeshBasicMaterial( { color: 0xffaa00, wireframe: true } ) );
    materials.push( new THREE.MeshDepthMaterial() );
    materials.push( new THREE.MeshLambertMaterial( { color: 0x666666, emissive: 0xff0000, shading: THREE.SmoothShading } ) );
    materials.push( new THREE.MeshPhongMaterial( { color: 0x000000, specular: 0x666666, emissive: 0xff0000, shininess: 10, shading: THREE.SmoothShading, opacity: 0.9, transparent: true } ) );
    materials.push( new THREE.MeshBasicMaterial( { map: texture, transparent: true } ) );

    // 13 : Standard WeDesignLive
    materials.push( new THREE.MeshPhongMaterial( { color: 0x9B2D30, transparent:true })  );
    // 14
    materials.push( new THREE.MeshPhongMaterial( { color: 0x9B2D30, transparent:false, shading: THREE.FlatShading })  );
    //opacity: 0.7,

    console.log("Choose material "+(materials.length-1));

    return materials;
}


/**
    Stuff to make things look pretty.
    Taken from http://threejs.org/examples/#webgl_materials
*/
function generateTexture() {
    var canvas = document.createElement( 'canvas' );
    //var canvas = document.getElementById( 'my_canvas' );
    canvas.width = 256;
    canvas.height = 256;
    var context = canvas.getContext( '2d' );
    var image = context.getImageData( 0, 0, 256, 256 );
    var x = 0, y = 0;
    for ( var i = 0, j = 0, l = image.data.length; i < l; i += 4, j ++ ) {
        x = j % 256;
        y = x == 0 ? y + 1 : y;
        image.data[ i ] = 255;
        image.data[ i + 1 ] = 255;
        image.data[ i + 2 ] = 255;
        image.data[ i + 3 ] = Math.floor( x ^ y );
    }
    context.putImageData( image, 0, 0 );
    return canvas;
}
