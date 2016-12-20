ImpliSolid
==========
A library for modelling and polygonisation of solids using the Implicit Surfaces approach.
The original purpose is to use for 3D printing.

But it can be used for research on implicit modeling, as well as visualisation and artistic purposes.

ImpliSolid is the solid modelling engine behind WeDesign.Live, a web-based solid designer for 3D printing.
WeDesign.Live uses MP5, the new collaborative standard for 3D printing. ImpliSolid can use format MP5 for modelling the shapes.

ImpliSolid is free software and open source under GPL3 license.

[plan.md](./plan.md)


Usage:

```javascript
var heart_shape = {"type":"Union","children":[{"type":"cylinder","matrix":[10,0,0,-3.73,0,10,0,-1.955,0,0,10,0,0,0,0,1]},{"type":"cube","matrix":[10,0,0,1.867,0,10,0,-1.732,0,0,10,0,0,0,0,1]},{"type":"cylinder","matrix":[10,0,0,1.869,0,10,0,3.688,0,0,10,0,0,0,0,1]}],"matrix":[10.0,0,0,82.637,0,10.0,0,126.373,0,0,10.0,5.0,0,0,0,1]};
var polygonization_settings;

IMPLICIT_WORKER.getLiveGeometry_from_json(
    heart_shape, polygonization_settings={},
    function (geom, shape_id){
        mesh = new THREE.Mesh( geom, new THREE.MeshBasicMaterial({}) );
        scene.add( mesh );
    }
);
```
