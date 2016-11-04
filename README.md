ImpliSolid
=========
**ImpliSolid** is a Solid Modelling library/engine based on *Implicit Surfaces* modelling (aka *F-REP*).

ImpliSolid uses academic research such as the polygonization algorithm published by Ohtake & Belyaev.

ImpliSolid uses very efficient vectorised calculations to provide instant polygonisation of Implicit Surfaces on your browser.

Currently two open-source projects that use this library:

* [mp5-private](http://github.com/sohale/mp5-private), i.e. the [WeDesign.Live](http://beta.wedesign.live) (incubated open-source)
* [mp5slicer](http://github.com/sohale/mp5slicer) A slicer for 3D printing (incubated open-source)

#### Architecture
The Library has three levels, each with a separate API. The final version will have 5 levels.
Levels 1 and 2 are independent of ThreeJS, hence can be used in NodeJS or as a WebWorker.
Only the Highest level API (level 3) uses Three.JS (For example see: mcc2_3js_r79.js ).

## Implemenation

ImpliSolid is implemented for JavaSCript, C++, Python.
### JavaScript
The JS part can run on **Browser**, or NodeJS, or as a Web Worker.

Implisolid was originally intended for use in Browsers without need to any serverside.

ImpliSolid will soon be available in **WebAssembly** (tutorial comming soon).

The **npm** will be published soon.

### C++
ImpliSolid is written in C++.

The C++ implementation can be used standalone for modern C++ compilers such as g++. It is compiled into JavaScript using *Emscripten*.
### Python
Two Python implementations are already available:
1. As a native Python implementation
2. As a Python binding.



##Available targets
List of platforms with scripts (to install and compile):

- JavaScript (asm.js, on browser)
- JavaScript npm.
- JavaScript, WebWorker version
- WebAssembly (comming soon)
- JavaScript (browser) using Docker
- Python binding (recommended)
- Python native (installs dependencies such as VTK)
- Clang/LLVM
- Native C++: g++

#Tutorials:
> - [Prerequisites](#prerequisites)
- [File organisation](#file-organisation)
- [Python](#python)
- [C++ to javascript](#c-to-javascript)
    - [Boost](#boost)
    - [Emscripten](#emscripten)
    - [How to compile](#how-to-compile)

### Deployment on web server:

This document describes how to set up your installation so that you can work on the solid modeler projects.

Specifically, instructions for installing MP5/WeDesign Solid Modeler on web server:

**Table of Contents**

> - [Prerequisites](#prerequisites)
- [File organisation](#file-organisation)
- [Python](#python)
- [C++ to javascript](#c-to-javascript)
    - [Boost](#boost)
    - [Emscripten](#emscripten)
    - [How to compile](#how-to-compile)


----------

Prerequisites
-------------
Basic requirements:

 - A Linux distribution of your choice (everything can work on an other os but these tutorials will focus on linux)
 - Git
 - Bash
 - NodeJS
 - A text editor
 - C++ compiler (g++ is a good choice and is often preinstalled)


Installation and proper configuration of those software is out of the scope of this document. Please refer to their documentation.

----------

File organisation
-------------------

As of today there are three main folders to this project :

- <b>Clean_code</b> contains a for now frozen implementation of the Dual/Primal Mesh Optimization algorithme for Polygonized Implicit Surfaces. You can learn more about this algorithme here : http://www.hyperfun.org/SM02ob.pdf. It is the main output of the solid modeler project.
- <b>implicit</b> contains still worked on python scripts related to the project. A big part of this file has been transfered to clean_code.
- <b>js_iteration_1</b> contains the c++ implementation of the projet, as well as its transcription to javascript and html.

-------------------

Python
-------------------

In order to work on and execute the python scripts you may need to install the following libraries :
- apdb : the debugging process is made easier
- vtk : It is mainly used to compute a marching cube mesh, the first step of the main algorithm described above. You may want to follow this tutorial to install vtk : http://www.vtk.org/Wiki/VTK/Building/Linux.
- numpy : mathematical functions
- traits
- mayavi : display. This tutorial can help you install mayavi and traits and check wether your installation was succesfull : http://docs.enthought.com/mayavi/mayavi/installation.html

For most of those libraries, a simple `sudo apt-get install python-LibraryName` will do the trick. If not you may want to try the `easy_install` command or the use of pip.

-------------------

C++ to javascript
-------------------
In this part of the project, Emscripten is used to convert c++ code into javascript code that is later used in html.
#### <b>Boost</b>
For the c++ code, we use the Boost library. It does not need to be built but still need to be installed. You can follow this tutorial to do so : http://www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html.

#### <b>Emscripten</b>
The installation of Emscripten often proves a little trickier than the other installations. Here are two tutorials you should follow <b>in the order proposed</b> : https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html#platform-notes-installation-instructions-portable-sdk  http://kripken.github.io/emscripten-site/docs/building_from_source/building_fastcomp_manually_from_source.html#building-fastcomp-from-source.

#### How to compile
The combined use of the Boost library and Emscripten makes compiling a little bit different. For this, we have created the `mc_name.sh` files. For exemple, I use the `mc_marc.sh` file to compile a c++ file into a js file (that is later called in the html file you will simply need to launch in your browser). In this `mc_name.sh` file you'll only find on uncommented line : it is used to compile your c++ file and goes like this :
`~/FolderComtainingYourEmsdk_PortableFolder/emsdk_portable/emscripten/master/em++ -I /LocationOfYourBoostInstalation/boost_1_57_0/ -s EXPORTED_FUNCTIONS="['_make_object', '_main']" -s NO_EXIT_RUNTIME=1 -s ASSERTIONS=1 -pedantic -std=c++14 mc2_sol.cpp -o mc2_sol.cpp.js`. In this exemple you'll compile the file mc2_sol.cpp into mc2_sol.cpp.js exporting its important functions (make_object and main).

-----------

##### Compiling C++ to javascript using Docker
You can use a pre installed socker container
```docker pull mp51/solidmodel  ```
(to be continued)
----------

