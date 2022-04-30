ImpliSolid
=========

[![Join the chat at https://gitter.im/implisolid/Lobby](https://badges.gitter.im/implisolid/Lobby.svg)](https://gitter.im/implisolid/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

**ImpliSolid** is a Geometric Modelling library suitable for solid omdlling engine based on *Implicit Surfaces* modelling (aka *F-REP*).
The main usecase for ImpliSolid is 3D printing.

ImpliSolid uses very efficientcalculations to provide instant polygonisation of Implicit Surfaces efficient eniough to run on your browser using CPU only.

The main strength is its ability to work efficiently with **sharp edges**.
It also uses **adaptive subdivition** for smooth and perfect curved surfaces.
These are achieved using relatively lower resolution meshes.

It uses "vectorised" numerical calulations to achieve higher speed by utilising Instruction Pipelining in modern CPUs.
This enables it to be useful on consumer and home computers.
<!-- This enables it to be useful on consumer and home computers on browser without GPU.-->

ImpliSolid use is not limited to browsers. It has implementations in C++, Python (native) and JavaScript.

| | |
|------:|:-------|
|  Implisolid homepage: | [homepage](https://sohale.github.io/implisolid/) |
|  Github:  | [sohale/implisolid](https://github.com/sohale/implisolid) |
| An interactive live demo:| [mp5 editor](http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html) |
<!-- |  An interactive editor: | [link defunct](https://api-project-1000362687695.appspot.com/mp5interactive/mp5_json_code.html) | -->


**`ImpliSolid`** uses academic research on the polygonization algorithm such as the algorithm by [Ohtake](https://www.u-tokyo.ac.jp/focus/en/people/people000639.html) & [Belyaev](https://scholar.google.co.uk/citations?user=UgOo39sAAAAJ&hl=en).

<!-- https://dl.acm.org/doi/10.1145/882262.882293 -->

Currently two open-source projects that use this library:

* [mp5-private](http://github.com/sohale/mp5-private), i.e. the [WeDesign.Live](http://beta.wedesign.live) (incubated)
* [mp5slicer](http://github.com/sohale/mp5slicer) A slicer for 3D printing (incubated)

## E2E demo
For single-click execution in your computer (tested on Ubuntu and MacOS), run:
```bash
git clone git@github.com:sohale/implisolid.git
cd implisolid/
bash ./scripts/e2e-test-builds.bash
```
This will run and end-to-end demo: Pulls the code, compiles the code for Emscripten. Then launches a web server and runs a demo on browser [like this](http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html). The interactive live demo instantly polygonises the objects and visualises them as the `mp5` file is edited.




## Implemenation

ImpliSolid is implemented for JavaSCript, C++, Python.
### JavaScript
Implicit can be compiled into JavaScript using *Emscripten*.

The Javascript compile target part can run on **Browser**, or NodeJS, or as a Web Worker.

Implisolid was originally intended for use in Browsers without need to any serverside.

ImpliSolid will soon be available in **WebAssembly** (tutorial comming soon).

The **npm** is available from [here](www.npmjs.com/package/@sohale/implisolid) (will be moved).

#### Architecture
The `IMPLICIT` javascript library is a wrapper for the library for browser: js + [threejs](https://threejs.org).
It has three levels, each with a separate API. The final version will have 5 levels.
Levels 1 and 2 are independent of ThreeJS, hence can be used in NodeJS or as a WebWorker.
Only the Highest level API (level 3) uses Three.JS (For example see: [mp5_json_code.html](js_iteration_2/examples/mp5interactive/mp5_json_code.html) and [demo1-deploy.sh](scripts/demos/demo1/demo1-deploy.sh), [mcc2_3js_r79.js](js_iteration_1/mcc2_3js_r79.html) or [this code fragment](docs/readme.md)).


* A variant version builds objects incrementally.
* A variant builds objects incrementally in Web Worker of browser (see [implisolid_worker.js](js_iteration_2/js/implisolid_worker.js)).

* The main target can be compiled using `scripts/build-emscripten.sh`
* A web-assembly version will be available soon.

### C++
ImpliSolid is written in modern C++.

The C++ implementation can be used standalone for modern C++ compilers such as `g++` and `clang`. The main target is compiled into JavaScript using *Emscripten*.
### Python
Two Python implementations are already available:
1. As a native Python implementation
2. As a Python binding that connects to a binary compilation of the C++ version.


## Available targets
List of platforms with scripts (to install and compile):

- JavaScript (browser, formerly called `asm.js`)
- JavaScript, WebWorker version
- WebAssembly (comming soon)
- JavaScript on `npm` (comming soon)
- JavaScript (browser) using Docker (comming soon)
- Python binding (recommended)
- Python native (requires dependencies such as [VTK](https://vtk.org/doc/nightly/html/md__builds_gitlab_kitware_sciviz_ci_Documentation_Doxygen_PythonWrappers.html) or [Mayavi](https://docs.enthought.com/mayavi/mayavi/))
- Clang / LLVM
- Native C++: `g++` (see [script](js_iteration_1/build_g++.sh))

## Tutorials:

### E2E demo
For single-click execution in your computer (tested on Ubuntu and MacOS), run:
```bash
git clone git@github.com:sohale/implisolid.git
cd implisolid/
bash ./scripts/e2e-test-builds.bash
```
This will run and end-to-end demo: Pulls the code, compiles the code for Emscripten. Then launches a web server and runs a demo on browser [like this](http://sohale.github.io/demos/implisolid-build/demo1/mp5_json_code.html). The interactive live demo instantly polygonises the objects and visualises them as the `mp5` file is edited.

### Deployment on web server

This document describes how to set up your installation so that you can work on the solid modeler projects.

Specifically, instructions for installing MP5/WeDesign Solid Modeler on web server:


#### Prerequisites
-------------
##### Basic requirements:

 - A Linux distribution of your choice (everything can work on an other os but these tutorials will focus on linux)
 - Git
 - Bash
 - NodeJS
 - A text editor
 - C++ compiler (g++ is a good choice and is often preinstalled)


Installation and proper configuration of above software is out of the scope of this document. Please refer to their documentation.

## File organisation
-------------------

As of 2016 there are three main folders to this project :

- <b>Clean_code</b> contains a for now frozen implementation of the Dual/Primal Mesh Optimization algorithme for Polygonized Implicit Surfaces. You can learn more about this algorithme here : http://www.hyperfun.org/SM02ob.pdf. It is the main output of the solid modeler project.
- <b>implicit</b> contains still worked on python scripts related to the project. A big part of this file has been transfered to `python_clean_code`.
- <b>js_iteration_1</b> contains the c++ implementation of the projet, as well as its transcription to javascript and html.

-------------------

### Python
-------------------

In order to work on and execute the python scripts you may need to install the following libraries :
- apdb : the debugging process is made easier
- vtk : It is mainly used to compute a marching cube mesh, the first step of the main algorithm described above. You may want to follow this tutorial to install vtk : http://www.vtk.org/Wiki/VTK/Building/Linux.
- numpy : mathematical functions
- traits
- mayavi : display. This tutorial can help you install mayavi and traits and check wether your installation was succesfull : http://docs.enthought.com/mayavi/mayavi/installation.html

For most of those libraries, a simple `sudo apt-get install python-LibraryName` will do the trick. If not you may want to try the `easy_install` command or the use of pip.

-------------------

### C++ to javascript
-------------------
See [scripts/e2e-test-builds.bash](scripts/e2e-test-builds.bash).
Simply run `bash ./scripts/e2e-test-builds.bash`


#### Older notes are below:


In this part of the project, Emscripten is used to convert c++ code into javascript code that is later used in html.
##### Boost
For the c++ code, we use the Boost library. It does not need to be built but still need to be installed. You can follow this tutorial to do so : http://www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html.

##### Emscripten
Emscripten can be used simply using Docker `emscripten/emsdk`. See abovementioned script.

Old note:
The installation of Emscripten often proves a little trickier than the other installations. Here are two tutorials you should follow <b>in the order proposed</b> : https://kripken.github.io/emscripten-site/docs/getting_started/downloads.html#platform-notes-installation-instructions-portable-sdk  http://kripken.github.io/emscripten-site/docs/building_from_source/building_fastcomp_manually_from_source.html#building-fastcomp-from-source.

##### How to compile
The combined use of the Boost library and Emscripten makes compiling a little bit different. For this, we have created the `mc_name.sh` files. For exemple, I use the `mc_marc.sh` file to compile a c++ file into a js file (that is later called in the html file you will simply need to launch in your browser). In this `mc_name.sh` file you'll only find on uncommented line : it is used to compile your c++ file and goes like this :
`~/FolderComtainingYourEmsdk_PortableFolder/emsdk_portable/emscripten/master/em++ -I /LocationOfYourBoostInstalation/boost_1_57_0/ -s EXPORTED_FUNCTIONS="['_make_object', '_main']" -s NO_EXIT_RUNTIME=1 -s ASSERTIONS=1 -pedantic -std=c++14 mc2_sol.cpp -o mc2_sol.cpp.js`. In this exemple you'll compile the file mc2_sol.cpp into mc2_sol.cpp.js exporting its important functions (make_object and main).

-----------

##### Compiling C++ to javascript using Docker

Simply see [build-emscripten.sh](scripts/build-emscripten.sh).

Old note:

You can use a pre installed docker container to compile implisolid:
```docker run -v /Users/Tiger/Documents/mp5-private/implisolid:/src -t mp51/solidmodel /bin/bash /src/js_iteration_1/build_mcc2_docker.sh -o```

----------
