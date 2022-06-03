
* simple_viewer1
   * The idea is [...]
   * See [examples/simple_viewer1](https://github.com/sohale/implisolid/tree/revival-sohale/examples/simple_viewer1) ([asisnow](https://github.com/sohale/implisolid/blob/d7856f21df8d470a1096fa146a79ba2b410547de/examples/simple_viewer1))

* example browser
   * Will have a panel summoning all examples (interactive demo, viewer, maybe a nodejs example runner, maybe a textual runner of unit tests)
   * See [it](https://github.com/sohale/implisolid/tree/revival-sohale/examples/example_browser) in progress.

* How to embed the architectural diagram in the .md files (readme). It does not show:
```html
<iframe src="https://docs.google.com/presentation/d/e/2PACX-1vTGPTl24WtyxHqmImwH2CqjKBnzvBe4tvjNdZGD1e1uXCJqTb8QgG8wfAIhzHfpNnxQ38Pt93GGiLAB/embed?start=false&loop=false&delayms=3000" frameborder="0" width="1534" height="759" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe>
```
* The node-js -based. Purposes:
   * unit tests
   * Sanity test of implisolid
   * Help decouple (make it indenpendent) from threejs, hence, better organistion.
   * (Somehow compaatible with webworker?)

* Use `autodiff` for symbolic mathematics in C++
   * Latest progress: compiles. [Snippep](https://github.com/sohale/implisolid/blob/d7856f21df8d470a1096fa146a79ba2b410547de/sandbox/autodiff/implicit-functions/cpp/autodiff-sample1.cpp) with a minimal usage example of aautodiff (full e2e script)
   * The idea is to export glsl functions (as an overrided method of the class?)

* A simple `autodiff`-based glsl visualiser (also see `symmpy`)
   * Maybe: online gcc-compile server?
   * Stash / sanadbox of developmental implicit functions (for people with repos? a package manager in style of: atom/npm/templlator)


* Create a demo folder with a glsl demo for some simple SDF
  * First, a simple glsl demo
  * 2, Run a simple SDF on that panel
  * Later on (vision/telos!): develop a sandbox for python-based.
  * Sets the stage for autodiff-basaed and sympy-based.

* A `sympy`-based glsl visualiser
  * Use `sympy` (Python) for exploring new implicit functions (As a branch of investigation.
  * Exports `.glsl`
  * Use separaate C++ autodiff for the other branch: C++-based sybolic mathematics. But Keep the sympy alive.
  * A demo page (view only)
  * (Maybe) a demo page that serves python?
  * autodiff-basaed and sympy-based ecosystems (See Stash/sandbox above).

* Use `@types/emscripten`

* A `tfjs` backend !

* Architctual overhaul:
## Architectural Structure
"Implisolid current architecture: Some proposed changes shown in purple:
![svg]( https://drive.google.com/uc?export=view&id=1VKJaGe-Hycb6qcfEx0rTbSrc6FtcT0G-   "Implisolid current architecture (needs change)" )

Required:
* decouple into independent interoperable parts
    * low level vectorised stuff
    * basic structures
    * the "simplical complexes" layer (includes the mesh).
* More structural changes:
* C++ namespaces for various above-mentioned parts.
* Folder structure to reflect the new structure.
* Clean up: A branch to keep (but remove from main branch) the contributions that are not used: mark, etc. Also the newer ones (already in a separate branch: Tiger and Shane Blobixx )
* Separate the js parts that use ThreeJS (half-way) from the rest.
* NodeJS-based unit tests need to be able to test the web-worker functionality.

