
## 9 March 2021
* produce_object_old2() is no longer used. IT was exported but the function wass already dremoved from mcc2.cpp
### intents:
* todo: mcc2.cpp should be moved to iteration2
* Remaining hpp in js_iteration1 are: ...
* Note about js_iteration1: Why some hpp files still remain there, should be move js, or cpp there?

* What is the preferred direcrtory structure for cpp versus js files?
  * We need to have the C++ parts standalone, if there is no js / Emscriptedn involved.
  * The interface parts need to be in a separate cpp file. mcc2.cpp, which will contain that, shoul dbe strcitly minimal.

* Use CMake

### Questions
* What was the point about `MORE_ABOUT_INFO`?
    * It is shown on console
* A different way of pulling the libraries `/lib`
