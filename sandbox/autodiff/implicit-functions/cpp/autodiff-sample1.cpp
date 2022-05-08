// based on emscr1.cpp
#include <iostream>

/* MAkes sure also to enable `-std=c++17` */
// autodiff include
#include <autodiff/forward/dual.hpp>
using namespace autodiff;

  // The single-variable function for which derivatives are needed
  dual f(dual x)
  {
      return 1 + x + x*x + 1/x + log(x);
  }

extern "C" {
    // void build_geometry(const char* shape_parameters_json, const char* mc_parameters_json);
    // void build_geometry_u(const char* shape_parameters_json, const char* mc_parameters_json, const char* call_specs);
    void about();
    void about2();
}
void about() {
    std::cout << "Hello, from emscr2! autodiff 1\n";
}
void about2() {

    dual x = 2.0;                                 // the input variable x
    dual u = f(x);                                // the output variable u

    double dudx = derivative(f, wrt(x), at(x));   // evaluate the derivative du/dx

    std::cout << "u = " << u << std::endl;        // print the evaluated output u
    std::cout << "du/dx = " << dudx << std::endl; // print the evaluated derivative du/dx

    std::cout << "Hello, from emscr2! autodiff 2\n";
}
int main() {
    std::cout << "main.\n";
}
