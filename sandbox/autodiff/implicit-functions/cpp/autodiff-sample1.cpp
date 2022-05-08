// basaed on emscr1.cpp
#include <iostream>

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
    std::cout << "Hello, from emscr2! autodiff 2\n";
}
int main() {
    std::cout << "Hello, from emscr2! autodiff 3\n";
}
