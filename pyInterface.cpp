#include <pybind11/pybind11.h>
#include "js_iteration_1/mcc2.cpp"
#include <string>
#include <array>

void*  getVerts(){
    float* verts = get_v_ptr();
    int arraysize;
    arraysize = get_v_size();
    std::array<float, arraysize> verts_11;
    verts_11.data = verts;
    return verts_11;

}

std::array<float, array_size>  getFaces(){
    float* faces = get_f_ptr();
    int array_size = get_f_size();
    std::array<float, array_size> faces_11;
    faces_11.data = faces;
    return faces_11;
}

void buildGeometry(std::string shape_parameters_json, std::string mc_parameters_json){
    char* shapeParams = shape_parameters_json.data;
    char* mcParams = mc_parameters_json.data;
    build_geometry(shapeParams, mcParams);

}


namespace py = pybind11;

PYBIND11_PLUGIN(pyModeler) {
    py::module m("pyModeler", "interface to implicit modeling functions");

    m.def("get_verts", &getVerts, "get a vertex list");
    m.def("finish_geometry", &finish_geometry, "finish geometry");
    m.def("get_faces", &getFaces, "get a vertex list");
    m.def("build_geometry", &buildGeometry, "get a vertex list");

#ifdef VERSION_INFO
    m.attr("__version__") = py::str(VERSION_INFO);
#else
    m.attr("__version__") = py::str("dev");
#endif

    return m.ptr();
}


