#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "js_iteration_1/mcc2.cpp"
#include <string>
#include <array>
#include <iostream>

namespace py = pybind11;
using namespace std;

py::array_t<float> getVerts(){
    float* verts = (float *)get_v_ptr();
    unsigned int array_size;
    array_size = get_v_size();


    return py::array(py::buffer_info(verts,
                    sizeof(float),
                    py::format_descriptor<float>::value,
                    2,
                    {array_size, 3},
                    {sizeof(float)*3, sizeof(float)}));

}

py::array_t<int>  getFaces(){

    int* faces = (int *)get_f_ptr();
    unsigned int array_size;
    array_size = get_f_size();

    return py::array(py::buffer_info(faces, sizeof(int),
                   py::format_descriptor<int>::value,
                   2,
                   {array_size, 3},
                   {sizeof(int)*3,sizeof(int) }));
}

void buildGeometry(std::string shape_parameters_json, std::string mc_parameters_json){
    const char* shapeParams = shape_parameters_json.data();
    const char* mcParams = mc_parameters_json.data();
    build_geometry(shapeParams, mcParams);

}


namespace py = pybind11;

PYBIND11_PLUGIN(pymplicit) {
    py::module m("pymplicit", "interface to implicit modeling functions");

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


