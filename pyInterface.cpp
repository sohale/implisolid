
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "js_iteration_1/mcc2.cpp"
#include <string>
#include <array>
#include <iostream>
#include <fstream>

namespace py = pybind11;
using namespace std;

py::array_t<REAL> getVerts(){

    ofstream myfile;
    myfile.open ("solimod.log", std::ios_base::app);

    REAL* verts = (REAL *)get_v_ptr();
    unsigned int array_size;
    array_size = get_v_size();

    myfile.close();


    return py::array(py::buffer_info(verts,
                    sizeof(REAL),
                    py::format_descriptor<REAL>::value,
                    2,
                    {array_size, 3},
                    {sizeof(REAL)*3, sizeof(REAL)}));

}

// tomorrow
py::array_t<int>  getFaces(){

    ofstream myfile;
    myfile.open ("solimod.log", std::ios_base::app);

    vertexindex_type* faces = (vertexindex_type *)get_f_ptr();
    unsigned int array_size;
    array_size = get_f_size();

    myfile.close();

    return py::array(py::buffer_info(faces, 
                   sizeof(vertexindex_type),
                   py::format_descriptor<vertexindex_type>::value,
                   2,
                   {array_size, 3},
                   {sizeof(vertexindex_type)*3,sizeof(vertexindex_type) }));
}

void buildGeometry(std::string shape_parameters_json, std::string mc_parameters_json){

    ofstream myfile;
    myfile.open ("solimod.log", std::ios_base::app);

    const char* shapeParams = shape_parameters_json.data();
    const char* mcParams = mc_parameters_json.data();
    build_geometry(shapeParams, mcParams);

    myfile.close();

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


