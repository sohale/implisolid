cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_C_COMPILER "clang-3.8")
set(CMAKE_CXX_COMPILER "clang++-3.8")
project(pymplicit)
add_subdirectory(/usr/include/pybind11 pybind11)



include_directories("/usr/include/boost/boost_1_61_0/")
include_directories("/usr/include/eigen")

pybind11_add_module(pymplicit pyInterface.cpp)

