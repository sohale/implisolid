cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_C_COMPILER "clang-3.8")
set(CMAKE_CXX_COMPILER "clang++-3.8")
add_compile_options(-stdlib=libc++)
add_compile_options(-std=c++14)
project(pymplicit)
add_subdirectory(pybind11)
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})


INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )


pybind11_add_module(pymplicit pyInterface.cpp)
