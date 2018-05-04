cmake_minimum_required(VERSION 3.3)
project(pyhepmc)

SET(CMAKE_CXX_STANDARD 14)

SET(PYTHON_EXECUTABLE python2 CACHE STRING "use Python-2.7")

include_directories(hepmc3/include)
file(GLOB hepmc3_sources hepmc3/src/*.cc hepmc3/src/Search/*.cc)
add_library(hepmc3 STATIC ${hepmc3_sources})
set_target_properties(hepmc3 PROPERTIES CXX_VISIBILITY_PRESET hidden)

add_subdirectory(pybind11)
pybind11_add_module(pyhepmc main.cpp)
target_link_libraries(pyhepmc PRIVATE hepmc3)