cmake_minimum_required (VERSION 3.0)
project (Molepy)


set(PYBIND11_PYTHON_VERSION 2.7)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(pybind11)

pybind11_add_module(CartSpher CartSpher.cc)
pybind11_add_module(kinetic kinetic.cc)
pybind11_add_module(os_hrr os_hrr.cc)
pybind11_add_module(overlap overlap.cc)
pybind11_add_module(potential potential.cc)

add_library(boys boys.cc)
add_library(os_vrr os_vrr.cc)

target_link_libraries(os_vrr PRIVATE boys)
target_link_libraries(potential PRIVATE boys)
target_link_libraries(os_hrr PRIVATE os_vrr)


set_target_properties(CartSpher kinetic os_hrr overlap potential
     PROPERTIES
     SUFFIX ".so"
     LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)
