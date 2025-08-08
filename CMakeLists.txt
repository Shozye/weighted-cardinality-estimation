cmake_minimum_required(VERSION 3.18)
project(weighted_cardinality_estimation LANGUAGES CXX)

find_package(pybind11 CONFIG REQUIRED)

add_library(murmurhash3 OBJECT
    lib/murmurhash3/MurmurHash3.cpp
)

target_include_directories(murmurhash3 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/murmurhash3
)

pybind11_add_module(_core
    src/weighted_cardinality_estimation/_bindings.cpp
    src/weighted_cardinality_estimation/exp_sketch.cpp
    src/weighted_cardinality_estimation/fast_exp_sketch.cpp
    $<TARGET_OBJECTS:murmurhash3>
)

install(TARGETS _core
        LIBRARY DESTINATION weighted_cardinality_estimation)

target_include_directories(_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/murmurhash3
    ${CMAKE_CURRENT_SOURCE_DIR}/src/weighted_cardinality_estimation
)


target_compile_features(_core PRIVATE cxx_std_17)
