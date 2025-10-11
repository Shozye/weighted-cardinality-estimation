cmake_minimum_required(VERSION 3.18)
project(weighted_cardinality_estimation LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_CLANG_TIDY "clang-tidy;-warnings-as-errors=*")

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

add_library(murmurhash3 OBJECT
    lib/murmurhash3/MurmurHash3.cpp
)

set_target_properties(murmurhash3 PROPERTIES CXX_CLANG_TIDY "")

target_include_directories(murmurhash3 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/murmurhash3
)

add_library(compact_vector INTERFACE)
target_include_directories(compact_vector SYSTEM INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/compact_vector
)

pybind11_add_module(_core
    src/weighted_cardinality_estimation/_bindings.cpp
    src/weighted_cardinality_estimation/exp_sketch.cpp
    src/weighted_cardinality_estimation/fast_exp_sketch.cpp
    src/weighted_cardinality_estimation/base_q_sketch.cpp
    src/weighted_cardinality_estimation/fast_q_sketch.cpp
    src/weighted_cardinality_estimation/q_sketch_dyn.cpp
    src/weighted_cardinality_estimation/utils.cpp
    src/weighted_cardinality_estimation/seeds.cpp
    src/weighted_cardinality_estimation/fisher_yates.cpp
    $<TARGET_OBJECTS:murmurhash3>
)

install(TARGETS _core
        LIBRARY DESTINATION weighted_cardinality_estimation)

target_link_libraries(_core PRIVATE compact_vector)

target_include_directories(_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/murmurhash3
    ${CMAKE_CURRENT_SOURCE_DIR}/src/weighted_cardinality_estimation
)


target_compile_features(_core PRIVATE cxx_std_17)
target_compile_options(_core       PRIVATE -Wall -Wextra -Wpedantic -Wreorder)
