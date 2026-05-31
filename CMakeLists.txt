cmake_minimum_required(VERSION 3.18)
project(weighted_cardinality_estimation LANGUAGES CXX)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Only enable clang-tidy for Debug builds, where it's expected to be installed.
# cibuildwheel performs Release builds, so this block will be skipped.
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    find_program(CLANG_TIDY_EXE clang-tidy)
    if(CLANG_TIDY_EXE)
        message(STATUS "Enabling clang-tidy for Debug build.")
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};-warnings-as-errors=*")
    else()
        message(WARNING "clang-tidy not found, CXX linter disabled.")
    endif()
endif()

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

add_library(murmurhash3 OBJECT
    lib/murmurhash3/MurmurHash3.cpp
)

set_target_properties(murmurhash3 PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_target_properties(murmurhash3 PROPERTIES CXX_CLANG_TIDY "")

target_include_directories(murmurhash3 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/murmurhash3
)

add_library(compact_vector INTERFACE)
target_include_directories(compact_vector SYSTEM INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/compact_vector
)
add_library(pcg_random INTERFACE)
target_include_directories(pcg_random SYSTEM INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/pcg_random
)

file(GLOB WCE_SOURCES
     CONFIGURE_DEPENDS
     "src/weighted_cardinality_estimation/*.cpp"
)

pybind11_add_module(_core
    ${WCE_SOURCES}
    $<TARGET_OBJECTS:murmurhash3>
)

install(TARGETS _core
        LIBRARY DESTINATION weighted_cardinality_estimation)

target_link_libraries(_core PRIVATE compact_vector)
target_link_libraries(_core PRIVATE pcg_random)

target_include_directories(_core PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/murmurhash3
    ${CMAKE_CURRENT_SOURCE_DIR}/src/weighted_cardinality_estimation
)

target_compile_features(_core PRIVATE cxx_std_17)
target_compile_options(_core PRIVATE -Wall -Wextra -Wpedantic -Wreorder)
