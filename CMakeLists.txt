project(particle_filter)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    )
    
add_library(
    sources
    src/particle_filter.cpp
    src/particle_filter.h
    src/helper_functions.h
    src/json.hpp
    src/map.h
    )

set(GOOGLETEST_ROOT googletest/googletest CACHE STRING "Google Test source root")

include_directories(
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/include
    )

set(GOOGLETEST_SOURCES
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest-all.cc
    ${PROJECT_SOURCE_DIR}/${GOOGLETEST_ROOT}/src/gtest_main.cc
    )

foreach(_source ${GOOGLETEST_SOURCES})
    set_source_files_properties(${_source} PROPERTIES GENERATED 1)
endforeach()

add_library(googletest ${GOOGLETEST_SOURCES})

add_executable(
    unit_tests
    src/test/main_tests.cpp
    src/test/particle_filter_tests.cpp
    )

add_dependencies(unit_tests googletest)

target_link_libraries(
    unit_tests
    googletest
    sources
    )

enable_testing()

add_test(unit ${PROJECT_BINARY_DIR}/unit_tests)



if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


add_executable(particle_filter src/main.cpp)

target_link_libraries(particle_filter sources z ssl uv uWS)
