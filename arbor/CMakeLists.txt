cmake_minimum_required(VERSION 3.9)
project(arbor-mitral LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 14)

find_package(arbor REQUIRED)
add_executable(single single.cpp)

target_link_libraries(single PRIVATE arbor::arbor arbor::arborenv)
target_include_directories(single PRIVATE common/cpp/include)

set_target_properties(single PROPERTIES OUTPUT_NAME single)

