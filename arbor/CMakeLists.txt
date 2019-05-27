cmake_minimum_required(VERSION 3.9)
project(arbor-mitral LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 14)

find_package(arbor REQUIRED)
add_executable(soma single-soma.cpp)

target_link_libraries(soma PRIVATE arbor::arbor arbor::arborenv)
target_include_directories(soma PRIVATE nlohmann)

set_target_properties(soma PROPERTIES OUTPUT_NAME soma)

