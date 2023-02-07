cmake_minimum_required(VERSION 3.24)
project(solver_library)

set(CMAKE_CXX_STANDARD 11)

enable_testing()

add_subdirectory(src)
add_subdirectory(test)