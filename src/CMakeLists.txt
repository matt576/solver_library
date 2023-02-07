set(LIBRARY_NAME ode)

# header files
set(MODULE_HEADERS 
    vector.h
    function.h
    solver.h)

# cpp files
set(MODULE_SOURCES
    vector.cpp
    solver.cpp)

# build the module library
add_library(${LIBRARY_NAME} ${MODULE_HEADERS} ${MODULE_SOURCES})
add_library(lib::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

target_include_directories(
    ${LIBRARY_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_executable(${LIBRARY_NAME}_vector vector.cpp)
add_executable(${LIBRARY_NAME}_solver solver.cpp)

# require C++11
target_compile_features(${LIBRARY_NAME} PUBLIC cxx_std_11)

set_target_properties(${LIBRARY_NAME} PROPERTIES LINKER_LANGUAGE CXX)
