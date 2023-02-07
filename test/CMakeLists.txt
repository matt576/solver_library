include(CTest)

# helper function to add a test based on a file called TEST_NAME.cpp
function(addTest TEST_NAME)
    set(${TEST_NAME} ${TEST_NAME})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_compile_features(${TEST_NAME} PUBLIC cxx_std_11)
    add_test(${TEST_NAME} ${TEST_NAME})
    target_link_libraries(${TEST_NAME} PRIVATE lib::ode)
endfunction(addTest)

addTest(test_vector)
