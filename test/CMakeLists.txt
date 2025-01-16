enable_testing()
include(CTest)

add_executable(test_potential test_potential.cpp)
target_link_libraries(test_potential PRIVATE pot Catch2::Catch2WithMain)

add_executable(test_parallel test_parallel.cpp)
target_link_libraries(test_parallel PRIVATE pot Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(test_potential)
