cmake_minimum_required(VERSION 3.24)

project(ellipticalpan_test CXX)
set (CMAKE_CXX_STANDARD 11)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.12.0.zip
)
FetchContent_MakeAvailable(googletest)
enable_testing()

add_executable(EllipticalPANTest test/*.cpp)
target_link_libraries(EllipticalPANTest gtest_main)

include(GoogleTest)
gtest_discover_tests(EllipticalPANTest)
