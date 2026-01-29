cmake_minimum_required(VERSION 3.14)
project(Roots VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(roots src/roots/roots.cpp)
target_include_directories(roots PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/roots)
target_compile_features(roots PUBLIC cxx_std_17)

add_executable(roots_demo src/roots/main.cpp)
target_link_libraries(roots_demo PRIVATE roots)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(test_roots tests/test_roots.cpp)
target_link_libraries(test_roots PRIVATE roots GTest::gtest_main)
add_test(NAME roots_test COMMAND test_roots)