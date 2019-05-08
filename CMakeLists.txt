
cmake_minimum_required(VERSION 3.10.2)

set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Base degridding project
project(degridder)
add_executable(degridder main.cpp degridder.c)
target_link_libraries(degridder m)

# Unit testing for degridding
project(tests)
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
add_executable(tests unit_testing.cpp degridder.c)
target_link_libraries(tests ${GTEST_LIBRARIES} pthread)
