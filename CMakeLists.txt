cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(Fast-SLMC CXX)
set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_executable(SLMC-demo test/splayTest.cpp)
