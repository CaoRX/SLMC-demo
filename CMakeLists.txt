cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(Fast-SLMC CXX)
set (CMAKE_CXX_STANDARD 17)

# set (CMAKE_VERBOSE_MAKEFILE on)

# set (LAPACK_LINK_FLAGS "-lcblas -llapacke -llapack -lblas")
# set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINK_FLAGS}")

# message("cmake_exe_linker_flags ${CMAKE_EXE_LINKER_FLAGS}")

# configure_file(${PROJECT_SOURCE_DIR}/config/QDM-config.cmake.in ${PROJECT_SOURCE_DIR}/config/QDM-config.cmake @ONLY)

include_directories(${PROJECT_SOURCE_DIR}/include)

# file(GLOB_RECURSE SOURCE_FILES
# 	${PROJECT_SOURCE_DIR}/src/*.cpp
# 	${PROJECT_SOURCE_DIR}/include/*.hpp
# )
# message(STATUS $ENV{SOURCE_FILES})
# add_library(SLMCLib STATIC ${SOURCE_FILES})

add_executable(SLMC-demo.out test/splayTest.cpp)
# target_link_libraries(splayTest.out SLMCLib ${lapackblas_libraries})