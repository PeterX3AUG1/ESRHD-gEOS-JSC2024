cmake_minimum_required(VERSION 3.13)
#project(Euler CXX C Fortran)
project(Euler CXX C)

set(CMAKE_C_STANDARD 99)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fms-extensions -funsigned-char -pedantic -W -Wall \
#                               -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings \
#                               -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations \
#                               -Wnested-externs -Wunused -Wfloat-equal -fopenmp")

include_directories(include)
include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${src}
)

file(GLOB_RECURSE src_files src/*.cpp)

add_executable(main
        ${src_files})

# openMP 配置
 FIND_PACKAGE(OpenMP REQUIRED)
 if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 endif ()