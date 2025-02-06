cmake_minimum_required(VERSION 3.29)
project(FEM_project C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

include_directories(.)

add_executable(FEM_project
        main.c
        meshing.c
        tab_mngmt.c
        unit_tests.c
        elem_eval.c)
