cmake_minimum_required(VERSION 3.29)
project(FEM_project C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(FEM_project
        main.c
        meshing.c
        tab_mngmt.c
        unit_tests.c)
