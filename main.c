/**
 * @file main.c
 * @brief Here we lunch every unit test
 *
 * @details The mesh files data to lunch test are located in the Mesh_files folder
 */

#include "meshing.h"
#include "unit_tests.h"
#include "elem_eval.h"
#include "tab_mngmt.h"
#include "stdio.h"
#include "stdlib.h"

int main(int argc, char* argv[])
{
    /// Test write_mesh
    /*
    int nrefdom = 0;
    const int nrefcot[] = {1, 2, 3, 4};

    // car1x1q_4 re_write test
    test_write_mesh(0.0,1.0,0.0,1.0,5,5,1,nrefdom,nrefcot,"Mesh_files/check_write_car1x1q_4","Mesh_files/car1x1q_4");

    // car1x1t_1 re_write test
    test_write_mesh(0.0,1.0,0.0,1.0,2,2,2,nrefdom,nrefcot,"Mesh_files/check_write_car1x1t_1","Mesh_files/car1x1t_1");

    // car1x1t_4 re_write test
    test_write_mesh(0.0,1.0,0.0,1.0,5,5,2,nrefdom,nrefcot,"Mesh_files/check_write_car1x1t_4","Mesh_files/car1x1t_4");

    // car3x3t_3 re_write test
    test_write_mesh(0.0,3.0,0.0,3.0,4,4,2,nrefdom,nrefcot,"Mesh_files/check_write_car3x3t_3","Mesh_files/car3x3t_3");
    */

    /// Test read_mesh
    /*
    char* read_mesh_file1 = "Mesh_files/car1x1q_4";
    char* check_read_mesh_output1 = "Mesh_files/check_read_car1x1q_4";
    test_read_mesh(read_mesh_file1,check_read_mesh_output1);

    char* read_mesh_file2 = "Mesh_files/car1x1t_1";
    char* check_read_mesh_output2 = "Mesh_files/check_read_car1x1t_1";
    test_read_mesh(read_mesh_file2,check_read_mesh_output2);

    char* read_mesh_file3 = "Mesh_files/car1x1t_4";
    char* check_read_mesh_output3 = "Mesh_files/check_read_car1x1t_4";
    test_read_mesh(read_mesh_file3,check_read_mesh_output3);

    char* read_mesh_file4 = "Mesh_files/car3x3t_3";
    char* check_read_mesh_output4 = "Mesh_files/check_read_car3x3t_3";
    test_read_mesh(read_mesh_file4,check_read_mesh_output4);

    // If all test are passed, nothing is printed
    */

    /// Test eval_K
    // -> please check expression of problem function (a_ij, bN, ...) in the problem_functions.h file before
    //printf("=============================================================================\n"
    //       "Test of eval_K for the simple case : 1 box that contains 2 triangle element :\n"
    //       "=============================================================================\n");
    //test_eval_K("Mesh_files/check_read_car1x1t_1");

    printf("=====================================================================================\n"
           "Test of eval_K for the simple case : 3x3 boxes that contain 2 triangle element each :\n"
           "=====================================================================================\n");
    test_eval_K("Mesh_files/check_read_car3x3t_3");

    return 0;
}
