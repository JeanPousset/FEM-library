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
#include <stdio.h>
#include <stdlib.h>

int nucas;

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
    /*
    printf("=============================================================================\n"
           "Test of eval_K for the simple case : 1 box that contains 2 triangle element :\n"
           "=============================================================================\n");
    test_eval_K("Mesh_files/check_read_car1x1t_1");
    */
    /*
    printf("=====================================================================================\n"
           "Test of eval_K for basic case : 3x3 boxes that contain 2 triangle element each :\n"
           "=====================================================================================\n");
    test_eval_K("Mesh_files/check_read_car3x3t_3");
    */

    /* // There isn't any correction
    printf("=====================================================================================\n"
           "Test of eval_K for basic case : 4x4 boxes that contain 1 quadrangle element each :\n"
           "=====================================================================================\n");
    test_eval_K("Mesh_files/check_read_car1x1q_4");


    printf("=====================================================================================\n"
           "Test of eval_K for basic case : 4x4 boxes that contain 2 triangle element each :\n"
           "=====================================================================================\n");
    test_eval_K("Mesh_files/check_read_car1x1t_4");
    */


    /*
    printf("=====================================================================================\n"
           "Test of eval_K for case exemple results d1t1_2 :\n"
           "=====================================================================================\n");
    test_eval_K("Mesh_files/d1t1_2");
    */
    nucas = 1;


    /// Assembly test
    //test_assembly("Mesh_files/d1t1_2");
    //test_assembly("Mesh_files/d1t1_4");
    //test_assembly("Mesh_files/d1q1_4");
    //test_assembly("Mesh_files/check_read_car1x1t_1");
    //test_assembly("Mesh_files/check_read_car3x3t_3");

    /// DMS -> OMS with BC condition into 2nd member test
    //test_DMS_to_OMS_with_BC("Mesh_files/d1t1_2");
    //test_DMS_to_OMS_with_BC("Mesh_files/d1t1_4");
    //test_DMS_to_OMS_with_BC("Mesh_files/d1q1_4");
    //test_DMS_to_OMS_with_BC("Mesh_files/car1x1t_1");
    //test_DMS_to_OMS_with_BC("Mesh_files/car3x3t_3");

    /// OMS -> profile
    //test_OMS_to_profile("Mesh_files/d1t1_2");
    //test_OMS_to_profile("Mesh_files/d1t1_4");
    //test_OMS_to_profile("Mesh_files/d1q1_4");
    //test_OMS_to_profile("Mesh_files/car1x1t_1");
    //test_OMS_to_profile("Mesh_files/car3x3t_3");

    /// Results



    nucas = 1;
    //test_solver("Mesh_files/d1t1_2");
    //test_solver("Mesh_files/d1t1_4");
    //test_solver("Mesh_files/d1q1_4");
    //test_solver("Mesh_files/car1x1t_4");


    // nucas = 1;
    //Poisson_Dh_ex1("Mesh_files/d1t1_4",0);


    /// Convergence test for Poisson with homogeneous Dirichlet boundaries
    /*
    nucas = 1; int export11 = 11;
    Poisson_Dh("Mesh_files/d1t1_2",export11);
    Poisson_Dh("Mesh_files/d1t1_4",export11);
    Poisson_Dh("Mesh_files/d1t1_8",export11);
    Poisson_Dh("Mesh_files/d1t1_16",export11);
    Poisson_Dh("Mesh_files/d1t1_32",export11);
    Poisson_Dh("Mesh_files/d1t1_64",export11);
    */

    /*
    nucas = 1; int export12 = 12;
    Poisson_Dh("Mesh_files/d1q1_2",export12);
    Poisson_Dh("Mesh_files/d1q1_4",export12);
    Poisson_Dh("Mesh_files/d1q1_8",export12);
    Poisson_Dh("Mesh_files/d1q1_16",export12);
    Poisson_Dh("Mesh_files/d1q1_32",export12);
    Poisson_Dh("Mesh_files/d1q1_64",export12);
    */

    /*
    nucas = 2; int export21 = 21;
    Poisson_Dh("Mesh_files/d1t1_2",export21);
    Poisson_Dh("Mesh_files/d1t1_4",export21);
    Poisson_Dh("Mesh_files/d1t1_8",export21);
    Poisson_Dh("Mesh_files/d1t1_16",export21);
    Poisson_Dh("Mesh_files/d1t1_32",export21);
    Poisson_Dh("Mesh_files/d1t1_64",export21);
    */

    /*
    nucas = 2; int export22 = 22;
    Poisson_Dh("Mesh_files/d1q1_2",export22);
    Poisson_Dh("Mesh_files/d1q1_4",export22);
    Poisson_Dh("Mesh_files/d1q1_8",export22);
    Poisson_Dh("Mesh_files/d1q1_16",export22);
    Poisson_Dh("Mesh_files/d1q1_32",export22);
    Poisson_Dh("Mesh_files/d1q1_64",export22);
    */

    /// Convergence test for -Laplacien(u) + u = f_Omega with Neumann boundaries

    /*
    nucas = 3; int export31 = 31;
    elliptic_N("Mesh_files/d1t1_2",export31);
    elliptic_N("Mesh_files/d1t1_4",export31);
    elliptic_N("Mesh_files/d1t1_8",export31);
    elliptic_N("Mesh_files/d1t1_16",export31);
    elliptic_N("Mesh_files/d1t1_32",export31);
    elliptic_N("Mesh_files/d1t1_64",export31);
    */

    /*
    nucas = 3; int export32 = 32;
    elliptic_N("Mesh_files/d1q1_2",export32);
    elliptic_N("Mesh_files/d1q1_4",export32);
    elliptic_N("Mesh_files/d1q1_8",export32);
    elliptic_N("Mesh_files/d1q1_16",export32);
    elliptic_N("Mesh_files/d1q1_32",export32);
    elliptic_N("Mesh_files/d1q1_64",export32);
    */

    /// ------------------------- DOMAIN 2  --------------------------------

    /*
    nucas = 1; int export211 = 211;
    Poisson_Dnh("Mesh_files/d2t1_2",export211);
    Poisson_Dnh("Mesh_files/d2t1_4",export211);
    Poisson_Dnh("Mesh_files/d2t1_8",export211);
    Poisson_Dnh("Mesh_files/d2t1_16",export211);
    Poisson_Dnh("Mesh_files/d2t1_32",export211);
    Poisson_Dnh("Mesh_files/d2t1_64",export211);
    */

    /*
    nucas = 1; int export212 = 212;
    Poisson_Dnh("Mesh_files/d2q1_2",export212);
    Poisson_Dnh("Mesh_files/d2q1_4",export212);
    Poisson_Dnh("Mesh_files/d2q1_8",export212);
    Poisson_Dnh("Mesh_files/d2q1_16",export212);
    Poisson_Dnh("Mesh_files/d2q1_32",export212);
    Poisson_Dnh("Mesh_files/d2q1_64",export212);
    */

    /*
    nucas = 2; int export221 = 221;
    Poisson_Dnh("Mesh_files/d2t1_2",export221);
    Poisson_Dnh("Mesh_files/d2t1_4",export221);
    Poisson_Dnh("Mesh_files/d2t1_8",export221);
    Poisson_Dnh("Mesh_files/d2t1_16",export221);
    Poisson_Dnh("Mesh_files/d2t1_32",export221);
    Poisson_Dnh("Mesh_files/d2t1_64",export221);


    nucas = 2; int export222 = 222;
    Poisson_Dnh("Mesh_files/d2q1_2",export222);
    Poisson_Dnh("Mesh_files/d2q1_4",export222);
    Poisson_Dnh("Mesh_files/d2q1_8",export222);
    Poisson_Dnh("Mesh_files/d2q1_16",export222);
    Poisson_Dnh("Mesh_files/d2q1_32",export222);
    Poisson_Dnh("Mesh_files/d2q1_64",export222);
    */

    /// Convergence test for -Laplacien(u) + u = f_Omega with Neumann boundaries
    // -> to complexe (false)

    /*
    nucas = 3; int export231 = 231;
    elliptic_NDnh("Mesh_files/d2t1_2",export231);
    elliptic_NDnh("Mesh_files/d2t1_4",export231);
    elliptic_NDnh("Mesh_files/d2t1_8",export231);
    elliptic_NDnh("Mesh_files/d2t1_16",export231);
    elliptic_NDnh("Mesh_files/d2t1_32",export231);
    elliptic_NDnh("Mesh_files/d2t1_64",export231);
    */

    /*
    nucas = 3; int export232 = 232;
    elliptic_NDnh("Mesh_files/d2q1_2",export232);
    elliptic_NDnh("Mesh_files/d2q1_4",export232);
    elliptic_NDnh("Mesh_files/d2q1_8",export232);
    elliptic_NDnh("Mesh_files/d2q1_16",export232);
    elliptic_NDnh("Mesh_files/d2q1_32",export232);
    elliptic_NDnh("Mesh_files/d2q1_64",export232);
    */


    return 0;
}
