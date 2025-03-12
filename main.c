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

    int i,j,k;

    // Mesh file reading
    char *mf_test1 = "Mesh_files/car1x1t_4";
    int type, n_nod, n_elem, n_nod_elem, n_edg_elem;
    float **nod_coords;
    int **nod_gNb;
    int **ref_edg;
    read_mesh(mf_test1,&type,&n_nod,&nod_coords,&n_elem,&nod_gNb,&n_nod_elem,&n_edg_elem,&ref_edg);

    // Edges conditions assignment
    int ref_interior = 0;
    int ref_Dh[1] = {1};   int n_Dh = 1;
    int ref_Dnh[1] = {4};  int n_Dnh = 1;
    int ref_NF[2] = {2,3}; int n_NF = 2;

    float **a_K = matF_alloc(n_nod_elem,DIM);

    // Elementary results memory allocation (they will be correctly initialize for every element K
    float **A_K = matF_alloc(n_nod_elem,n_nod_elem);
    float *l_K = malloc(n_nod_elem*sizeof(float));
    int *nodes_D = malloc(n_nod_elem*sizeof(int));
    float *uD_aK = malloc(n_nod_elem*sizeof(float));

    for (k=0; k<n_elem; k++){

        // Initialize elementary result and Dirichlet containers with correct values
        for(i=0; i<n_nod_elem; i++){
            l_K[i] = 0;
            nodes_D[i] = 1;
            uD_aK[i]= 0; // (3) maybe useless because element with 0 values won't be used / accessed ?
            for(j=0; j<n_nod_elem; j++){
                A_K[i][j] = 0;
            }
        }

        // Get nodes coordinates of element K (a_K)
        //(1) Is it ok to not make a deep copy ? yes bc constant ?
        //      maybe better to make a deep copy to have everything aligned ?
        // for (j=0; j<n_nod_elem; j++) a_K[j] = nod_coords[nod_gNb[i][j]-1]; // -1 because global number begin at 1
        for (i=0; i<n_nod_elem; i++){
            for(j=0; j<DIM; j++){
                a_K[i][j] = nod_coords[nod_gNb[k][i]-1][j]; // -1 because global number begin at 1
                //printf("a_K(i,j) = %f\n",a_K[i][j]);
            }
        }
        // (3) I have to put &() to get the line of i of ref_edg[i] ??
        //for(i=0;i<n_edg_elem;i++) printf("ref_edg_%d of elem %d : %d\n",i,k,ref_edg[k][i]);
        eval_K(ref_interior,ref_Dh,ref_Dnh,ref_NF,n_Dh,n_Dnh,n_NF,type,n_nod_elem,(const float**)a_K,n_edg_elem,ref_edg[k],A_K,l_K,nodes_D,uD_aK);

        imp_eval_K(k+1,type,n_nod_elem,A_K,l_K,nodes_D,uD_aK);

        // (4) How to handle case where a nod is in a corner and is of 2 different edge conditions at the time ??
        //      -> because for now if I'm not wrong the old condition is replaced by the last condition in nodes_D
        //       (*) Is it ok with the hypothesis of continuous conditions ?


    }

    // free memory
    free_mat(a_K);
    free_mat(A_K);
    free(l_K);
    free(nodes_D);
    free(uD_aK);

    return 0;
}
