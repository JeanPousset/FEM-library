/**
 * @file unit_tests.c
 * @brief test-functions for unit tests
 * @details Contains test-functions that check if functions are correctly working
 */

#include "unit_tests.h"
#include "meshing.h"
#include "tab_mngmt.h"
#include "elem_eval.h"
#include "assembly.h"
#include "forfun.h"
#include <stdio.h>
#include <stdlib.h>


/// Check if a 2x2 matrix is correctly inverted
int test_inv2x2(void) {
    int i, j, res = 1;
    float **M = matF_alloc(2, 2);
    float **M_inv = matF_alloc(2, 2);
    float **M_cor = matF_alloc(2, 2);
    float det;
    M[0][0] = 3;
    M[0][1] = 2;
    M[1][0] = 2;
    M[1][1] = 2;
    M_cor[0][0] = 1;
    M_cor[0][1] = -1;
    M_cor[1][0] = -1;
    M_cor[1][1] = 1.5;

    det = inv_2x2(M, M_inv);
    if (det != 2.0) res = 0;     // Checking if the matrix is correctly inverted
    for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) if (M_cor[i][j] != M_inv[i][j]) res = 0;
    if (res == 1) { printf(" 'inv2x2' function passed test with success \n"); }

    // memory clean
    free_mat(M);
    free_mat(M_inv);
    free_mat(M_cor);
    return res;
}


/// Function that lunch a diff command on the terminal to compare 2 file contents.
void diff_file(const char *f1, const char *f2) {
    char diff_command[512];
    // The file compare command is different in Windows OS
#ifdef _WIN32
    snprintf(diff_command, sizeof(diff_command), "FC %s %s", f1, f2);
#else
    snprintf(diff_command, sizeof(diff_command), "diff %s %s", f1, f2);
#endif
    FILE *fp = popen(diff_command, "r");
    if (fp == NULL) perror("Error when opening a pipe");

    char answer_buffer[256];
    while (fgets(answer_buffer, sizeof(answer_buffer), fp) != NULL) printf("%s", answer_buffer);
    pclose(fp);
}


///    Test the write_mesh function implemented in "meshing.c"
void test_write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int nrefdom, const int *nrefcot,
                     char *mesh_file_path, char *correct_mesh_file) {
    write_mesh(a, b, c, d, n1, n2, t, nrefdom, nrefcot, mesh_file_path);
    diff_file(mesh_file_path, correct_mesh_file);
}

/// Test the read_mesh form "meshing.c/h", read an ouput mesh file and rewrite its content into the output mesh file
void test_read_mesh(char *input_mesh, char *out_check_mesh) {
    int type, n_nodes, n_elem, n_nod_elem, n_edges;
    float **p_coords;
    int **p_gNb_node;
    int **p_refEdg;

    // reading input_meshfile
    read_mesh(input_mesh, &type, &n_nodes, &p_coords, &n_elem, &p_gNb_node, &n_nod_elem, &n_edges, &p_refEdg);
    // Opening ouput meshfile in writting mode
    FILE *mesh_f = fopen(out_check_mesh, "w");    // mode 'w' for writing (delete everything there was before)
    if (mesh_f == (FILE *) NULL) {
        printf("Error when trying to open the output (mesh) file : '%s' in test_read_mesh function", out_check_mesh);
        perror("^^ IO - Error ^^");
    }

    // Writing output (test) meshfile
    int i, j;
    fprintf(mesh_f, "%d", n_nodes);
    // coordinates
    for (i = 0; i < n_nodes; i++) fprintf(mesh_f, "\n%f %f", p_coords[i][0], p_coords[i][1]);
    // nodes global number and edges references number

    fprintf(mesh_f, "\n%d %d %d %d", n_elem, type, n_nod_elem, n_edges);
    for (i = 0; i < n_elem; i++) {
        fprintf(mesh_f, "\n");
        for (j = 0; j < n_nod_elem; j++) fprintf(mesh_f, "%d ", p_gNb_node[i][j]);
        for (j = 0; j < n_edges; j++) fprintf(mesh_f, "%d ", p_refEdg[i][j]);
    }
    fprintf(mesh_f, "\n"); // One list new line to have exactly the same shape as other mesh files
    fclose(mesh_f);

    // clean tab memory
    free_mat(p_coords);
    free_mat(p_gNb_node);
    free_mat(p_refEdg);

    // Using diff command line to test if there is a difference between the original and the file
    diff_file(input_mesh, out_check_mesh);
}


/// Print results of A_K, l_K (2nd member), Dirichlet conditions in each node and values of non-homogeneous Dirichlet nodes
void print_eval_K(int K, int type, int n_nod_elem, float **A_K, float *l_K, int *nodes_D, float *uD_aK) {
    int i, j;
    printf("\n");
    printf(" ELEMENT=%3d    DE TYPE=%5d    NB NOEUDS=%2d\n", K, type, n_nod_elem);
    printf(" NuDElem   uDElem    SMbrElem    MatElem\n");
    for (i = 0; i < n_nod_elem; i++) {
        printf(" %6d %10.4e %10.4e", nodes_D[i], uD_aK[i], l_K[i]);
        for (j = 0; j <= i; j++) { printf(" %10.4e", A_K[i][j]); }
        printf("\n");
    }
}

/// Test the eval_K function on every mesh element and call print_eval_K so we can compare to the exact solution
void test_eval_K(const char *mesh_file) {
    int i, j, k;

    // Mesh file reading
    int type, n_nod, n_elem, n_nod_elem, n_edg_elem;
    float **nod_coords;
    int **nod_gNb;
    int **ref_edg;
    read_mesh(mesh_file, &type, &n_nod, &nod_coords, &n_elem, &nod_gNb, &n_nod_elem, &n_edg_elem, &ref_edg);

    // Edges conditions assignment
    int ref_interior = 0;
    int ref_Dh[1] = {1};
    int n_Dh = 1;
    int ref_Dnh[1] = {4};
    int n_Dnh = 1;
    int ref_NF[2] = {2, 3};
    int n_NF = 2;

    // Elementary results memory allocation (they will be correctly initialize for every element K
    float **A_K = matF_alloc(n_nod_elem, n_nod_elem);
    float *l_K = malloc(n_nod_elem * sizeof(float));
    int *nodes_D = malloc(n_nod_elem * sizeof(int));
    float *uD_aK = malloc(n_nod_elem * sizeof(float));
    float **a_K = malloc(n_nod_elem * sizeof(float *)); // Special case : dynamic tab of float pointers

    for (k = 0; k < n_elem; k++) {

        for (i = 0; i < n_nod_elem; i++) {
            // Initialize elementary result and Dirichlet containers with correct values
            l_K[i] = 0;
            nodes_D[i] = 1;
            uD_aK[i] = 0; // (JP:) maybe useless because element with 0 values won't be used / accessed ?
            for (j = 0; j < n_nod_elem; j++) {
                A_K[i][j] = 0;
            }
            // get nodes coordinates of K (we could also use selectPts)
            a_K[i] = nod_coords[nod_gNb[k][i] - 1]; /// @warning : only a simple pointers copy (not a deep copy)
        }

        eval_K(ref_interior, ref_Dh, ref_Dnh, ref_NF, n_Dh, n_Dnh, n_NF, type, n_nod_elem, (const float **) a_K,
               n_edg_elem, ref_edg[k], A_K, l_K, nodes_D, uD_aK);

        print_eval_K(k + 1, type, n_nod_elem, A_K, l_K, nodes_D, uD_aK);

        // (4) How to handle case where a nod is in a corner and is of 2 different edge conditions at the time ??
        //      -> because for now if I'm not wrong the old condition is replaced by the last condition in nodes_D
        //       (*) It's ok with the hypothesis of continuous conditions for Dirichlet conditions

    }
    // free memory
    free(a_K);      // Only an array of pointers -> not really a matrix
    free_mat(A_K);
    free(l_K);
    free(nodes_D);
    free(uD_aK);
}


void test_assembly(const char *mesh_file) {

    // Mesh file reading
    int type, n_nod, n_elem, n_nod_elem, n_edg_elem;
    float **nod_coords;
    int **nod_gNb;
    int **ref_edg;
    read_mesh(mesh_file, &type, &n_nod, &nod_coords, &n_elem, &nod_gNb, &n_nod_elem, &n_edg_elem, &ref_edg);

    // Edges conditions assignment
    int ref_interior = 0;
    int ref_Dh[1] = {1};
    int n_Dh = 1;
    int ref_Dnh[1] = {4};
    int n_Dnh = 1;
    int ref_NF[2] = {2, 3};
    int n_NF = 2;

    int n_lowA = NON_0_ROW * n_nod;
    // Second member in the linear system of the variational formulation's discretization
    float *sec_mb = calloc(n_nod, sizeof(float)); // Avant
    // array to know if nodes are on a Dirichlet edge
    int *D_gb_nod = malloc(n_nod * sizeof(int));
    // Contains indices of the first non-zero element of the line i of the bottom triangular part of A_sparse
    float *D_values = malloc(n_nod * sizeof(float));
    int *first_non0_row = calloc(n_nod, sizeof(int));
    float *A_sparse = calloc(n_nod + n_lowA, sizeof(float));
    // column indices of non-zeros element of triangular bottom part of A_sparse
    int *col_ind = malloc(n_lowA * sizeof(float));
    // Position in triangular bottom part of A_sparse of the next element in the same line
    int *nextInRow = calloc(n_lowA, sizeof(float));

    assembly(type, n_nod, n_elem, n_nod_elem, n_edg_elem, (const float **) nod_coords, (const int **) nod_gNb,
             (const int **) ref_edg, ref_interior, ref_Dh, ref_Dnh, ref_NF, n_Dh, n_Dnh, n_NF, n_lowA, sec_mb, D_gb_nod,
             D_values, first_non0_row, A_sparse, col_ind, nextInRow);
    // display
    affsmd_(&n_nod, first_non0_row, col_ind, nextInRow, A_sparse, sec_mb, D_gb_nod, D_values);


    // Clean memory
    free_mat(nod_coords);
    free_mat(nod_gNb);
    free_mat(ref_edg);
    free(sec_mb);
    free(D_gb_nod);
    free(D_values);
    free(first_non0_row);
    free(A_sparse);
    free(col_ind);
    free(nextInRow);
}

///@brief Test assembly then conversion DMS -> OMS + BC taken into account
void test_DMS_to_OMS_with_BC(const char *mesh_file) {

    // Mesh file reading
    int type, n_nod, n_elem, n_nod_elem, n_edg_elem;
    float **nod_coords;
    int **nod_gNb;
    int **ref_edg;
    read_mesh(mesh_file, &type, &n_nod, &nod_coords, &n_elem, &nod_gNb, &n_nod_elem, &n_edg_elem, &ref_edg);

    // For conversion DMS to OMS
    int length_A_sparse, i, n_row;

    // Edges conditions assignment
    int ref_interior = 0;
    int ref_Dh[1] = {1};
    int n_Dh = 1;
    int ref_Dnh[1] = {4};
    int n_Dnh = 1;
    int ref_NF[2] = {2, 3};
    int n_NF = 2;

    // * --------- Memory allocation for DMS of A --------- *
    int n_lowA = NON_0_ROW * n_nod;
    // Second member in the linear system of the variational formulation's discretization
    float *sec_mb = calloc(n_nod, sizeof(float));
    // array to know if nodes are on a Dirichlet edge
    int *BC_nod = malloc(n_nod * sizeof(int));
    // Contains indices of the first non-zero element of the line i of the bottom triangular part of A_sparse
    float *D_values = malloc(n_nod * sizeof(float));
    int *first_non0_row_DMS = calloc(n_nod, sizeof(int));
    float *A_DMS = calloc(n_nod + n_lowA, sizeof(float));
    // column indices of non-zeros element of triangular bottom part of A_sparse
    int *col_ind_DMS = malloc(n_lowA * sizeof(float));
    // Position in triangular bottom part of A_sparse of the next element in the same line
    int *nextInRow = calloc(n_lowA, sizeof(float));

    length_A_sparse = assembly(type, n_nod, n_elem, n_nod_elem, n_edg_elem, (const float **) nod_coords,
                               (const int **) nod_gNb,
                               (const int **) ref_edg, ref_interior, ref_Dh, ref_Dnh, ref_NF, n_Dh, n_Dnh, n_NF, n_lowA,
                               sec_mb, BC_nod, D_values, first_non0_row_DMS, A_DMS, col_ind_DMS, nextInRow);


    // * --------- Memory allocation for OMS of A --------- *
    n_row = n_row_with_BC(n_nod,BC_nod);    // Count how many nodes are Dirichlet's
    float *A_OMS = malloc(length_A_sparse * sizeof(float)); // A (sparse) matrix in an OMS
    int *first_non0_row_OMS = calloc(n_nod, sizeof(int));
    int *col_ind_OMS = malloc((length_A_sparse - n_nod) * sizeof(float));
    float *sec_mb_BC = malloc((n_nod - n_row) * sizeof(float)); // 2nd member with BC conditions (smaller than previous)

    // Transforms into OMS and take Neumann nodes into 2nd member
    DMS_to_OMS_with_BC(A_DMS, first_non0_row_DMS, col_ind_DMS, nextInRow, sec_mb, BC_nod, D_values, n_nod, A_OMS,
                       first_non0_row_OMS, col_ind_OMS, sec_mb_BC);

    // Display result
    affsmo_(&n_nod, first_non0_row_OMS, col_ind_OMS, A_OMS, sec_mb_BC);

    // Clean memory for DMS
    free_mat(nod_coords);
    free_mat(nod_gNb);
    free_mat(ref_edg);
    free(sec_mb);
    free(BC_nod);
    free(D_values);
    free(first_non0_row_DMS);
    free(A_DMS);
    free(col_ind_DMS);
    free(nextInRow);

    // Clean memory for OMS
    free(A_OMS);
    free(first_non0_row_OMS);
    free(col_ind_OMS);
    free(sec_mb_BC);
}

int linear_system(const char *mesh_file, float **A_pr, int **profile, float **sec_mb_BC) {
    ///@warning in these function we allocate the result INSIDE !
    // Mesh file reading
    int type, n_nod, n_elem, n_nod_elem, n_edg_elem;
    float **nod_coords;
    int **nod_gNb;
    int **ref_edg;
    read_mesh(mesh_file, &type, &n_nod, &nod_coords, &n_elem, &nod_gNb, &n_nod_elem, &n_edg_elem, &ref_edg);

    // For conversion DMS to OMS
    int length_A_DMS,n_row, i;

    // Edges conditions assignment
    int ref_interior = 0;
    int ref_Dh[1] = {1};
    int n_Dh = 1;
    int ref_Dnh[1] = {4};
    int n_Dnh = 1;
    int ref_NF[2] = {2, 3};
    int n_NF = 2;

    // * --------- Assembly into a DMS for A --------- *

    int n_lowA = NON_0_ROW * n_nod;
    // Second member in the linear system of the variational formulation's discretization
    float *sec_mb = calloc(n_nod, sizeof(float));
    // array to know if nodes are on a Dirichlet edge
    int *BC_nod = malloc(n_nod * sizeof(int));
    // Contains indices of the first non-zero element of the line i of the bottom triangular part of A_sparse
    float *D_values = malloc(n_nod * sizeof(float));
    int *first_non0_row_DMS = calloc(n_nod, sizeof(int));
    float *A_DMS = calloc(n_nod + n_lowA, sizeof(float));
    // column indices of non-zeros element of triangular bottom part of A_sparse
    int *col_ind_DMS = malloc(n_lowA * sizeof(float));
    // Position in triangular bottom part of A_sparse of the next element in the same line
    int *nextInRow = calloc(n_lowA, sizeof(float));

    length_A_DMS = assembly(type, n_nod, n_elem, n_nod_elem, n_edg_elem, (const float **) nod_coords,
                     (const int **) nod_gNb,
                     (const int **) ref_edg, ref_interior, ref_Dh, ref_Dnh, ref_NF, n_Dh, n_Dnh, n_NF, n_lowA,
                     sec_mb, BC_nod, D_values, first_non0_row_DMS, A_DMS, col_ind_DMS, nextInRow);


    // * --------- DMS to OMS with BC inside 2nd member --------- *

    // Memory allocation
    n_row = n_row_with_BC(n_nod,BC_nod);    // Count how many nodes are Dirichlet's
    /// @note OMS size is too big, we'll get the minimum size after transforming it into profile storage
    float *A_OMS = malloc(length_A_DMS * sizeof(float));
    int *first_non0_row_OMS = calloc(n_row, sizeof(int));
    int *col_ind_OMS = malloc(length_A_DMS * sizeof(float)); // same here, too big
    *sec_mb_BC = malloc(n_row * sizeof(float)); // 2nd member with BC conditions (smaller than previous)

    // Transforms into OMS and take Neumann nodes into 2nd member
    DMS_to_OMS_with_BC(A_DMS, first_non0_row_DMS, col_ind_DMS, nextInRow, sec_mb, BC_nod, D_values, n_nod, A_OMS,
                       first_non0_row_OMS, col_ind_OMS, sec_mb_BC);

    // Clean memory for DMS
    free_mat(nod_coords);
    free_mat(nod_gNb);
    free_mat(ref_edg);
    free(sec_mb);
    free(BC_nod);
    free(D_values);
    free(first_non0_row_DMS);
    free(A_DMS);
    free(col_ind_DMS);
    free(nextInRow);


    // * --------- Transforms OMS into profile --------- *

    // Memory allocation for profile storage (here we know the exact size)
    int n_lowA_pr = length_profile(n_row, first_non0_row_OMS, col_ind_OMS);
    *A_pr = calloc(n_lowA_pr + n_row, sizeof(float)); // Must be initialized to 0
    *profile = malloc(n_row * sizeof(int));

    // Transformation into profile
    OMS_to_profile(A_OMS, n_row, first_non0_row_OMS, col_ind_OMS, *A_pr, *profile);

    // Clean memory for OMS
    free(A_OMS);
    free(first_non0_row_OMS);
    free(col_ind_OMS);

    return n_row;
}

///@brief Gets the profile storage of A and prints it
void test_OMS_to_profile(const char *mesh_file) {
    int param_imp = 0;
    int n_row;
    float *A_pr;
    int *profile;
    float *sec_mb_BC;
    n_row = linear_system(mesh_file, &A_pr, &profile, &sec_mb_BC);

    // Print matrix in a profile storage
    impmpr_(&param_imp,&n_row, profile, A_pr, &A_pr[n_row]);

    // Clean memory for A_pr and 2nd member that were allocated in linear_system
    free(A_pr);
    free(profile);
    free(sec_mb_BC);
}


