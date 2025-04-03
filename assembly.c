#include "assembly.h"
#include "forfun.h"
#include "tab_mngmt.h"
#include "elem_eval.h"
#include "unit_tests.h"
#include <stdio.h>
#include <stdlib.h>

void assembly(
        int type,
        int n_nod,
        int n_elem,
        int n_nod_elem,
        int n_edg_elem,
        const float **nod_coords,
        const int **nod_gNb,
        const int **ref_edg,
        int ref_interior,      /// @param ref_interior reference number for edges inside the domain
        const int *ref_Dh,     /// @param ref_Dh array that contains number of domain boundaries carrying an homogeneous Dirichlet condition
        const int *ref_Dnh,    ///
        const int *ref_NF,     ///
        int n_Dh,              ///
        int n_Dnh,             ///
        int n_NF,               ///
        int n_lowA,
        float *sec_mb,       /// @warning : have to be initialized to 0
        int *D_gb_nod,
        float *D_values,
        int *first_non0_row,
        float *A_sparse,     /// @warning : have to be initialized to 0
        int *col_ind,
        int *nextInRow
) {
    // Load mesh
    int k, i, j;   // iterators
    int I, J, I_, J_; // global numbers of nodes

    // Elementary results memory allocation (they will be correctly initialize for every element K)
    float **A_K = matF_alloc(n_nod_elem, n_nod_elem);
    float *l_K = malloc(n_nod_elem * sizeof(float));
    int *nodes_D = malloc(n_nod_elem * sizeof(int));
    float *uD_aK = malloc(n_nod_elem * sizeof(float));
    float **a_K = malloc(n_nod_elem * sizeof(float *)); // Special case : dynamic tab of float pointers

    /*--------Containers allocation for assembly----------------------------------------------------------------------*/

    // int for an arithmetic pointer operation to know to the next free memory emplacement in A_sparse
    int next_free_A = 1;

    /*--------Assembly------------------------------------------------------------------------------------------------*/
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

        // elementary evaluations
        eval_K(ref_interior, ref_Dh, ref_Dnh, ref_NF, n_Dh, n_Dnh, n_NF, type, n_nod_elem, (const float **) a_K,
               n_edg_elem, ref_edg[k], A_K, l_K, nodes_D, uD_aK);

        // Assembly
        for (i = 0; i < n_nod_elem; i++) {

            // 2nd member assembly
            I = nod_gNb[k][i]-1;
            sec_mb[I] += l_K[i];

            // Dirichlet node handle
            if (nodes_D[i] == 1) {
                D_gb_nod[I] = I+1;    // non-Dirichlet node
            } else if (nodes_D[i] == 0) {
                D_gb_nod[I] = 0;    // homogeneous Dirichlet node
            } else {
                D_gb_nod[I] = -(I+1);   // non-homogeneous Dirichlet node
                D_values[I] = uD_aK[i];
            }

            // Matrix assembly
            for (j = 0; j <= i; j++) {
                J = nod_gNb[k][j]-1;
                if (J==I){ // case of a diagonal term
                    A_sparse[J] += A_K[i][j];
                }
                else { // reorganize so that J <= I
                    if (J > I) {
                        J_ = I+1;
                        I_ = J+1;
                    } else {
                        I_ = I+1;
                        J_ = J+1;
                    } /// @warning indices in Fortran begin to 1 -> we had + 1
                    // Check we haven't reached A_sparse end of memory zone
                    if (next_free_A >= n_lowA) {
                        perror("Warning overflow of sparse matrix array in 'assembly' function (assembly.c) ");
                    }
                    assmat_(&I_, &J_, &(A_K[i][j]), first_non0_row, col_ind, nextInRow, A_sparse + n_nod, &next_free_A);
                }
            }
        }

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