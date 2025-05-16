#include "assembly.h"
#include "forfun.h"
#include "tab_mngmt.h"
#include "elem_eval.h"
#include <stdio.h>
#include <stdlib.h>

/// @brief evaluates elementary matrices and 2nd members for each element and then assemblies everything into a messy morse stored matrix and a 2nd member vector
int assembly(
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
        int *BC_nod,         /// @param BC_nod boundary condition identifier
        float *D_values,
        int *first_non0_row,
        float *A_DMS,     /// (messy Morse Storage) @warning : have to be initialized to 0
        int *col_ind,
        int *nextInRow        /// @param nextInRow next non-zero element in the same row (on the left or right side)
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

    // Initialize BC_nod to I+1 (Neumann)
    for (I = 0; I < n_nod; I++) { BC_nod[I] = I + 1; }

    /*--------Containers allocation for assembly----------------------------------------------------------------------*/

    // int for an arithmetic pointer operation to know to the next free memory emplacement in A_DMS
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
            I = nod_gNb[k][i] - 1;
            sec_mb[I] += l_K[i];

            // Boundary conditions (Neumann-Dirichlet) node handle
            if (nodes_D[i] == 0) { // homogeneous Dirichlet node
                BC_nod[I] = 0;
            } else if (nodes_D[i] == -1) { // non-homogeneous Dirichlet node
                BC_nod[I] = -(I + 1);
                D_values[I] = uD_aK[i];
            }
            // nothing to do for Neumann nodes

            // Matrix assembly
            for (j = 0; j <= i; j++) {
                J = nod_gNb[k][j] - 1;
                if (J == I) { // case of a diagonal term
                    A_DMS[J] += A_K[i][j];
                } else { // reorganize so that J <= I
                    if (J > I) {
                        J_ = I + 1;
                        I_ = J + 1;
                    } else {
                        I_ = I + 1;
                        J_ = J + 1;
                    } /// @warning indices in Fortran begin to 1 -> we had + 1
                    // Check we haven't reached A_DMS end of memory zone
                    if (next_free_A >= n_lowA) {
                        perror("Warning overflow of sparse matrix array in 'assembly' function (assembly.c) ");
                    }
                    assmat_(&I_, &J_, &(A_K[i][j]), first_non0_row, col_ind, nextInRow, A_DMS + n_nod,
                            &next_free_A);
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

    return n_nod+next_free_A; // Length of A_DMS array;
}


/// @brief transform DMS to OMS and takes Dirichlet nodes into condition
void DMS_to_OMS_with_BC(
        const float *A_DMS,
        const int *first_non0_row_DMS,
        const int *col_ind_DMS,
        const int *nextInRow,
        const float *sec_mb,
        const int *BC_nod,
        const float *D_values,
        int n_nod,
        float *A_OMS,
        int *first_non0_row_OMS,
        int *col_ind_OMS,
        float *sec_mb_BC
) {
    cdesse_(&n_nod, first_non0_row_DMS, col_ind_DMS, nextInRow, A_DMS, sec_mb, BC_nod, D_values, first_non0_row_OMS,
            col_ind_OMS, A_OMS, sec_mb_BC);
}


/// @brief compute array's length of profile storage of A
int length_profile(int n_row, const int *first_non0_row_OMS, const int *col_ind_OMS){
    int res = 0;
    for(int i=0; i<n_row-1; i++){
        if(first_non0_row_OMS[i] < first_non0_row_OMS[i+1]) {
            res += i - (col_ind_OMS[first_non0_row_OMS[i] - 1] - 1) + 1;
        }
    }
    res++; //  fictional element
    return res;
}

/// @brief transforms OM storage into profile storage
void OMS_to_profile(
        const float *A_OMS,
        int n_row,     /// @param n_row number of rows in the matrix A_OMS (!= n_nod)
        const int *first_non0_row,
        const int *col_ind_OMS,
        int n_lowA_pr,  /// @param n_lowA_pr length of low_Apr array
        float *A_pr,   /// @param A_pr matrix A stored in the profile way (supposed to be initialized with 0)
        int *profile   /// @param profile array of indices in low part of A_pr of the first non-zero element of each row
) {
    int i, l;

    for(i=0; i<n_lowA_pr; i++) {A_pr[n_row+i] = 0.0; } // Initialisation to 0
    profile[0] = first_non0_row[0];
    A_pr[0] = A_OMS[0]; // diagonal
    for (l = 1; l < n_row-1; l++) {
        A_pr[l] = A_OMS[l]; // diagonal
        if (first_non0_row[l - 1] == first_non0_row[l]) {
            profile[l] = profile[l - 1];
        } else {
            profile[l] = profile[l - 1] + l - col_ind_OMS[first_non0_row[l - 1] - 1] + 1; // -1 + 1
        }
        for (i = 0; i < first_non0_row[l + 1] - first_non0_row[l]; i++) {
            A_pr[n_row + profile[l] - 1 + col_ind_OMS[first_non0_row[l] - 1 + i] -
                 col_ind_OMS[first_non0_row[l] - 1]] = A_OMS[n_row + first_non0_row[l] - 1 + i];
        }
    }
    // last row : only fictional element
    A_pr[n_row-1] = A_OMS[n_row-1]; // diagonal
    profile[n_row-1] = n_lowA_pr;
}