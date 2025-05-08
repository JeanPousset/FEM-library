#ifndef FEM_PROJECT_ASSEMBLY_H
#define FEM_PROJECT_ASSEMBLY_H

#define NON_0_ROW 12

///@brief
///@return length of the A_DMS array;
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
        float *sec_mb,
        int *D_gb_nod,
        float *D_values,
        int *first_non0_row,
        float *A_sparse,
        int *col_ind,
        int *nextInRow
);

///@brief count number of row in A_OMS after BC are taken in condition
int n_row_with_BC(int n_nod,const int *BC_nod);


/// @brief take the messy (Disordered) matrix, takes Dirichlet nodes into condition and give the Ordered MS (ODS) matrix
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
);

/// @brief compute array's length of profile storage of A
int length_profile(int n_row, const int *first_non0_row_OMS, const int *col_ind_OMS);


/// @brief transforms OM storage into profile storage
void OMS_to_profile(
        const float *A_OMS,
        int n_row,     /// @param n_row number of rows in the matrix A_OMS (!= n_nod)
        const int *first_non0_row,
        const int *col_ind_OMS,
        float *A_pr,   /// @param A_pr matrix A stored in the profile way (supposed to be initialized with 0)
        int *profile   /// @param profile array of indices in low part of A_pr of the first non-zero element of each row
);


#endif //FEM_PROJECT_ASSEMBLY_H
