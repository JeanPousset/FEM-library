#ifndef FEM_PROJECT_ASSEMBLY_H
#define FEM_PROJECT_ASSEMBLY_H

#define NON_0_ROW 12

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
        float *sec_mb,
        int *D_gb_nod,
        float *D_values,
        int *first_non0_row,
        float *A_sparse,
        int *col_ind,
        int *nextInRow
);


#endif //FEM_PROJECT_ASSEMBLY_H
