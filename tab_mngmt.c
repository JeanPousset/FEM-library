/**
 * @file tab_mngmt.c
 * @brief Functions that handle memory for dynamically allocated matrices
 */

#include <stdlib.h>

/// Allocate a 2-dimensional tab (matrix) of float dynamically with value 0.0 for each element
float **matF_alloc0(int dim1, int dim2) {
    float **ptr;

    ptr = malloc(dim1*sizeof(float *));
    if (ptr != NULL) {
        int i, taille_ligne = dim2*sizeof(float);
        float *tmp = calloc(dim1,taille_ligne); // only this line change
        if (tmp != NULL) {
            for (i=0; i<dim1; i++) {
                ptr[i] = tmp;
                tmp += dim2;
            }
        }
        else {
            free(ptr);
            ptr = NULL;
        }
    }
    return(ptr);
}

/// Allocate a 2-dimensional tab (matrix) of float dynamically
float **matF_alloc(int dim1, int dim2) {
    float **ptr;

    ptr = malloc(dim1*sizeof(float *));
    if (ptr != NULL) {
        int i, taille_ligne = dim2*sizeof(float);
        float *tmp = malloc(dim1*taille_ligne);
        if (tmp != NULL) {
            for (i=0; i<dim1; i++) {
                ptr[i] = tmp;
                tmp += dim2;
            }
        }
        else {
            free(ptr);
            ptr = NULL;
        }
    }
    return(ptr);
}

/// Allocate a 2-dimensional tab (matrix) of integer dynamically
int **matI_alloc(int dim1, int dim2) {
    int **ptr;

    ptr = malloc(dim1*sizeof(int *));
    if (ptr != NULL) {
        int i, taille_ligne = dim2*sizeof(int);
        int *tmp = malloc(dim1*taille_ligne);
        if (tmp != NULL) {
            for (i=0; i<dim1; i++) {
                ptr[i] = tmp;
                tmp += dim2;
            }
        }
        else {
            free(ptr);
            ptr = NULL;
        }
    }
    return(ptr);
}


/// Free memory that has been allocated by matF_alloc / matI_alloc
void free_mat(void *ptr) {
    void **ptrT=ptr;
    free(ptrT[0]);
    free(ptr);
}
