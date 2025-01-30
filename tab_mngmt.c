#include <stdlib.h>

// Allocate a 2 dimensionnal tab dynamically
float **alloctab(int dim1, int dim2) {
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

int **alloctabI(int dim1, int dim2) {
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

#include <stdlib.h>


// Free memory that has been allocated by alloctab / alloctabI
void freetab(void *ptr) {
    void **ptrT=ptr;
    free(ptrT[0]);
    free(ptr);
}
