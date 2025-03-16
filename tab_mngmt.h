/**
 * @file tab_mngmt.h
 * @brief Header for functions that handle memory for dynamically allocated matrices
 */

#ifndef FEM_PROJECT_TAB_MNGMT_H
#define FEM_PROJECT_TAB_MNGMT_H

/**
 * @brief  This function allocates memory to store a matrix of dimensions dim1 x dim2 of type float
 *
 * @details The function allocates an array of pointers (ptr)
  each of whose dim1 elements points to a zone of dim2 elements.

  The function returns NULL if an error occurs during allocation.

  The allocated memory can be freed by calling freetab.
  Note:
    The useful elements of the matrix are arranged consecutively in memory,
    in dim1 blocks of dim2 elements each. It follows that this set can
    be transmitted to a Fortran procedure via the value of ptr[0], as a single
    array of dimensions (dim2, dim1). We therefore have the following addressing correspondence
    (ptr[i][j] and tab(j,i) designate the same element):
      - C language: ptr[i][j], i=0,dim1-1, j=0,dim2-1,
    Matrix elements are then referred to using the two-index notation indicated above.

 * @param dim1 number of rows
 * @param dim2 number of columns
 *
 * @return double-pointer to handle the matrix
 */
float **matF_alloc(int dim1, int dim2);

/**
 * @brief  This function allocates memory to store an integer matrix of dimensions dim1 x dim2 (element type : int)
 *
 * @details Same as matF_alloc (@see matF_alloc) but for int
 *
 * @param dim1 number of rows
 * @param dim2 number of columns
 *
 * @return double-pointer to handle the matrix
 */
int **matI_alloc(int dim1, int dim2);


/**
 * @brief  This function allocates memory to store a null (each element = 0) matrix of dimensions dim1 x dim2 of type float
 *
 * @details Same as matF_alloc (@see matF_alloc) but every element of the matrix is initialized to 0
 *
 * @param dim1 number of rows
 * @param dim2 number of columns
 *
 * @return double-pointer to handle the matrix
 */
float **matF_alloc0(int dim1, int dim2);



/**
 * @brief free memory that has been allocated by matF_alloc, matI_alloc or matF_alloc0
 *
 * @warning Not use for 1-dimensional container !
 *
 * @param ptr double-pointer that handles the matrix to be cleaned
 */
void free_mat(void *ptr);

#endif //FEM_PROJECT_TAB_MNGMT_H
