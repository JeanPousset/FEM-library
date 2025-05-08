/**
 * @file unit_tests.h
 * @brief Header for unit tests
 * @details Contains test-functions that check if functions are correctly working
 */

#ifndef FEM_PROJECT_UNIT_TEST_H
#define FEM_PROJECT_UNIT_TEST_H


/**
 * @brief Check if a 2x2 matrix is correctly inverted
 */
int test_inv2x2(void);

/**
 * @brief Test the read_mesh implemented in "meshing.c/h", read an ouput mesh file and rewrite its content into the output mesh file
 *
 * @details Function that lunch a diff command on the terminal to compare 2 file contents.
    If there isn't any difference, nothing is printed to the screen.
    If you see something, file contents are different
 */
void test_read_mesh(char *input_mesh, char *out_check_mesh);


/**
 * @brief Compares the result of write_mesh (mesh_file_path) with a correct one (correct_mesh_file) of the same meshing
 */
void test_write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int nrefdom, const int *nrefcot,char *mesh_file_path, char *correct_mesh_file);


/**
 * @brief Function that lunch a diff command on the terminal to compare 2 file contents.
 * @details If there isn't any difference, nothing is printed to the screen.
    If you see something, file contents are different
*/
void diff_file(const char* f1,const char* f2);


/**
 * @brief Print results of A_K, l_K (2nd member), Dirichlet conditions in each node and values of non-homogeneous Dirichlet nodes
 *
 * @param K          element number
 * @param type       type of element (1:quadrangle, 2: triangle)
 * @param n_nod_elem number of nodes per element
 * @param A_K        matrix of the elementary values of A_ij_K for each couple of nodes (i,j)  [n_nod_elem * n_nod_elem]
 * @param l_K        vector of elementary values of l_K_i for each node i [n_nod_elem]
 * @param nodes_D    array where element i is equal to -1 if node i of K is Dirichlet, 0 if Dirichlet homogeneous and -1 otherwise
 * @param uD_aK      array that contains values of the non-homogeneous Dirichlet nodes
 */
void print_eval_K(int K, int type, int n_nod_elem, float **A_K, float *l_K,
              int *nodes_D, float *uD_aK);


/**
 * @brief Test the eval_K function on every mesh element and call print_eval_K so we can compare it to the exact solution
 */
void test_eval_K(const char* mesh_file);

/**
 * @brief Test the assembly of the global matrix and 2nd member from the evaluation on each element K
 */

void test_assembly(const char *mesh_file);

/**
 * @brief Test assembly then conversion DMS -> OMS + BC taken into account
 */
void test_DMS_to_OMS_with_BC(const char *mesh_file);

/**
 * @brief Assemblies, then transform into profile storage
 * @warning in these function we allocate the result INSIDE !
 * @return number of rows of the linear system matrix
 */
int linear_system(const char *mesh_file, float **A_pr, int **profile, float **sec_mb_BC);

/**
 * @brief Gets the profile storage of A and prints it
 */
void test_OMS_to_profile(const char *mesh_file);

#endif //FEM_PROJECT_UNIT_TEST_H
