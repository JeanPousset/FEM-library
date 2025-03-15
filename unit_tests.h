#ifndef FEM_PROJECT_UNIT_TEST_H
#define FEM_PROJECT_UNIT_TEST_H


// Test the inv of a matrix
int test_inv2x2(void);

// Test the read_mesh implemented in "meshing.c/h", read an ouput mesh file and rewrite its content into the output mesh file
void test_read_mesh(char *input_mesh, char *out_check_mesh);
/*
    Function that lunch a diff command on the terminal to compare 2 file contents.
    If there isn't any difference, nothing is printed to the screen.
    If you see something, file contents are different
*/

void test_write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int nrefdom, const int *nrefcot,char *mesh_file_path, char *correct_mesh_file);
/*
    Test the write_mesh function implemented in "meshing.c"
    It compares the result of write_mesh (mesh_file_path) with a correct one (correct_mesh_file) of the same meshing
*/

void diff_file(const char* f1,const char* f2);
/*
    Function that lunch a diff command on the terminal to compare 2 file contents.
    If there isn't any difference, nothing is printed to the screen.
    If you see something, file contents are different
*/

void print_eval_K(int K, int type, int n_nod_elem, float **A_K, float *l_K,
              int *nodes_D, float *uD_aK);
/*
    Print results of A_K, l_K (2nd member), Dirichlet conditions in each node and values of non-homogeneous Dirichlet nodes
[IN]
    - K     : element number
    - type  : type of element (e.g. 2 for triangle)
    - n_nod_elem : number of nodes per element
    - **A_K      : Matrix of the elementary values of A_ij_K for each couple of nodes (i,j)  [n_nod_elem * n_nod_elem]
    - **l_K      : vector of elementary values of l_K_i for each node i [n_nod_elem]
    - *nodes_D   : array where element i is equal to -1 if node i of K is Dirichlet, 0 if Dirichlet homogeneous and -1 otherwise
    - *uD_aK     : array that contains values of the non-homogeneous Dirichlet nodes
*/

// Test the eval_K function on every mesh element and call print_eval_K so we can compare to the exact solution
void test_eval_K(const char* mesh_file);

#endif //FEM_PROJECT_UNIT_TEST_H
