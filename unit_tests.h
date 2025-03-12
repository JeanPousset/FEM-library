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

void imp_eval_K(int K, int type, int n_nod_elem, float **A_K, float *l_K,
              int *nodes_D, float *uD_aK);
/************************************************************************
  Imprime les resultats de la matrice et du second membre elementaires
  ainsi que les conditions Dirichlet en chaque noeud
  et les valeurs des conditions Dirichlet non homogene

*** Arguments ***
   K        : Numero de l'element
   typEl    : Numero de type de l'element
   nbneel   : Nombre de noeuds de l'element
   MatElem  : Matrice elementaire de dimensions (nbneel,nbneel)
   SMbrElem : Second membre elementaire de dimension nbneel
   NuDElem  : Tableau de reperage des noeuds porteurs de conditions de Dirichlet
   uDElem   : Tableau des valeurs de blocage
              pour les noeuds Dirichlet non homogene
************************************************************************/
#endif //FEM_PROJECT_UNIT_TEST_H
