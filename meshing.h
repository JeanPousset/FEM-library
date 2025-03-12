#ifndef FEM_PROJECT_MESHING_H
#define FEM_PROJECT_MESHING_H


/*
    Fill the 2 dim tab refEdg with reference number for each edge of each element
*/
void etiqAr(int typel, int n1, int n2, int nrefdom, const int *nrefcot, int nbtel, int nbaret, int **refEdg);

void write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int ref_interior, const int* ref_sides, char* mesh_out_path);

/*
    Read a mesh file and assign parameters values, coordinates, and edges references numbers :
    Every parameter to be assigned and the two 2 dimensionnal tabs (coords, gNb_node, refEdg) are passed as pointers .
    !!
    pcoords, p_refEdg and p_nod_gNb are supposed to be empty pointers that aren't allocated yet
    !!
*/
int read_mesh(char* input_mesh, int *type, int* nodes_gN, float*** p_coords, int* n_elem, int*** p_nod_gNb, int *n_nod_elem, int* n_edg_elem, int*** p_refEdg);
#endif //FEM_PROJECT_MESHING_H
