#ifndef FEM_MESHING_H
#define FEM_MESHING_H

/*
    Create a discretisation of [a,b] with n points
    as a dynamic allocated array of size (n)
*/
float* linspace(float a, float b, int n);

/*
    Fill the 2 dim tab refEdg with reference number for each edge of each element
*/
void etiqAr(int typel, int n1, int n2, int nrefdom, const int *nrefcot, int nbtel, int nbaret, int **refEdg);

void write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int nrefdom, const int* nrefcot, char* mesh_out_path);

/*
    Read a mesh file and assign parameters values, coordinates, and edges references numbers :
    Every parameter to be assigned and the two 2 dimensionnal tabs (coords, gNb_node, refEdg) are passed as pointers .
    !!
    pcoords, p_refEdg and p_gNb_node are supposed to be empty pointers that aren't allocated yet
    !!
*/
int read_mesh(char* input_mesh, int *type, int* n_nodes, float*** p_coords, int* n_elem, int*** p_gNb_node, int *n_nod_elem, int* n_edges, int*** p_refEdg);
#endif //FEM_MESHING_H
