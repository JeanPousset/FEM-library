/**
 * @file meshing.h
 * @brief Header of functions to handle mesh files (read & write)
 * @warning The mesh files have to be in a specific format
 * @warning The domain is supposed to be 2D and rectangular
 */
#ifndef FEM_PROJECT_MESHING_H
#define FEM_PROJECT_MESHING_H

/**
 * @brief Fill the 2 dim tab refEdg with reference number for each edge of each element
 *
 * @param type         type of the element (1: quadrangle, 2:triangle)
 * @param n1           number of node on the x-axis of the domain
 * @param n2           number of node on the x-axis of the domain
 * @param ref_interior reference number for edges inside the domain (usually 0)
 * @param ref_sides    array of reference numbers for each domain sides
 * @param n_elem       number of element of the mesh
 * @param n_edg_elem   number of edges per element
 * @param refEdg       result : matrix of reference number for each edge of each element
 */
void etiqAr(int type, int n1, int n2, int ref_interior, const int *ref_sides, int n_elem, int n_edg_elem, int **refEdg);


/**
 * @brief Write a mesh file based on domain discretization parameters
 *
 * @param a             x_coordinate of the domain's bottom right corner
 * @param b             x_coordinate of the domain's bottom left corner
 * @param c             y_coordinate of the domain's bottom right corner
 * @param d             y_coordinate of the domain's top right corner
 * @param n1            number of node on the x-axis of the domain
 * @param n2            number of node on the y-axis of the domain
 * @param t             type of the element (1: quadrangle, 2:triangle)
 * @param ref_interior  reference number for edges inside the domain (usually 0)
 * @param ref_sides     array of reference numbers for each domain sides
 * @param mesh_out_path string for the relative path of the mesh file to write
 */
void write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int ref_interior, const int* ref_sides, char* mesh_out_path);


/**
 * @brief Read a mesh file and assign parameters values, coordinates, and edges references numbers
 *
 * @warning Every parameter to be assigned and the two 2 dimensional tabs (coords, gNb_node, refEdg) are passed as pointers .
 * @warning p_coords, p_refEdg and p_nod_gNb are supposed to be empty pointers that aren't allocated yet
 *
 * @param input_mesh mesh file's relative path
 * @param type       (pointer to) type of the element (1: quadrangle, 2:triangle)
 * @param n_nod      (pointer to) total number of nodes
 * @param p_coords   (pointer to) nodes' coordinates
 * @param n_elem     (pointer to) number of element
 * @param p_nod_gNb  (pointer to) matrix that contains nodes' global number for each element
 * @param n_nod_elem (pointer to) number of nodes per element
 * @param n_edg_elem (pointer to) number of edges per element
 * @param p_refEdg   (pointer to) matrix that contains edges reference number for each element
 *
 * @return 0: successfully reading, -1 : error when opening the mesh file
 */
int read_mesh(const char *input_mesh, int *type, int* n_nod, float*** p_coords, int* n_elem, int*** p_nod_gNb, int *n_nod_elem, int* n_edg_elem, int*** p_refEdg);

#endif //FEM_PROJECT_MESHING_H
