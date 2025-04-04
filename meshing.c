/**
 * @file meshing.c
 * @brief Functions to handle mesh files (read & write)
 * @warning The mesh files have to be in a specific format
 * @warning The domain is supposed to be 2D and rectangular
 */
#include "meshing.h"
#include "tab_mngmt.h"
#include <stdio.h>


/// Fill the 2 dim tab refEdg with reference number for each edge of each element
void etiqAr(int type, int n1, int n2, int ref_interior, const int *ref_sides,int n_elem, int n_edg_elem, int **refEdg)
{
    int i,j;

    // Fill everything at r0
    for (i=0; i<n_elem; i++) {
        for (j=0; j<n_edg_elem; j++){
            refEdg[i][j] = ref_interior;
        }
    }

    switch (type) {

        case 1 : // Quadrangle
            // Lower and upper bound of the domain
            for (i=0; i<n1-1; i++){
                refEdg[i][3] = ref_sides[0];
                refEdg[(n1-1)*(n2-2)+i][1] = ref_sides[2];
            }
            // Right and left bound of the domain
            for (i=0; i < n2-1 ; i++){
                refEdg[i*(n1-1)][2] = ref_sides[3];
                refEdg[i*(n1-1)+n1-2][0] = ref_sides[1];
            }
            break;

        case 2 : // Triangle
            // Lower and upper bound of the domain
            for (i=0; i < 2*(n1-1); i+=2){
                refEdg[i][2] = ref_sides[0];
                refEdg[2*(n1-1)*(n2-2)+i+1][2] = ref_sides[2];
            }
            // Right and left bound of the domain
            for (i=0; i < 2*(n2-1)-1 ; i+=2){
                refEdg[i*(n1-1)][1] = ref_sides[3];
                refEdg[i*(n1-1)+2*(n1-1)-1][1] = ref_sides[1];
            }
            break;

        default:
            printf("Unknown type t = %d for 'etiqAr' function",type);
            perror("Error : Unvalid parameter");
    }
}

/// Write a mesh file based on domain discretization parameters
void write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int ref_interior, const int* ref_sides, char* mesh_out_path)
{
    int i,j,k;
    int m,p,q;
    int** refEdg;

    // Opening output file in writing mode
    FILE *mesh_f = fopen(mesh_out_path, "w");	// mode 'w' for writing (delete everything there was before)
    if (mesh_f == (FILE*)NULL) {
        printf("Error when trying to open the output (mesh) file : '%s' in write_mesh function", mesh_out_path);
        perror("^^ IO - Error ^^");
    }

    // write nb of nodes
    fprintf(mesh_f,"%d\n",n1*n2);
    // Write nodes coordinates x_i y_i
    for(j=0; j<n2; j++)
        for (i=0; i<n1; i++)
            fprintf(mesh_f,"%f %f\n", i*(b-a)/(n1-1), j*(d-c)/(n2-1));


    // Next we depend on geometry type :
    switch (t){

        case 1 : // Quadrangle
            p = 4;
            q = 4;
            m = (n1-1)*(n2-1);
            fprintf(mesh_f, "%d %d %d %d",m,t,p,q);

            // fill table of reference number for each edge of each element
            refEdg = matI_alloc(m,q);
            etiqAr(t,n1,n2,ref_interior,ref_sides,m,q,refEdg);

            // Write vertices
            for(j=0; j<n2-1; j++){
                for (i=0; i<n1-1; i++) { // for each quadrande
                    fprintf(mesh_f, "\n%d ", j * n1 + i + 2);   // vertice 1
                    fprintf(mesh_f, "%d ", (j+1) * n1 + i + 2); // vertice 2
                    fprintf(mesh_f, "%d ", (j+1) * n1 + i + 1); // vertice 3
                    fprintf(mesh_f, "%d ", j * n1 + i + 1);      // vertice 4
                    // edges references number
                    for (k=0; k<q; k++) fprintf(mesh_f, "%d ", refEdg[j*(n1-1)+i][k]);
                }
            }
            free_mat(refEdg);
            break;

        case 2 : // Triangle
            p = 3;
            q = 3;
            m = 2*(n1-1)*(n2-1);
            fprintf(mesh_f, "%d %d %d %d",m,t,p,q);

            // fill table of reference number for each edge of each element
            refEdg = matI_alloc(m,q);
            etiqAr(t,n1,n2,ref_interior,ref_sides,m,q,refEdg);

            // Write vertices
            for(j=0; j<n2-1; j++){
                for (i=0; i<n1-1; i++) {
                    // 1st triangle
                    fprintf(mesh_f, "\n%d ", j * n1 + i + 2);   // vertice 1
                    fprintf(mesh_f, "%d ", (j+1) * n1 + i + 1); // vertice 2
                    fprintf(mesh_f, "%d ", j * n1 + i + 1);      // vertice 3
                    // edges references number for 1st triangle
                    for (k=0; k<q; k++) fprintf(mesh_f, "%d ", refEdg[2*j*(n1-1)+2*i][k]);
                    // 2nd triangle
                    fprintf(mesh_f, "\n%d ", (j+1) * n1 + i + 1); // vertice 1
                    fprintf(mesh_f, "%d ", j * n1 + i + 2);       // vertice 2
                    fprintf(mesh_f, "%d ", (j+1) * n1 + i + 2);    // vertice 3
                    // edges references number for 2nd triangle
                    for (k=0; k<q; k++) fprintf(mesh_f, "%d ", refEdg[2*j*(n1-1)+2*i+1][k]);
                }
            }
            free_mat(refEdg);
            break;
        default:
            printf("Unknown type t = %d for 'write_mesh' function",t);
            perror("Error : Unvalid parameter");
    }
    fprintf(mesh_f, "\n"); // One list new line to have exactly the same shape as other mesh files
    fclose(mesh_f); // Close mesh output file
}

/// Read a mesh file and assign parameters values, coordinates, node global numbers, and edges references numbers
int read_mesh(const char *input_mesh, int *type, int* n_nod, float*** p_coords, int* n_elem, int*** p_nod_gNb, int *n_nod_elem, int* n_edg_elem, int*** p_refEdg)
{
    int i,j;

    // Opening meshfile in read mode
    FILE *mesh_f = fopen(input_mesh, "r");	// mode 'w' for reading
    if (mesh_f == (FILE*)NULL) {
        printf("Error when trying to open the input (mesh) file : '%s' in read_mesh function", input_mesh);
        perror("^^ IO - Error ^^");
        return -1;
    }

    fscanf(mesh_f,"%d",n_nod);        // number of nodes (n)
    *p_coords = matF_alloc(*n_nod,2); // allocate a n x 2 tab for coordinates

    // fill coordinates x_i y_i
    for (i=0; i<*n_nod; i++) fscanf(mesh_f,"%f %f",&(*p_coords)[i][0], &(*p_coords)[i][1]);

    // assign nb element (m), type (t), nb of geometrical node per element (p), nb of egdes per element (q);
    fscanf(mesh_f,"%d %d %d %d",n_elem,type,n_nod_elem,n_edg_elem);

    // fill global node number reference for each element and edges reference number
    *p_nod_gNb = matI_alloc(*n_elem,*n_nod_elem);
    *p_refEdg = matI_alloc(*n_elem,*n_edg_elem);
    for (i=0; i<*n_elem; i++){
        for (j=0; j<*n_nod_elem; j++) fscanf(mesh_f,"%d ",&(*p_nod_gNb)[i][j]);
        for (j=0; j<*n_edg_elem; j++) fscanf(mesh_f,"%d ",&(*p_refEdg)[i][j]);
    }
    fclose(mesh_f); // Close mesh output file
    return 0;
}
