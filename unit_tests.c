#include "unit_tests.h"
#include "meshing.h"
#include <stdio.h>

// Test the read_mesh form "meshing.c/h", read an ouput mesh file and rewrite its content into the output mesh file
void test_read_mesh(char *input_mesh, char *out_check_mesh) {
    int type, n_nodes, n_elem, n_nod_elem, n_edges;
    float **p_coords;
    int **p_gNb_node;
    int **p_refEdg;

    // reading input_meshfile
    int res = read_mesh(input_mesh,&type,&n_nodes,&p_coords,&n_elem,&p_gNb_node,&n_nod_elem,&n_edges,&p_refEdg);
    // Opening ouput meshfile in writting mode
    FILE *mesh_f = fopen(out_check_mesh, "w");	// mode 'w' for writing (delete everything there was before)
    if (mesh_f == (FILE*)NULL) {
        printf("Error when trying to open the output (mesh) file : '%s' in test_read_mesh function", out_check_mesh);
        perror("^^ IO - Error ^^");
    }

    // Writing output (test) meshfile
    int i,j;
    fprintf(mesh_f,"%d",n_nodes);
    // coordinates
    for (i=0; i<n_nodes; i++) fprintf(mesh_f,"\n%f %f",p_coords[i][0],p_coords[i][1]);
    // nodes global number and edges references number
    for (i=0; i<n_elem; i++){
        fprintf(mesh_f,"\n");
        for (j=0; j<n_nod_elem; j++) fprintf(mesh_f,"%d ",p_gNb_node[i][j]);
        for (j=0; j<n_edges; j++) fprintf(mesh_f,"%d ",p_refEdg[i][j]);
    }
    fclose(mesh_f);
}
