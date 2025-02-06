#include "unit_tests.h"
#include "meshing.h"
#include <stdio.h>

/*
    Function that lunch a diff command on the terminal to compare 2 file contents.
    If there isn't any difference, nothing is printed to the screen.
    If you see something, file contents are different
*/
void diff_file(const char* f1,const char* f2)
{
    char diff_command[512];
    // The file compare command is different in Windows OS
#ifdef _WIN32
    snprintf(diff_command, sizeof(diff_command), "FC %s %s", f1, f2);
#else
    snprintf(diff_command, sizeof(diff_command), "diff %s %s", f1, f2);
#endif
    FILE *fp = popen(diff_command, "r");
    if (fp == NULL) perror("Error when opening a pipe");

    char answer_buffer[256];
    while (fgets(answer_buffer, sizeof(answer_buffer), fp) != NULL) printf("%s", answer_buffer);
    pclose(fp);
}


/*
    Test the write_mesh function implemented in "meshing.c"
    It compares the result of write_mesh (mesh_file_path) with a correct one (correct_mesh_file) of the same meshing
*/
void test_write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int nrefdom, const int *nrefcot,char *mesh_file_path, char *correct_mesh_file) {
    write_mesh(a, b, c, d, n1, n2, t, nrefdom, nrefcot, mesh_file_path);
    diff_file(mesh_file_path,correct_mesh_file);
}

// Test the read_mesh form "meshing.c/h", read an ouput mesh file and rewrite its content into the output mesh file
void test_read_mesh(char *input_mesh, char *out_check_mesh) {
    int type, n_nodes, n_elem, n_nod_elem, n_edges;
    float **p_coords;
    int **p_gNb_node;
    int **p_refEdg;

    // reading input_meshfile
    read_mesh(input_mesh,&type,&n_nodes,&p_coords,&n_elem,&p_gNb_node,&n_nod_elem,&n_edges,&p_refEdg);
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

    fprintf(mesh_f, "\n%d %d %d %d",n_elem,type,n_nod_elem,n_edges);
    for (i=0; i<n_elem; i++){
        fprintf(mesh_f,"\n");
        for (j=0; j<n_nod_elem; j++) fprintf(mesh_f,"%d ",p_gNb_node[i][j]);
        for (j=0; j<n_edges; j++) fprintf(mesh_f,"%d ",p_refEdg[i][j]);
    }
    fprintf(mesh_f, "\n"); // One list new line to have exactly the same shape as other mesh files
    fclose(mesh_f);

    // Using diff command line to test if there is a difference between the original and the file
    diff_file(input_mesh,out_check_mesh);
}

