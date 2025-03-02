#include "unit_tests.h"
#include "meshing.h"
#include "tab_mngmt.h"
#include "elem_eval.h"
#include <stdio.h>



// Test the inv of a matrix !! Not completed -> memory leaks
int test_inv2x2(void)
{
    int i,j;
    float** M = matF_alloc(2,2);
    float** M_inv = matF_alloc(2,2);
    float** M_cor = matF_alloc(2,2);
    float det;
    M[0][0] = 3; M[0][1] = 2;
    M[1][0] = 2; M[1][1] = 2;
    M_cor[0][0] = 1;  M_cor[0][1] = -1;
    M_cor[1][0] = -1; M_cor[1][1] = 1.5;

    det = inv_2x2(M,M_inv);

    if (det != 2.0) return 0;
    // Checking if the matrix is correctly inverted
    for(i=0; i<2; i++) for(j=0; j<2; j++) if (M_cor[i][j] != M_inv[i][j]) return 0;
    printf(" 'inv2x2' function passed test with success \n");
    return 1;
}











// Function that lunch a diff command on the terminal to compare 2 file contents.
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



//    Test the write_mesh function implemented in "meshing.c"
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

    // clean tab memory
    freetab(p_coords);
    freetab(p_gNb_node);
    freetab(p_refEdg);

    // Using diff command line to test if there is a difference between the original and the file
    diff_file(input_mesh,out_check_mesh);
}

