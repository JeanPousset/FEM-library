#ifndef FEM_MESHING_H
#define FEM_MESHING_H

/*
    Create a discretisation of [a,b] with n points
    as a dynamic allocated array of size (n)
*/
float* linspace(float a, float b, int n);


void write_mesh(float a, float b, float c, float d, int n1, int n2, int t, const int* ref, char* mesh_file_out);


#endif //FEM_MESHING_H
