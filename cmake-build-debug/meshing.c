#include "meshing.h"
#include <stdio.h>
#include <stdlib.h>

// Create a discretisation of [a,b] with n points
float* linspace(float a, float b, int n)
{
    int i;
    float* res = malloc(n*sizeof(float));
    for(i=0;i<n;i++) res[i] = (a-b)/n * i;
    return res;
}

void write_mesh(float a, float b, float c, float d, int n1, int n2, int t, const int* ref, char* mesh_out_path)
{
    int i,j,k;
    int m,p,q;

    // Opening output file in writing mode
    FILE *mesh_f = fopen(mesh_out_path, "w");	// mode 'w' for writing (delete everything there was before)
    if (mesh_f == (FILE*)NULL) {
        printf("Error when trying to open the output (mesh) file : '%s' in write_mesh_from_file function", mesh_file_out);
        perror("^^ IO - Error ^^");
    }

    // write value of n (nb of nodes)
    fprintf(mesh_f,"%d\n",n1*n2);

    // Write nodes coordinates x_i y_i
    for(j=0; j<n2; j++){
        for (i=0; i<n1; i++){
            fprintf(mesh_f,"%f %f\n",(b-a)/n1*i,(d-c)/n2*j);
        }
    }

    // Next we depend on geometry type :
    switch (t){

        case 1 : // Quadrande
            p = 4;
            q = 4;
            m = (n1-1)*(n2-2);
            fprintf(mesh_f, "%d %d %d %d\n",m,t,p,q);

            // Write vertices
            for(j=0; j<n2-1; j++){
                for (i=0; i<n1-1; i++) { // for each quadrande
                    fprintf(mesh_f, "%d ", i * n1 + j * n2 + 2);       // vertice 1
                    fprintf(mesh_f, "%d ", i * (n1 + 1) + j * n2 + 2); // vertice 2
                    fprintf(mesh_f, "%d ", i * (n1 + 1) + j * n2 + 1); // vertice 3
                    fprintf(mesh_f, "%d ", i * n1 + j * n2 + 1);       // vertice 4
                }
            }
            break;

        case 2 : // Triangle
            p = 3;
            q = 3;
            m = 2*(n1-1)*(n2-2);
            fprintf(mesh_f, "%d %d %d %d\n",m,t,p,q);

            // Write vertices
            for(j=0; j<n2-1; j++){
                for (i=0; i<n1-1; i++) {
                    // 1st triangle
                    fprintf(mesh_f, "%d ", i * n1 + j * n2 + 2);       // vertice 1
                    fprintf(mesh_f, "%d ", i * (n1 + 1) + j * n2 + 1); // vertice 2
                    fprintf(mesh_f, "%d ", i * n1 + j * n2 + 1);       // vertice 3
                    // 2nd triangle
                    fprintf(mesh_f, "%d ", i * (n1 + 1) + j * n2 + 1); // vertice 1
                    fprintf(mesh_f, "%d ", i * n1 + j * n2 + 2);       // vertice 2
                    fprintf(mesh_f, "%d ", i * (n1 + 1) + j * n2 + 2); // vertice 3
                }
            }
            break;
        default:
            printf("Unknown type t = %d for 'write_mesh' function");
            perror("Error : Unvalid parameter");
    }

    fclose(mesh_f); // Close mesh output file
}
