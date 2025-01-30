#include "meshing.h"
#include "unit_tests.h"

int main(int argc, char* argv[])
{
    /* // Test write_mesh
    char* mesh_file = "mesh_f1.txt";

    float a = 0.0;
    float b = 1.0;
    float c = 0.0;
    float d = 1.0;

    int n1 = 5;
    int n2 = 5;
    int t = 1; // Quadrangle

    int nrefdom = 0;
    const int nrefcot[] = {1, 2, 3, 4};
    write_mesh(a,b,c,d,n1,n2,t,nrefdom,nrefcot,mesh_file);
    */

    // Test write_mesh
    char* real_mesh_file = "../car1x1q_4";
    char* check_read_mesh_output = "../check_read_car1x1q_4";
    test_read_mesh(real_mesh_file,check_read_mesh_output);

    return 0;


}