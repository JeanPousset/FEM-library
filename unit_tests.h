#ifndef FEM_PROJECT_UNIT_TEST_H
#define FEM_PROJECT_UNIT_TEST_H

// Test the read_mesh implemented in "meshing.c/h", read an ouput mesh file and rewrite its content into the output mesh file
void test_read_mesh(char *input_mesh, char *out_check_mesh);

// Test the read_mesh implemented in "meshing.c". Compare the result with a correct mesh file of the same meshing
void test_write_mesh(float a, float b, float c, float d, int n1, int n2, int t, int nrefdom, const int *nrefcot,char *mesh_file_path, char *correct_mesh_file);


// Function that lunch a diff command on the terminal to compare 2 file contents.
void diff_file(const char* f1,const char* f2);

#endif //FEM_PROJECT_UNIT_TEST_H
