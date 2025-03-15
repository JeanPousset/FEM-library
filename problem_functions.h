#ifndef FEM_PROJECT_PROBLEM_FUNCTIONS_H
#define FEM_PROJECT_PROBLEM_FUNCTIONS_H

float a00(const float* x) { return 0;}

float a11(const float* x) { return 1;}

float a12(const float* x) { return 0;}

float a22(const float* x) { return 1;}

float bN(const float* x) { return 0;}

float fN(const float* x) { return 1;}

float f_OMEGA(const float* x) { return 1;}

float uD(const float* x) { return 100*x[0] + x[1];}

#endif //FEM_PROJECT_PROBLEM_FUNCTIONS_H
