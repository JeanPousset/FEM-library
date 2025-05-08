/**
 * @file problem_functions.h
 * @brief Implementations of mathematical function associated to the physical problem to be addressed
 */

#ifndef FEM_PROJECT_PROBLEM_FUNCTIONS_H
#define FEM_PROJECT_PROBLEM_FUNCTIONS_H

#include <math.h>
#include <stdio.h>

float a00(const float* x);

float a11(const float* x);

float a12(const float* x);

float a22(const float* x);

float bN(const float* x);

float fN(const float* x);

float f_OMEGA(const float* x);

float uD(const float* x);

/// @brief exact solution function
extern int nucas;
float exact_sol(const float *x);


#endif //FEM_PROJECT_PROBLEM_FUNCTIONS_H
