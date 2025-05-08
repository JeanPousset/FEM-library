#ifndef FEM_PROJECT_RESULTS_H
#define FEM_PROJECT_RESULTS_H

/// @brief compute the array of exact solution for each point given in pts_coords
void eval_exact_sol(int n_row, const float **pts_coords, float *Uex);

/// @brief compute finite element solution by solving Linear problem A*U = 2nd_mb with Cholesky
void solve_cholesky(const float *A_pr, const int *profile, int n_row, int n_lowA_pr, const float *sec_mb_BC, float *u);



#endif //FEM_PROJECT_RESULTS_H
