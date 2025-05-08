#include "results.h"
#include "problem_functions.h"
#include "forfun.h"
#include <stdlib.h>

/// @brief compute the array of exact solution for each point given in pts_coords
void eval_exact_sol(int n_row, const float **pts_coords, float *Uex) {
    for (int I = 0; I < n_row; I++) {
        Uex[I] = exact_sol(pts_coords[I]);
    }
}

/// @brief compute finite element solution by solving Linear problem A*U = 2nd_mb with Cholesky
void solve_cholesky(const float *A_pr, const int *profile, int n_row, int n_lowA_pr, const float *sec_mb_BC, float *u) {
    float eps; // singularity threshold of A
    // Memory allocation
    float *L = malloc(n_lowA_pr * sizeof(float)); // low part (strictly) of L
    float *d = malloc(n_row * sizeof(float));     // L's diagonal

    // * -------------- Cholesky decomposition (A = LL^T) -------------- *
    ltlpr_(&n_row, profile, A_pr, d, &(A_pr[n_row]), L, &eps); // Cholesky decomposition
    /// @note Cholesky decomposition keeps the profile -> we use the same array : profile for L

    // * -------------- Solve descent triangular problem : L^T*v = 2nd_mb -------------- *
    float *v = malloc(n_row * sizeof(float)); // Memory allocation of the temporary solution
    rsprl_(&n_row, profile, d, L, sec_mb_BC, v);

    // * -------------- Solve ascent triangular problem : L*u = v -------------- *
    rspru_(&n_row, profile, d, L, v, u);
    // Memory clean
    free(L);
    free(d);
    free(v);
}