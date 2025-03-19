/**
 * @file elem_eval.h
 * @brief Header for functions that evaluate discretization of variational formulation : A_K_ij and l_K_i for the element K
 *
 * @warning In these functions, tab allocation is supposed to be done outside the function with correct dimension.
 * @warning We work in a 2-dimensional space.
 * @warning We suppose we approximate the element geometry at the degree 1 (N_NOD_EDG : 2) if degree > 1 then it becomes
 *          a variable that is equal to degree + 1
 * @warning We suppose dynamical tab of sum of contribution are initialized to 0.
 */


#include "elem_eval.h"
#include "tab_mngmt.h"
#include "problem_functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/// Evaluates every base function w_i i in [1,q] for a given point x_hat in the coordinates of the reference element
void evalFbase(int t, const float *x_hat, float *w_hat_x_hat) {
    switch (t) {
        case 1: // Quadrangle
            w_hat_x_hat[0] = x_hat[0] * (1 - x_hat[1]);
            w_hat_x_hat[1] = x_hat[0] * x_hat[1];
            w_hat_x_hat[2] = (1 - x_hat[0]) * x_hat[1];
            w_hat_x_hat[3] = (1 - x_hat[0]) * (1 - x_hat[1]);
            break;
        case 2: // Triangle
            w_hat_x_hat[0] = x_hat[0];
            w_hat_x_hat[1] = x_hat[1];
            w_hat_x_hat[2] = 1 - x_hat[0] - x_hat[1];
            break;
        case 3: // Segment
            w_hat_x_hat[0] = x_hat[0];
            w_hat_x_hat[1] = (1 - x_hat[0]);
            break;
        default:
            printf("Unknown type t = %d for 'evalFbase' function", t);
            perror("Error : Unvalid parameter");
    }
}

/// Evaluates every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalDFbase(int t, const float *x_hat, float **Dw_hat_x_hat) {
    switch (t) {
        case 1: // Quadrangle
            Dw_hat_x_hat[0][0] = 1 - x_hat[1];
            Dw_hat_x_hat[0][1] = x_hat[0] * (1 - x_hat[1]);
            Dw_hat_x_hat[1][0] = x_hat[1];
            Dw_hat_x_hat[1][1] = x_hat[0];
            Dw_hat_x_hat[2][0] = -x_hat[1];
            Dw_hat_x_hat[2][1] = 1 - x_hat[0];
            Dw_hat_x_hat[3][0] = x_hat[1] - 1;
            Dw_hat_x_hat[3][1] = x_hat[0] - 1;
            break;
        case 2: // Triangle
            Dw_hat_x_hat[0][0] = 1;
            Dw_hat_x_hat[0][1] = 0;
            Dw_hat_x_hat[1][0] = 0;
            Dw_hat_x_hat[1][1] = 1;
            Dw_hat_x_hat[2][0] = -1;
            Dw_hat_x_hat[2][1] = -1;
            break;
        case 3: // Segment (! Matrix of size 2*1)
            Dw_hat_x_hat[0][0] = 1;
            Dw_hat_x_hat[1][0] = -1;
            break;
        default:
            printf("Unknown type t = %d for 'evalFbase' function", t);
            perror("Error : Unvalid parameter");
    }
}

/// Inverse a 2 by 2 matrix, store the result into another matrix and return the determinant
float inv_2x2(float **M, float **M_inv) {
    float det = M[0][0] * M[1][1] - M[0][1] * M[1][0]; // compute determinant
    if (det < 10e-9) { // Throw an error if the matrix isn't invertible
        printf("determinant value that pose a problem : %f\n", det);
        perror("Matrix is almost singular or not invertible (inv_2x2) call");
        return (1);
    } else {
        float inv_det = 1 / det;
        M_inv[0][0] = inv_det * M[1][1];
        M_inv[0][1] = -inv_det * M[0][1];
        M_inv[1][0] = -inv_det * M[1][0];
        M_inv[1][1] = inv_det * M[0][0];
        return det;
    }
}

/// Return the quadrature order corresponding to the element type
int quad_order(int t) {
    switch (t) {
        case 1: // Quadrangle
            return 9;
        case 2: // Triangle
        case 3: // Segment (fallthrought)
            return 3;
        default:
            printf("Unknown type t = %d for 'quad_order' function", t);
            perror("Error : Unvalid parameter");
            return 0;
    }
}

/// Compute quadrature points and weights corresponding to the element type
void wp_quad(int t, float **x_quad_hat, float *weights) {
    int i;
    float weight_1;
    float weight_2;
    switch (t) {

        case 1: // Quadrangle
            // weights (pre-compute to avoid useless identical division)
            weight_1 = 1.0 / 36;
            weight_2 = 1.0 / 9;
            for (i = 0; i < 4; i++) weights[i] = weight_1;
            for (i = 4; i < 8; i++) weights[i] = weight_2;
            weights[8] = 4 * weight_1;
            // quadrature points (x_i_hat)
            x_quad_hat[0][0] = 1;
            x_quad_hat[0][1] = 0;   // x1
            x_quad_hat[1][0] = 1;
            x_quad_hat[1][1] = 1;   // x2
            x_quad_hat[2][0] = 0;
            x_quad_hat[2][1] = 1;   // x3
            x_quad_hat[3][0] = 0;
            x_quad_hat[3][1] = 0;   // x4
            x_quad_hat[4][0] = 1;
            x_quad_hat[4][1] = 0.5; // x5
            x_quad_hat[5][0] = 0.5;
            x_quad_hat[5][1] = 1;   // x6
            x_quad_hat[6][0] = 0;
            x_quad_hat[6][1] = 0.5; // x7
            x_quad_hat[7][0] = 0.5;
            x_quad_hat[7][1] = 0;   // x8
            x_quad_hat[8][0] = 0.5;
            x_quad_hat[8][1] = 0.5; // x9
            break;

        case 2: // Triangle
            // weights
            weight_1 = 1.0 / 6;
            for (i = 0; i < 3; i++) weights[i] = weight_1;
            // quadrature points (x_i_hat)
            x_quad_hat[0][0] = 0.5;
            x_quad_hat[0][1] = 0.5;  // x1
            x_quad_hat[1][0] = 0;
            x_quad_hat[1][1] = 0.5;  // x2
            x_quad_hat[2][0] = 0.5;
            x_quad_hat[2][1] = 0;    // x3
            break;

        case 3: // Segment
            // weights
            weight_1 = 1.0 / 6;
            weights[0] = weight_1;
            weights[1] = weight_1;
            weights[2] = 4 * weight_1;
            // quadrature points (x_i_hat)
            x_quad_hat[0][0] = 1;    // x1
            x_quad_hat[1][0] = 0;    // x2
            x_quad_hat[2][0] = 0.5;  // x3
            break;

        default:
            printf("Unknown type t = %d for 'quad_order' function", t);
            perror("Error : Unvalid parameter");
    }
}

/// Evaluates transformation Fk of a point of the reference element x_hat
void transFK(int n_nod_elem, const float **a_K, const float *w_hat_x_hat, float *Fk_x_hat) {
    Fk_x_hat[0] = 0;
    Fk_x_hat[1] = 0;
    for (int i = 0; i < n_nod_elem; i++) {
        //printf("a_K %f\n",a_K[i][0]);
        //printf("w_hat %f\n",w_hat_x_hat[i]);
        Fk_x_hat[0] += a_K[i][0] * w_hat_x_hat[i];
        Fk_x_hat[0] += a_K[i][1] * w_hat_x_hat[i];
    }
}

/// Evaluates the jacobian matrix of the transformation Fk of a point of the reference element x_hat
void jacobFK(int n_nod_elem, int d, const float **a_K, const float **Dw_hat_x_hat, float **jacob_Fk_x_hat) {
    int i;
    if (d == 2) { // quadrangle, triangle
        // initialisation
        jacob_Fk_x_hat[0][0] = 0;
        jacob_Fk_x_hat[0][1] = 0;
        jacob_Fk_x_hat[1][0] = 0;
        jacob_Fk_x_hat[1][1] = 0;
        for (i = 0; i < n_nod_elem; i++) {
            jacob_Fk_x_hat[0][0] += a_K[i][0] * Dw_hat_x_hat[i][0];
            jacob_Fk_x_hat[0][1] += a_K[i][0] * Dw_hat_x_hat[i][1];
            jacob_Fk_x_hat[1][0] += a_K[i][1] * Dw_hat_x_hat[i][0];
            jacob_Fk_x_hat[1][1] += a_K[i][1] * Dw_hat_x_hat[i][1];
        }
    } else if (d == 1) { // segment
        jacob_Fk_x_hat[0][0] = 0;
        jacob_Fk_x_hat[1][0] = 0;
        for (i = 0; i < n_nod_elem; i++) {
            jacob_Fk_x_hat[0][0] += a_K[i][0] * Dw_hat_x_hat[i][0];
            jacob_Fk_x_hat[1][0] += a_K[i][1] * Dw_hat_x_hat[i][0];
        }
    } else {
        printf("Unknown dimension d = %d for 'jacobFK' function \n", d);
        perror("Error : Unvalid parameter");
    }
}

/// Give the local number of the nodes at each side of the edge (we suppose we approximate the element geometry at the order 1)
void vertices_Edge(int edge_nb, int t, int *local_nodes) {
    switch (t) {
        case 1: // quadrangle
            local_nodes[0] = edge_nb;
            local_nodes[1] = edge_nb % 4 + 1;
            break;
        case 2: // triangle
            local_nodes[0] = edge_nb;
            local_nodes[1] = edge_nb % 3 + 1;
            break;
        default:
            printf("Unknown type t = %d for 'quad_order' function ."
                   "\n Carefull this function can not be called for a segment", t);
            perror("Error : Unvalid parameter");
    }
}

/// select n_pts in the coordSet and copy them into sel_coord
void selectPts(int n_pts, const int pts_nb[], float *coordSet[], float *select_coord[]) {
    for (int i = 0; i < n_pts; i++)
        select_coord[i] = coordSet[pts_nb[i]]; // (10) -1 or not because numbers begins at 1 in pts_nb
    // currently we dont put pts nb but pts indices
}

/// Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)w_i(x_k)
void q_contrib_gW(int n_nod_elem, const float *w_x, float diff, float g_x, float *sum_contrib) {
    for (int i = 0; i < n_nod_elem; i++) sum_contrib[i] += diff * g_x * w_x[i];
}

/// Compute quadrature contribution of the quadrature point x for each base function w_i for the integral of the shape w_i*g
void q_contrib_gWW(int n_nod_elem, const float *w_x, float diff, float g_x, float **sum_contrib) {
    int i, j;
    float base, contrib;
    base = diff * g_x;
    for (i = 0; i < n_nod_elem; i++) {
        contrib = base * w_x[i];
        for (j = 0; j < n_nod_elem; j++) {
            sum_contrib[i][j] += contrib * w_x[j];
        }
    }
}

/**Compute quadrature contribution of the quadrature point x for each base function w_i for the integral of the shape
    sum (alpha,beta = 1 -> 2) g_{alpha,beta}(x)*d_alpha(w_i)*d_beta(w_j) */
void q_contrib_gdWdW(int n_nod_elem, const float **Dw_x, float diff, const float **g_x, float **sum_contrib) {
    int i, j, alpha, beta;
    float base, contrib;

    for (alpha = 0; alpha < 2; alpha++) {
        for (beta = 0; beta < 2; beta++) {
            base = diff * g_x[alpha][beta];
            for (i = 0; i < n_nod_elem; i++) {
                contrib = base * Dw_x[i][alpha];
                for (j = 0; j < n_nod_elem; j++) {
                    sum_contrib[i][j] += contrib * Dw_x[j][beta];
                }
            }
        }
    }
}

/// Evaluates the integrals on the element K with quadrature formula */
void intElem(int t, int n_quad_pts, int n_nod_elem, const float **a_K, const float **x_hat_quad, const float *weights,
             float **A_K_elem, float *l_K_elem) {
    int k, i, j;
    float detJF, diff;

    // Memory allocation for object reset in each quadrature point
    float *x_quad = malloc(2 * sizeof(float));
    float *w_hat_x_hat = malloc(n_nod_elem * sizeof(float));
    float **Dw_hat_x_hat = matF_alloc(n_nod_elem, DIM);
    float **Dw_x_quad = matF_alloc(n_nod_elem, DIM);
    float **JFk_x_hat = matF_alloc(DIM, DIM);
    float **inv_JFk_x_hat = matF_alloc(DIM, DIM);
    float **a_alpha_beta = matF_alloc(DIM, DIM);

    for (k = 0; k < n_quad_pts; k++) {

        // evaluation of bases function et gradient of bases functions in the quadrature point x_hat_quad
        evalFbase(t, x_hat_quad[k], w_hat_x_hat);
        evalDFbase(t, x_hat_quad[k], Dw_hat_x_hat);

        // Get the quadrature points coordinates x_quad_k with F_K
        transFK(n_nod_elem, (const float **) a_K, w_hat_x_hat, x_quad);

        // Evaluation Dw(x_quad)/Dx_alpha, alpha in {1,2} and detJFK
        jacobFK(n_nod_elem, DIM, (const float **) a_K, (const float **) Dw_hat_x_hat, JFk_x_hat);
        detJF = inv_2x2(JFk_x_hat, inv_JFk_x_hat);
        diff = fabsf(detJF) * weights[k]; // differential elem * weight for the quadrature formula
        for (j = 0; j < 2; j++) {
            for (i = 0; i < n_nod_elem; i++) {
                Dw_x_quad[i][j] = Dw_hat_x_hat[i][0] * inv_JFk_x_hat[0][j] + Dw_hat_x_hat[i][1] * inv_JFk_x_hat[1][j];
            }
        }
        // pre-compute a_{alpha,beta}(x_quad)
        a_alpha_beta[0][0] = a11(x_quad);
        a_alpha_beta[0][1] = a12(x_quad);
        a_alpha_beta[1][0] = a12(x_quad);
        a_alpha_beta[1][1] = a22(x_quad);

        // Contribution of the quadrature point for the quadrature formula
        q_contrib_gW(n_nod_elem, w_hat_x_hat, diff, f_OMEGA(x_quad), l_K_elem);
        q_contrib_gWW(n_nod_elem, (const float *) w_hat_x_hat, diff, a00(x_quad), A_K_elem);
        q_contrib_gdWdW(n_nod_elem, (const float **) Dw_hat_x_hat, diff, (const float **) a_alpha_beta, A_K_elem);
    }
    // free memory
    free(x_quad);
    free(w_hat_x_hat);
    free_mat(Dw_hat_x_hat);
    free_mat(Dw_x_quad);
    free_mat(JFk_x_hat);
    free_mat(inv_JFk_x_hat);
    free_mat(a_alpha_beta);
}


/// Evaluates the integrals on the edges of element K (Kp = K prime) with quadrature formula
void
intEdge(int n_quad_pts, const float **a_Kp, const float **x_hat_quad, const float *weights, float **A_K_edg,
        float *l_K_edg
) {

    int k;
    float Kp_length, diff;

    // Memory allocation for object reset in each quadrature point
    float *x_quad = malloc(2 * sizeof(float));
    float *w_hat_x_hat = malloc(N_NOD_EDG * sizeof(float));
    float **Dw_hat_x_hat = matF_alloc(N_NOD_EDG, DIM);
    float **JFk_x_hat = matF_alloc(DIM, 1);

    for (k = 0; k < n_quad_pts; k++) {

        // evaluation of bases function et gradient of bases functions in the quadrature point x_hat_quad
        evalFbase(SEG_TYPE, x_hat_quad[k], w_hat_x_hat);
        evalDFbase(SEG_TYPE, x_hat_quad[k], Dw_hat_x_hat);

        transFK(N_NOD_EDG, a_Kp, w_hat_x_hat, x_quad); // Get the quadrature points coordinates x_quad_k
        jacobFK(N_NOD_EDG, 1, a_Kp, (const float **) Dw_hat_x_hat, JFk_x_hat);

        // evaluation of segment length to get the differential term in the quadrature formula
        Kp_length = sqrtf(JFk_x_hat[0][0] * JFk_x_hat[0][0] + JFk_x_hat[1][0] * JFk_x_hat[1][0]);
        diff = Kp_length * weights[k]; // differential coef * weight for the quadrature formula

        // Contribution of the quadrature point for the quadrature formula
        q_contrib_gW(N_NOD_EDG, (const float *) w_hat_x_hat, diff, fN(x_quad), l_K_edg);
        q_contrib_gWW(N_NOD_EDG, (const float *) w_hat_x_hat, diff, bN(x_quad), A_K_edg);
    }

    // free memory
    free(x_quad);
    free(w_hat_x_hat);
    free_mat(Dw_hat_x_hat);
    free_mat(JFk_x_hat);
}


/// Evaluates A_K_ij and second member l_K_i form the discretization of the variational formulation for the element K
void
eval_K(int ref_interior, const int *ref_Dh, const int *ref_Dnh, const int *ref_NF, int n_Dh, int n_Dnh, int n_NF, int t,
       int n_nod_elem, const float **a_K, int n_edg_elem, const int *ref_edg_K, float **A_K, float *l_K, int *nodes_D,
       float *uD_aK
) {
    int i, j, k, l;
    int flag_categorized;     // boolean to know if the edge has already been categorized
    int n_quad_pts_elem, n_quad_pts_edg;

    // ----------- Integral on the element K -----------

    n_quad_pts_elem = quad_order(t);
    // Memory allocation
    float **x_quad_hat_elem = matF_alloc(n_quad_pts_elem, 2);
    float *weights_elem = malloc(n_quad_pts_elem * sizeof(float));
    float **A_K_elem = matF_alloc0(n_nod_elem, n_nod_elem);  // Every value set to 0
    float *l_K_elem = calloc(n_nod_elem, sizeof(float));     // Every value set to 0

    wp_quad(t, x_quad_hat_elem, weights_elem); // Evaluates weights and quadrature pts for the quadrature on the element
    intElem(t, n_quad_pts_elem, n_nod_elem, (const float **) a_K, (const float **) x_quad_hat_elem, weights_elem,
            A_K_elem, l_K_elem);

    // Update results :
    for (i = 0; i < n_nod_elem; i++) {
        l_K[i] += l_K_elem[i];
        for (j = 0; j < n_nod_elem; j++) A_K[i][j] += A_K_elem[i][j];
    }
    // free Memory
    free_mat(x_quad_hat_elem);
    free(weights_elem);
    free_mat(A_K_elem);
    free(l_K_elem);


    // ----------- Integral on a potential edge element K

    // Memory update (clean and new)
    n_quad_pts_edg = quad_order(SEG_TYPE);
    float **x_quad_hat_edg = matF_alloc(n_quad_pts_edg, DIM);
    float *weights_edg = malloc(n_quad_pts_edg * sizeof(float));
    int *edg_nod_ind = malloc(N_NOD_EDG * sizeof(int));
    float **A_K_edg = matF_alloc(N_NOD_EDG, N_NOD_EDG);     // Values to be re-initialize to 0 for each Neumann edge
    float *l_K_edg = malloc(N_NOD_EDG * sizeof(float));     // Values to be re-initialize to 0 for each Neumann edge
    // Special case for this one : it's a dynamic tab of float pointers
    float **edg_nodes_coords = malloc(N_NOD_EDG * sizeof(float *));
    // (5) -> Is it a good way to handle it ? Because if I don't do this they will be a pb when memory clean because on the simple copy of selectPts

    for (i = 0; i < n_edg_elem; i++) {

        if (ref_edg_K[i] == ref_interior) { flag_categorized = 1; } // Most current case : edge is inside the domain
        else {
            flag_categorized = 0;                   // edge isn't categorized yet
            vertices_Edge(i + 1, t, edg_nod_ind);
            for (j = 0; j < N_NOD_EDG; j++) edg_nod_ind[j]--;
            // We retrieve indices (2) ? Shouldn't do it directly in vertices edges and change the function's name

            for (j = 0; j < n_Dh & flag_categorized == 0; j++) {      // Dirichlet homogeneous edge
                if (ref_edg_K[i] == ref_Dh[j]) {
                    for (k = 0; k < N_NOD_EDG; k++) nodes_D[edg_nod_ind[k]] = 0;
                    flag_categorized = 1;
                }
            }
            for (j = 0; j < n_Dnh & flag_categorized == 0; j++) {     // Dirichlet non-homogeneous edge
                if (ref_edg_K[i] == ref_Dnh[j]) {
                    for (k = 0; k < N_NOD_EDG; k++) {
                        nodes_D[edg_nod_ind[k]] = -1;
                        uD_aK[edg_nod_ind[k]] = uD(a_K[edg_nod_ind[k]]);
                    }
                    flag_categorized = 1;
                }
            }
            for (j = 0; j < n_NF & flag_categorized == 0; j++) {      // Neumann or Fourier edge
                if (ref_edg_K[i] == ref_NF[j]) {
                    selectPts(N_NOD_EDG, edg_nod_ind, (float **) a_K, edg_nodes_coords);/// @warning lost "const" on a_K
                    wp_quad(SEG_TYPE, x_quad_hat_edg, weights_edg);

                    // Re-intialize A_K_edg and l_K_edg to 0;
                    for (k = 0; k < N_NOD_EDG; k++) {
                        l_K_edg[k] = 0;
                        for (l = 0; l < N_NOD_EDG; l++) A_K_edg[k][l] = 0;
                    }

                    // Update results (we have the contribution for each node at the end of the edge (2 nodes)
                    intEdge(n_quad_pts_edg, (const float **) edg_nodes_coords, (const float **) x_quad_hat_edg,
                            weights_edg, A_K_edg, l_K_edg);
                    for (k = 0; k < N_NOD_EDG; k++) {
                        l_K[edg_nod_ind[k]] += l_K_edg[k];
                        for (l = 0; l < N_NOD_EDG; l++) A_K[edg_nod_ind[k]][edg_nod_ind[l]] += A_K_edg[k][l];
                    }
                    flag_categorized = 1;
                }
            }
        }
    }

    // free memory
    free_mat(x_quad_hat_edg);
    free(l_K_edg);
    free(edg_nod_ind);
    free(edg_nodes_coords); // (5) Here we only free an dynamic array of float pointers -> different from a matrix
    free_mat(A_K_edg);
    free(weights_edg);
}
