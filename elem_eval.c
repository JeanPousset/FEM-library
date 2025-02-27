#include "elem_eval.h"
#include <stdio.h>

// !! In these functions, tab allocation is supposed to be done outside the function with correct dimension


// evaluate every base function w_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalFbase(int t, const float *x_hat, float *w_hat)
{
    switch (t) {
        case 1: // Quadrangle
            w_hat[0] = x_hat[0]*(1-x_hat[1]);
            w_hat[1] = x_hat[0]*x_hat[1];
            w_hat[2] = (1-x_hat[0])*x_hat[1];
            w_hat[3] = (1-x_hat[0])*(1-x_hat[1]);
            break;
        case 2: // Triangle
            w_hat[0] = x_hat[0];
            w_hat[1] = x_hat[1];
            w_hat[2] = (1-x_hat[0])*(1-x_hat[1]);
            break;
        case 3: // Segment
            w_hat[0] = x_hat[0];
            w_hat[1] = (1-x_hat[0]);
            break;
        default:
            printf("Unknown type t = %d for 'evalFbase' function",t);
            perror("Error : Unvalid parameter");
    }
}

// evaluate every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalDFbase(int t, const float *x_hat, float **Dw_hat)
{
    switch (t) {
        case 1: // Quadrangle
            Dw_hat[0][1] = 1-x_hat[1];    Dw_hat[0][1] = x_hat[0]*(1-x_hat[1]);
            Dw_hat[1][0] = x_hat[1];      Dw_hat[1][1] = x_hat[0];
            Dw_hat[2][0] = -x_hat[1];     Dw_hat[2][1] = 1-x_hat[0];
            Dw_hat[3][0] = x_hat[1]-1;    Dw_hat[3][1] = x_hat[0]-1;
            break;
        case 2: // Triangle
            Dw_hat[0][0] = 1;             Dw_hat[0][1] = 0;
            Dw_hat[1][0] = 0;             Dw_hat[0][1] = 0;
            Dw_hat[2][0] = x_hat[1]-1;    Dw_hat[0][1] = x_hat[0]-1;
            break;
        case 3: // Segment (! Matrix of size 2*1)
            Dw_hat[0][0] = -1;
            Dw_hat[1][0] = -1;
            break;
        default:
            printf("Unknown type t = %d for 'evalFbase' function",t);
            perror("Error : Unvalid parameter");
    }
}

// Inverse a 2 by 2 matrix, store the result into another matrix and return the determinant
float inv_2x2(float **M, float **M_inv)
{
    float det = M[0][0]*M[1][1]-M[0][1]*M[1][0]; // compute determinant
    // Throw an error if the matrix isn't invertible
    if (det == 0.0){
        perror("Matrix is not invertible (inv_2x2) call");
        return(1);
    }else{
        float inv_det = 1/det;
        M_inv[0][0] = inv_det*M[1][1];       M_inv[0][1] = -inv_det*M[0][1];
        M_inv[1][0] = -inv_det*M[1][0];  M_inv[1][1] = inv_det*M[0][0];
        return det;
    }
}

// Return the quadrature order corresponding to the element type
int quad_order(int t)
{
    switch (t) {
        case 1: // Quadrangle
            return 9;
        case 2: // Triangle
        case 3: // Segment (fallthrought)
            return 3;
        default:
            printf("Unknown type t = %d for 'quad_order' function",t);
            perror("Error : Unvalid parameter");
            return 0;
    }
}

// Compute quadrature points and weights corresponding to the element type
void wp_quad(int t, float **x_quad, float *weight)
{
    int i;
    switch (t) {

        float weight_1;
        float weight_2;
        case 1: // Quadrangle
            // weights (pre-compute to avoid useless identical division)
            weight_1 = 1.0 / 36;
            weight_2 = 1.0/9;
            for (i=0; i<4; i++) weight[i] = weight_1;
            for (i=4; i<8; i++) weight[i] = weight_2;
            weight[8] = 4*weight_1;
            // quadrature points (x_i_hat)
            x_quad[0][0] = 1;     x_quad[0][1] = 0;   // x1
            x_quad[1][0] = 1;     x_quad[1][1] = 1;   // x2
            x_quad[2][0] = 0;     x_quad[2][1] = 1;   // x3
            x_quad[3][0] = 0;     x_quad[3][1] = 0;   // x4
            x_quad[4][0] = 1;     x_quad[4][1] = 0.5; // x5
            x_quad[5][0] = 0.5;   x_quad[5][1] = 1;   // x6
            x_quad[6][0] = 0;     x_quad[6][1] = 0.5; // x7
            x_quad[7][0] = 0.5;   x_quad[7][1] = 0;   // x8
            x_quad[8][0] = 0.5;   x_quad[8][1] = 0.5; // x9
            break;

        case 2: // Triangle
            // weights
            weight_1 = 1.0/6;
            for (i=0; i<3; i++) weight[i] = weight_1;
            // quadrature points (x_i_hat)
            x_quad[0][0] = 0.5;   x_quad[0][1] = 0.5;  // x1
            x_quad[1][0] = 0;     x_quad[1][1] = 0.5;  // x2
            x_quad[2][0] = 0.5;   x_quad[2][1] = 0;    // x3
            break;

        case 3: // Segment
            // weights
            weight_1 = 1.0/6;
            weight[0] = weight_1;
            weight[0] = weight_1;
            weight[0] = 4*weight_1;
            // quadrature points (x_i_hat)
            x_quad[0][0] = 1;    // x1
            x_quad[1][0] = 0;    // x2
            x_quad[2][0] = 0.5;  // x3
            break;

        default:
            printf("Unknown type t = %d for 'quad_order' function",t);
            perror("Error : Unvalid parameter");
    }
}

// evaluate transformation Fk of a point of the reference element x_hat
void transFk(int q, float **a_K, const float *w_hat_x_hat, float *Fk_x_hat)
{
    Fk_x_hat[0] = 0;
    Fk_x_hat[1] = 0;
    for (int i=0; i<q; i++){
        Fk_x_hat[0] += a_K[i][0] * w_hat_x_hat[i];
        Fk_x_hat[0] += a_K[i][1] * w_hat_x_hat[i];
    }
}

// evaluate the jacobian matrix of the transformation Fk of a point of the reference element x_hat
void jacobFk(int p, int d, float **a_K, float **Dw_x_hat, float **jacob_Fk_x_hat)
{
    int i;
    if (d==2) { // quadrangle, triangle
        // initialisation
        jacob_Fk_x_hat[0][0] = 0; jacob_Fk_x_hat[0][1] = 0;
        jacob_Fk_x_hat[1][0] = 0; jacob_Fk_x_hat[1][1] = 0;
        for(i=0; i<p; i++){
            jacob_Fk_x_hat[0][0] += a_K[i][0] * Dw_x_hat[i][0]; jacob_Fk_x_hat[0][1] += a_K[i][0] * Dw_x_hat[i][1];
            jacob_Fk_x_hat[1][0] += a_K[i][1] * Dw_x_hat[i][0]; jacob_Fk_x_hat[1][1] += a_K[i][1] * Dw_x_hat[i][1];
        }
    } else if (d==1){ // segment
        jacob_Fk_x_hat[0][0] = 0;
        jacob_Fk_x_hat[1][0] = 0;
        for(i=0; i<p; i++){
            jacob_Fk_x_hat[0][0] += a_K[i][0] * Dw_x_hat[i][0];
            jacob_Fk_x_hat[1][0] += a_K[i][1] * Dw_x_hat[i][0];
        }
    } else {
        printf("Unknown dimension d = %d for 'jacobFk' function \n",d);
        perror("Error : Unvalid parameter");
    }
}

// give the local number of the nodes at each side of the edge
void vertices_Edge(int edge_nb, int t, int* local_nodes)
{
    switch(t){
        case 1: // quadrangle
            local_nodes[0] = edge_nb ;
            local_nodes[1] = edge_nb % 4 + 1;
            break;
        case 2: // triangle
            local_nodes[0] = edge_nb ;
            local_nodes[1] = edge_nb % 3 + 1;
            break;
        default:
            printf("Unknown type t = %d for 'quad_order' function ."
                   "\n Carefull this function can not be called for a segment",t);
            perror("Error : Unvalid parameter");
    }
}

// select n_pts in the coordSet and copy them into sel_coord
void selectPts(int n_pts, const int pts_nb[], float *coordSet[], float *select_coord[])
{
    for(int i=0; i<n_pts; i++) select_coord[i] = coordSet[pts_nb[i]-1]; // -1 because numbers begins at 1 in pts_nb
}

