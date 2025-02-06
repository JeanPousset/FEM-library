#include "elem_eval.h"
#include <stdio.h>

// !! In these functions, tab allocation is supposed to be done outside the function with correct dimension


// evalute every base function w_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalFbase(int t, float *w_hat, const float *x_hat)
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
            w_hat[0] = x_hat[0]*(1-x_hat[1]);
            w_hat[1] = (1-x_hat[0])*(1-x_hat[1]);
            break;
        default:
            printf("Unknown type t = %d for 'evalFbase' function",t);
            perror("Error : Unvalid parameter");
    }
}

// evalute every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalDFbase(int t, float **Dw_hat, const float *x_hat)
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
        case 3: // Segment
            Dw_hat[0][0] = 1-x_hat[1];    Dw_hat[0][1] = -x_hat[0];
            Dw_hat[1][0] = x_hat[1]-1;    Dw_hat[1][1] = x_hat[0]-1;
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
    M_inv[0][0] = M[1][1]/det;   M_inv[0][1] = -M[0][1]/det;
    M_inv[1][0] = -M[1][0]/det;  M_inv[1][1] = M[0][0]/det;
    return det;
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

        case 1: // Quadrangle
            // weights
            for (i=0; i<4; i++) weight[i] = 1.0/36;
            for (i=4; i<8; i++) weight[i] = 1.0/9;
            weight[8] = 4.0/9;
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
            for (i=0; i<3; i++) weight[i] = 1.0/6;
            // quadrature points (x_i_hat)
            x_quad[0][0] = 0.5;   x_quad[0][1] = 0.5;  // x1
            x_quad[1][0] = 0;     x_quad[1][1] = 0.5;  // x2
            x_quad[2][0] = 0.5;   x_quad[2][1] = 0;    // x3
            break;

        case 3: // Segment
            // weights
            weight[0] = 1.0/6;
            weight[0] = 1.0/6;
            weight[0] = 2.0/3;
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
