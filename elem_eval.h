#ifndef FEM_PROJECT_ELEM_EVAL_H
#define FEM_PROJECT_ELEM_EVAL_H

// !! In these functions, tab allocation is supposed to be done outside the function with correct dimension


/*
    evalute every base function w_i i in [1,q] for a given point x_hat in the coordinates of the reference element
*/
void evalFbase(int t, float *w_hat_i, const float *x_hat);


/*
    evalute every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
*/
void evalDFbase(int t, float **Dw_hat, const float *x_hat);

/*
    Inverse a 2 by 2 matrix, store the result into an other matrix and return the determinant
*/
float inv_2x2(float **M, float **M_inv);

/*
    Return the quadrature order corresponding to the element type.
    Needed before calling wp_quad to allocate coordinates x_hat and weight tabs
*/
int quad_order(int t);

/*
    Compute quadrature points and weights corresponding to the element type
*/
void wp_quad(int t, float **x_quad, float *weight);


#endif //FEM_PROJECT_ELEM_EVAL_H
