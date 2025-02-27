#ifndef FEM_PROJECT_ELEM_EVAL_H
#define FEM_PROJECT_ELEM_EVAL_H

// !! In these functions, tab allocation is supposed to be done outside the function with correct dimension



void evalFbase(int t, const float *x_hat, float *w_hat);
/*
    evalute every base function w_i i in [1,q] for a given point x_hat in the coordinates of the reference element
*/


void evalDFbase(int t, const float *x_hat, float **Dw_hat);
/*
    evalute every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
    The dimension of the result (Dw_hat) is q x dim (q is the number of nodes here element degree is 1 so q = p = number of vertices)
    Here dim = 2 (2D) except for the segment where dim = 1
*/


float inv_2x2(float **M, float **M_inv);
/*
    Inverse a 2 by 2 matrix, store the result into an other matrix and return the determinant
*/


int quad_order(int t);
/*
    Return the quadrature order corresponding to the element type.
    Needed before calling wp_quad to allocate coordinates x_hat and weight tabs
*/

void wp_quad(int t, float **x_quad, float *weight);
/*
    Compute quadrature points and weights corresponding to the element type
*/


void transFk(int q, float **a_K, const float *w_hat_x_hat, float *Fk_x_hat);
/*
    evaluate transformation Fk of a point of the reference element x_hat
    - p : number of nodes (=p because we take element to 1st degree)
    - a_K contains the coordinates of  element K 's nodes
    - w_hat_x_hat is the vector that contains pre-compute base function value of the point x_hat
    - Fk_x_hat is the result
*/


void jacobFk(int p, int d, float **a_K, float **Dw_x_hat, float **jacob_Fk_x_hat);
/*
    evaluate the jacobian matrix of the transformation Fk of a point of the reference element x_hat
    - p              : number of nodes (=p because we take element to 1st degree)
    - a_K            : contains the coordinates of  element K 's nodes
    - w_hat_x_hat    : the vector that contains pre-compute base function value of the point x_hat
    - jacob_Fk_x_hat : result
*/


void vertices_Edge(int edge_nb, int t, int* local_nodes);
/*
    Give the local number of the nodes at each side of the edge
    - edge_nb     : local number of the edge for which we want the end nodes
    - t           : type of the element
    - local_nodes : result (array of size 2)
*/


void selectPts(int n_pts, const int pts_nb[], float *coordSet[], float *select_coord[]);
/*
    Select n_pts in the coordSet and copy them into sel_coord.
        coordSet and sel_coord are arrays of pointers to array that contains coordinates of a point
        ! Warning, it's only a copy of the pointers to the coordinates, not a deep copy !
    - n_pts        : nb of pts to copy
    - pts_nb       : numbers of points to copy. Can be the global number (resp local number)
    - coordSet     : set of coordinates of points in which we're going to extract some points. Can be all nodes coordinates (resp nodes coordinates of 1 element)
    - select_coord : results
*/


#endif //FEM_PROJECT_ELEM_EVAL_H
