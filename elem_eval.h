#ifndef FEM_PROJECT_ELEM_EVAL_H
#define FEM_PROJECT_ELEM_EVAL_H

// !! In these functions, tab allocation is supposed to be done outside the function with correct dimension

/// Hypothesis
#define DIM 2       // We work in a 2-dimensional space
#define N_NOD_EDG 2 // (we suppose we approximate the element geometry at the order 1)




void evalFbase(int t, const float *x_hat, float *w_hat_x_hat);
/*
    evaluate every base function w_i i in [1,q] for a given point x_hat in the coordinates of the reference element
*/


void evalDFbase(int t, const float *x_hat, float **Dw_hat);
/*
    evaluate every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
    The dimension of the result (Dw_hat) is q x dim (q is the number of nodes here element degree is 1 so q = p = number of vertices)
    Here dim = 2 (2D) except for the segment where dim = 1
*/


float inv_2x2(float **M, float **M_inv);
/*
    Inverse a 2 by 2 matrix, store the result into another matrix and return the determinant
*/


int quad_order(int t);
/*
    Return the quadrature order corresponding to the element type.
    Needed before calling wp_quad to allocate coordinates x_hat and weight tabs
*/

void wp_quad(int t, float **x_quad_hat, float *weight);
/*
    Compute quadrature points and weights corresponding to the element type
*/


void transFK(int q, const float **a_K, const float *w_hat_x_hat, float *Fk_x_hat);
/*
    evaluate transformation Fk of a point of the reference element x_hat
    - q : number of nodes (=p because we take element to 1st degree)
    - a_K contains the coordinates of  element K 's nodes
    - w_hat_x_hat is the vector that contains pre-compute base function value of the point x_hat
    - Fk_x_hat is the result
*/


void jacobFK(int p, int d, const float **a_K, const float **Dw_hat_x_hat, float **jacob_Fk_x_hat);
/*
    evaluate the jacobian matrix of the transformation Fk of a point of the reference element x_hat
    - p               : number of nodes (=p because we take element to 1st degree)
    - d               : dimension 2 (quadrangle/triangle) or 1 (edge)
    - a_K             : contains the coordinates of  element K 's nodes
    - Dw_hat_x_hat    : the vector that contains pre-compute base function value of the point x_hat
    - jacob_Fk_x_hat  : result
*/


void vertices_Edge(int edge_nb, int t, int* local_nodes);
/*
    Give the local number of the nodes at each side of the edge
    (we suppose we approximate the element geometry at the order 1)
    - edge_nb     : local number of the edge for which we want the end nodes
    - t           : type of the element
    - local_nodes : result (array of size 2)
*/


void selectPts(int n_pts, const int pts_nb[], float *coordSet[], float *select_coord[]);
/*
    Select n_pts in the coordSet and copy them into select_coord.
        coordSet and sel_coord are arrays of pointers to array that contains coordinates of a point
        ! Warning, it's only a copy of the pointers to the coordinates, not a deep copy !
    - n_pts        : nb of pts to copy
    - pts_nb       : numbers of points to copy. Can be the global number (resp local number)
    - coordSet     : set of coordinates of points in which we're going to extract some points. Can be all nodes coordinates (resp nodes coordinates of 1 element)
    - select_coord : results
*/

void q_contrib_gW(int n_nod_elem, const float *w_x, float diff, float g_x, float *sum_contrib);
/*
    Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)w_i(x_k) for
    i in [1,n_node_elem]
[IN]
    - n_node_elem  : number of nodes per element
    - *w_x         : w_i(x_k) values of bases function in the quadrature point (length: n_node_elem)
    - diff         : (dx)*|det_JF|  differential element (dx) multiplied by weight for the quadrature formula
    - g_x          : g(x_k) value of g function (of the shape) in the quadrature point x_k
[OUT]
    - sum_contrib  : vector (length nb quad pts) of sum of the other contribution that have already been computed
*/

void q_contrib_gWW(int n_nod_elem, const float *w_x, float diff, float g_x, float **sum_contrib);
/*
    Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)w_i(x_k)w_j(x_k)
    for (i,j) in [1,node_elem]*[1,node_elem];
[IN]
    - n_node_elem  : number of nodes per element
    - *w_x         : w_i/j(x_k) values of bases function in the quadrature point (length: n_node_elem)
    - diff         : (dx)*|det_JF| -> differential element (dx) multiplied by weight for the quadrature formula
    - g_x          : g(x_k) value of g function (of the shape) in the quadrature point x_k
[OUT]
    - sum_contrib  : matrix (size: n_node_elem * n_node_elem) of sum of the other contribution that have already been computed
*/

void q_contrib_gdWdW(int n_nod_elem, const float **Dw_x, float diff, const float **g_x, float **sum_contrib);
/*
    Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)Dw_i(x_k)Dw_j(x_k)
    for (i,j) in [1,node_elem]*[1,node_elem];
[IN]
    - n_node_elem  : number of nodes per element
    - *Dw_x        : Dw_i/j(x_k) -> values of bases function gradient in the quadrature point (size: n_node_elem*n_node_elem)
    - diff         : (dx)*|det_JF| -> differential element (dx) multiplied by weight for the quadrature formula
    - g_x          : g(x_k) value of g function (of the shape) in the quadrature point x_k
[OUT]
    - sum_contrib  : matrix (size: n_node_elem * n_node_elem) of sum of the other contribution that have already been computed
*/


void eval_K(
        int ref_interior,
        const int *ref_Dh,
        const int *ref_Dnh,
        const int *ref_NF,
        int n_Dh,
        int n_Dnh,
        int n_NF,
        int t,
        int n_nod_elem,
        const float **a_K,
        int n_edg_elem,
        const int *ref_edg_K,
        float **A_K,
        float *l_K,
        int *nodes_D,
        float *uD_aK
);

#endif //FEM_PROJECT_ELEM_EVAL_H
