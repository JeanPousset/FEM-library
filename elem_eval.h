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
#ifndef FEM_PROJECT_ELEM_EVAL_H
#define FEM_PROJECT_ELEM_EVAL_H

#define DIM 2
#define N_NOD_EDG 2
#define SEG_TYPE 3


/**
 *    @brief Evaluates every base function w_i i in [1,q] for a given point x_hat in the coordinates of the reference element
 *
 *    @param t type of the element (1: quadrangle, 2:triangle, 3:segment)
 *    @param x_hat coordinates in the reference element K_hat of the point to evaluate
 *    @param w_hat_x_hat evaluations of each base function in the point x_hat
*/
void evalFbase(int t, const float *x_hat, float *w_hat_x_hat);


/**
 *    @brief Evaluates every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
 *
 *    @details The dimension of the result (Dw_hat) is q x dim (q is the number of nodes here element degree is 1 so q = p = number of vertices). Here dim = 2 (2D) except for the segment where dim = 1
 *
 *    @param t type of the element (1: quadrangle, 2:triangle, 3:segment)
 *    @param x_hat coordinates in the reference element K_hat of the point to evaluate
 *    @param Dw_hat evaluations of each gradient of base function in the point x_hat
*/
void evalDFbase(int t, const float *x_hat, float **Dw_hat);


/**
 *    @brief Inverse a 2 by 2 matrix, store the result into another matrix and return the determinant
 *
 *    @param M      matrix whose inverse is desired
 *    @param M_inv  inverse of M
*/
float inv_2x2(float **M, float **M_inv);


/**
 *    @brief Returns the quadrature order corresponding to the element type.
 *
 *    @details Needed before calling wp_quad to allocate coordinates x_hat and weight tabs
 *
 *    @param t type of the element (1: quadrangle, 2:triangle, 3:segment)
 *    @return quadrature order corresponding to the type (quadrangle: 9, triangle or segment:3)
*/
int quad_order(int t);


/**
 *      @brief Compute quadrature points and weights corresponding to the element type
 *
 *      @param t           type of the element (1: quadrangle, 2:triangle, 3:segment)
 *      @param x_quad_hat  coordinates of the quadrature points in the reference element K_hat
 *      @param weights     weight associated to each quadrature point
*/
void wp_quad(int t, float **x_quad_hat, float *weights);


/**
 *   @brief Evaluates transformation Fk of a point of the reference element x_hat
 *
 *   @param n_nod_elem  number of nodes (=p because we take element to 1st degree)
 *   @param a_K         contains the coordinates of  element K 's nodes
 *   @param w_hat_x_hat vector that contains pre-compute base function value of the point x_hat
 *   @param Fk_x_hat    result
 */
void transFK(int n_nod_elem, const float **a_K, const float *w_hat_x_hat, float *Fk_x_hat);


/**
 *    @brief Evaluates the jacobian matrix of the transformation Fk of a point of the reference element x_hat
 *
 *    @param n_nod_elem     number of nodes per element (=p because we take element to 1st degree)
 *    @param d              dimension 2 (quadrangle/triangle) or 1 (edge)
 *    @param a_K            contains the coordinates of  element K 's nodes
 *    @param Dw_hat_x_hat   the vector that contains pre-compute base function value of the point x_hat
 *    @param jacob_Fk_x_hat result
 */
void jacobFK(int n_nod_elem, int d, const float **a_K, const float **Dw_hat_x_hat, float **jacob_Fk_x_hat);


/**
 *   @brief Give the local number of the nodes at each side of the edge
 *   @warning we suppose we approximate the element geometry at the order 1
 *
 *   @param edge_nb     local number of the edge for which we want the end nodes
 *   @param t           type of the element
 *   @param local_nodes result (array of size 2)
 */
void vertices_Edge(int edge_nb, int t, int *local_nodes);


/**
 *    @brief Select n_pts in the coordSet and copy them into select_coord.
 *
 *    @warning coordSet and sel_coord are arrays of pointers to array that contains coordinates of a point ! Warning, it's only a copy of the pointers to the coordinates, not a deep copy !
 *
 *    @param n_pts        nb of pts to copy
 *    @param pts_nb       numbers of points to copy. Can be the global number (resp local number)
 *    @param coordSet     set of coordinates of points in which we're going to extract some points. Can be all nodes coordinates (resp nodes coordinates of 1 element)
 *    @param select_coord results
 */
void selectPts(int n_pts, const int pts_nb[], float *coordSet[], float *select_coord[]);


/**
 *    @brief Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)w_i(x_k) for i in [1,n_node_elem]
 *
 *    @param n_nod_elem  number of nodes per element
 *    @param w_x         w_i(x_k) values of bases function in the quadrature point (length: n_node_elem)
 *    @param diff        (dx)*|det_JF|  differential element (dx) multiplied by weight for the quadrature formula
 *    @param g_x         g(x_k) value of g function (of the shape) in the quadrature point x_k
 *    @param sum_contrib vector (length nb quad pts) of sum of the other contribution that have already been computed
 */
void q_contrib_gW(int n_nod_elem, const float *w_x, float diff, float g_x, float *sum_contrib);


/**
 *   @brief Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)w_i(x_k)w_j(x_k) for (i,j) in [1,node_elem]*[1,node_elem];
 *
 *   @param n_nod_elem  number of nodes per element
 *   @param w_x         w_i/j(x_k) values of bases function in the quadrature point (length: n_node_elem)
 *   @param diff         (dx)*|det_JF| -> differential element (dx) multiplied by weight for the quadrature formula
 *   @param g_x          g(x_k) value of g function (of the shape) in the quadrature point x_k
 *   @param sum_contrib  matrix (size: n_node_elem * n_node_elem) of sum of the other contribution that have already been computed
 */
void q_contrib_gWW(int n_nod_elem, const float *w_x, float diff, float g_x, float **sum_contrib);


/**
 *    @brief Compute the contribution of the quadrature points x_k for the integral's quadrature of the shape g(x_k)Dw_i(x_k)Dw_j(x_k) for (i,j) in [1,node_elem]*[1,node_elem];
 *
 *    @param n_nod_elem  number of nodes per element
 *    @param Dw_x        Dw_i/j(x_k) -> values of bases function gradient in the quadrature point (size: n_node_elem*n_node_elem)
 *    @param diff        (dx)*|det_JF| -> differential element (dx) multiplied by weight for the quadrature formula
 *    @param g_x         g(x_k) value of g function (of the shape) in the quadrature point x_k
 *    @param sum_contrib matrix (size: n_node_elem * n_node_elem) of sum of the other contribution that have already been computed
 */
void q_contrib_gdWdW(int n_nod_elem, const float **Dw_x, float diff, const float **g_x, float **sum_contrib);


/**
 * @brief Evaluates the integrals on the element K with quadrature formula
 *
 * @param t          type of the element (1: quadrangle, 2:triangle)
 * @param n_quad_pts number of points in the quadrature formula
 * @param n_nod_elem number of nodes in the element K
 * @param a_K        nodes coordinates of element K
 * @param x_hat_quad coordinates of quadrature points in the reference element K_hat (!= K)
 * @param weights    array of weights associated to quadrature points
 * @param A_K_elem   matrix : part of discretization of A_K_ij that are integrals on K
 * @param l_K_elem   vector : part of discretization of l_K_i (2nd member) that are integrals on K
 */
void intElem(int t, int n_quad_pts, int n_nod_elem, const float **a_K, const float **x_hat_quad, const float *weights,
             float **A_K_elem, float *l_K_elem);

/**
 * @brief Evaluates the integrals on the edges of element K (Kp = K prime) with quadrature formula
 *
 * @param n_quad_pts number of points in the quadrature formula for an integral over a segment
 * @param a_Kp       nodes coordinates of edge Kp
 * @param x_hat_quad coordinates of quadrature points in the reference element Kp_hat (!= Kp)
 * @param weights    array of weights associated to quadrature points
 * @param A_K_edg    matrix : part of discretization of A_K_ij that are integrals on Kp
 * @param l_K_edg    vector : part of discretization of l_K_i (2nd member) that are integrals on K
 */
void intEdge(int n_quad_pts, const float **a_Kp, const float **x_hat_quad, const float *weights, float **A_K_edg,
             float *l_K_edg);


/**
 * @brief Evaluates A_K_ij and second member l_K_i form the discretization of the variational formulation for the element K
 *
 * @param ref_interior reference number for edges inside the domain
 * @param ref_Dh       array that contains number of domain boundaries carrying an homogeneous Dirichlet condition
 * @param ref_Dnh      array that contains number of domain boundaries carrying a non-homogeneous Dirichlet condition
 * @param ref_NF       array that contains number of domain boundaries carrying a Neumann or Fourier condition
 * @param n_Dh         size of ref_Dh
 * @param n_Dnh        size of ref_Dnh
 * @param n_NF         size of ref_NF
 * @param t            type of the element (1: quadrangle, 2:triangle)
 * @param n_nod_elem   number of nodes in the element K
 * @param a_K          nodes coordinates of element K
 * @param n_edg_elem   number of edges in the element K
 * @param ref_edg_K    array that contains reference number for each edge of K
 * @param A_K          matrix that contains values of A_ij for each couple of nodes (i,j) in [n_nod_elem]^2 of the element K
 * @param l_K          vector that contains values of l_i for each node i in [n_nod_elem] of the element K
 * @param nodes_D      array where the element i in [0,n_nod_elem] =0 if nod_i carries an homogeneous Dirichlet condition, -1 i nod_i carries a non-homogenous Dirichlet condition and 1 otherwise
 * @param uD_aK        array where the element i contains the value uD(a_K_i) if the nod_i carries a non-homogenous Dirichlet condition, 0 otherwise
 */
void
eval_K(int ref_interior, const int *ref_Dh, const int *ref_Dnh, const int *ref_NF, int n_Dh, int n_Dnh,
       int n_NF, int t,
       int n_nod_elem, const float **a_K, int n_edg_elem, const int *ref_edg_K, float **A_K, float *l_K,
       int *nodes_D,
       float *uD_aK
);


#endif //FEM_PROJECT_ELEM_EVAL_H
