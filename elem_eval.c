#include "elem_eval.h"
#include "tab_mngmt.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// !! In these functions, tab allocation is supposed to be done outside the function with correct dimension

// Problem functions :
float a00(const float* x) { return 1;}

float a11(const float* x) { return 1;}

float a12(const float* x) { return 0;}

float a22(const float* x) { return 1;}

float bN(const float* x) { return 0;}

float fN(const float* x) { return 1;}

float f_OMEGA(const float* x) { return 0;}

float uD(const float* x) { return 100*x[0] + x[2];}


// evaluate every base function w_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalFbase(int t, const float *x_hat, float *w_hat_x_hat)
{
    switch (t) {
        case 1: // Quadrangle
            w_hat_x_hat[0] = x_hat[0]*(1-x_hat[1]);
            w_hat_x_hat[1] = x_hat[0]*x_hat[1];
            w_hat_x_hat[2] = (1-x_hat[0])*x_hat[1];
            w_hat_x_hat[3] = (1-x_hat[0])*(1-x_hat[1]);
            break;
        case 2: // Triangle
            w_hat_x_hat[0] = x_hat[0];
            w_hat_x_hat[1] = x_hat[1];
            w_hat_x_hat[2] = 1-x_hat[0]-x_hat[1];
            break;
        case 3: // Segment
            w_hat_x_hat[0] = x_hat[0];
            w_hat_x_hat[1] = (1-x_hat[0]);
            break;
        default:
            printf("Unknown type t = %d for 'evalFbase' function",t);
            perror("Error : Unvalid parameter");
    }
}

// evaluate every gradient of base functions:  Dw_i (i in [1,q]) for a given point x_hat in the coordinates of the reference element
void evalDFbase(int t, const float *x_hat, float **Dw_hat_x_hat)
{
    switch (t) {
        case 1: // Quadrangle
            Dw_hat_x_hat[0][0] = 1-x_hat[1];    Dw_hat_x_hat[0][1] = x_hat[0]*(1-x_hat[1]);
            Dw_hat_x_hat[1][0] = x_hat[1];      Dw_hat_x_hat[1][1] = x_hat[0];
            Dw_hat_x_hat[2][0] = -x_hat[1];     Dw_hat_x_hat[2][1] = 1-x_hat[0];
            Dw_hat_x_hat[3][0] = x_hat[1]-1;    Dw_hat_x_hat[3][1] = x_hat[0]-1;
            break;
        case 2: // Triangle
            Dw_hat_x_hat[0][0] = 1;             Dw_hat_x_hat[0][1] = 0;
            Dw_hat_x_hat[1][0] = 0;             Dw_hat_x_hat[0][1] = 1;
            Dw_hat_x_hat[2][0] = -1;            Dw_hat_x_hat[0][1] = -1;
            break;
        case 3: // Segment (! Matrix of size 2*1)
            Dw_hat_x_hat[0][0] = -1;
            Dw_hat_x_hat[1][0] = -1;
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
    if (det < 10e-7){
        perror("Matrix is almost singular or not invertible (inv_2x2) call");
        return(1);
    }else{
        float inv_det = 1/det;
        M_inv[0][0] = inv_det*M[1][1];   M_inv[0][1] = -inv_det*M[0][1];
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
void wp_quad(int t, float **x_quad_hat, float *weight)
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
            x_quad_hat[0][0] = 1;     x_quad_hat[0][1] = 0;   // x1
            x_quad_hat[1][0] = 1;     x_quad_hat[1][1] = 1;   // x2
            x_quad_hat[2][0] = 0;     x_quad_hat[2][1] = 1;   // x3
            x_quad_hat[3][0] = 0;     x_quad_hat[3][1] = 0;   // x4
            x_quad_hat[4][0] = 1;     x_quad_hat[4][1] = 0.5; // x5
            x_quad_hat[5][0] = 0.5;   x_quad_hat[5][1] = 1;   // x6
            x_quad_hat[6][0] = 0;     x_quad_hat[6][1] = 0.5; // x7
            x_quad_hat[7][0] = 0.5;   x_quad_hat[7][1] = 0;   // x8
            x_quad_hat[8][0] = 0.5;   x_quad_hat[8][1] = 0.5; // x9
            break;

        case 2: // Triangle
            // weights
            weight_1 = 1.0/6;
            for (i=0; i<3; i++) weight[i] = weight_1;
            // quadrature points (x_i_hat)
            x_quad_hat[0][0] = 0.5;   x_quad_hat[0][1] = 0.5;  // x1
            x_quad_hat[1][0] = 0;     x_quad_hat[1][1] = 0.5;  // x2
            x_quad_hat[2][0] = 0.5;   x_quad_hat[2][1] = 0;    // x3
            break;

        case 3: // Segment
            // weights
            weight_1 = 1.0/6;
            weight[0] = weight_1;
            weight[0] = weight_1;
            weight[0] = 4*weight_1;
            // quadrature points (x_i_hat)
            x_quad_hat[0][0] = 1;    // x1
            x_quad_hat[1][0] = 0;    // x2
            x_quad_hat[2][0] = 0.5;  // x3
            break;

        default:
            printf("Unknown type t = %d for 'quad_order' function",t);
            perror("Error : Unvalid parameter");
    }
}

// evaluate transformation Fk of a point of the reference element x_hat
void transFK(int q, float **a_K, const float *w_hat_x_hat, float *Fk_x_hat)
{
    Fk_x_hat[0] = 0;
    Fk_x_hat[1] = 0;
    for (int i=0; i<q; i++){
        Fk_x_hat[0] += a_K[i][0] * w_hat_x_hat[i];
        Fk_x_hat[0] += a_K[i][1] * w_hat_x_hat[i];
    }
}

// evaluate the jacobian matrix of the transformation Fk of a point of the reference element x_hat
void jacobFK(int p, int d, const float **a_K, const float **Dw_hat_x_hat, float **jacob_Fk_x_hat) {
    int i;
    if (d==2) { // quadrangle, triangle
        // initialisation
        jacob_Fk_x_hat[0][0] = 0; jacob_Fk_x_hat[0][1] = 0;
        jacob_Fk_x_hat[1][0] = 0; jacob_Fk_x_hat[1][1] = 0;
        for(i=0; i<p; i++){
            jacob_Fk_x_hat[0][0] += a_K[i][0] * Dw_hat_x_hat[i][0]; jacob_Fk_x_hat[0][1] += a_K[i][0] * Dw_hat_x_hat[i][1];
            jacob_Fk_x_hat[1][0] += a_K[i][1] * Dw_hat_x_hat[i][0]; jacob_Fk_x_hat[1][1] += a_K[i][1] * Dw_hat_x_hat[i][1];
        }
    } else if (d==1){ // segment
        jacob_Fk_x_hat[0][0] = 0;
        jacob_Fk_x_hat[1][0] = 0;
        for(i=0; i<p; i++){
            jacob_Fk_x_hat[0][0] += a_K[i][0] * Dw_hat_x_hat[i][0];
            jacob_Fk_x_hat[1][0] += a_K[i][1] * Dw_hat_x_hat[i][0];
        }
    } else {
        printf("Unknown dimension d = %d for 'jacobFK' function \n",d);
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

// Compute quadrature contribution of the quadrature point x for each base function w_i for the integral of the shape w_i*g on an edge
void q_contrib_gW(int n_nod_elem, const float *w_x, float diff, float g_x, float *sum_contrib)
{
    for(int i=0; i<n_nod_elem; i++) sum_contrib[i] += diff * g_x * w_x[i];
}

// Compute quadrature contribution of the quadrature point x for each base function w_i for the integral of the shape w_i*g
void q_contrib_gWW(int n_nod_elem, const float *w_x, float diff, float g_x, float **sum_contrib)
{
    int i,j;
    float base, contrib;
    base = diff * g_x;
    for(i=0; i<n_nod_elem; i++){
        contrib = base * w_x[i];
        for(j=0; j<n_nod_elem; j++){
            sum_contrib[i][j] += contrib * w_x[j];
        }
    }
}

// Compute quadrature contribution of the quadrature point x for each base function w_i for the integral of the shape
// sum (alpha,beta = 1 -> 2) g_{alpha,beta}(x)*d_alpha(w_i)*d_beta(w_j)
void q_contrib_gdWdW(int n_nod_elem, const float **Dw_x, float diff, const float **g_x, float **sum_contrib)
{
    int i,j,alpha,beta;
    float base, contrib;

    for(alpha=0; alpha<2; alpha++){
        for(beta=0; beta<2; beta++){
            base = diff * g_x[alpha][beta];
            for(i=0; i<n_nod_elem; i++){
                contrib = base * Dw_x[i][alpha];
                for(j=0; j<n_nod_elem; j++){
                    sum_contrib[i][j] = contrib * Dw_x[j][beta];
                }
            }
        }
    }
}

// Evaluate the integrals on the element K with quadrature formula
void intElem(
        int t,
        int n_quad_pts,
        int n_nod_elem,
        const float **a_K,
        const float** x_hat_quad,
        const float *weights,
        float (*a00)(const float*),
        float (*a11)(const float*),
        float (*a12)(const float*),
        float (*a22)(const float*),
        float (*f_OMEGA)(const float*),
        float** A_K_elem,
        float* l_K_elem
        ) // I didn't put const before name of the function to not overload function signature
{
    int k,i,j;
    float detJF,diff;

    // Memory allocation for object reset in each quadrature point
    float *x_quad = malloc(2*sizeof(float));
    float *w_hat_x_hat = malloc(n_nod_elem*sizeof(float));
    float **Dw_hat_x_hat = matF_alloc(n_nod_elem,2);
    float **Dw_x_quad = matF_alloc(n_nod_elem,2);
    float **JFk_x_hat = matF_alloc(2,2);
    float **inv_JFk_x_hat = matF_alloc(2,2);
    float **a_alpha_beta = matF_alloc(2,2);

    // Get node coordinates of the element (a_K) -> (1) maybe get them before ??
    // -> yes because when a_K is freed, the coordinates are deleted because it's only a pointers copy

    for(k=0; k<n_quad_pts; k++) {
        transFK(n_nod_elem,(const float**)a_K, w_hat_x_hat, x_quad); // Get the quadrature points coordinates x_quad_k with F_K

        // evaluation of bases function et gradient of bases functions in the quadrature point x_hat_quad
        evalFbase(t,x_hat_quad[k],w_hat_x_hat);
        evalDFbase(t,x_hat_quad[k],Dw_hat_x_hat);

        // Evaluation Dw(x_quad)/Dx_alpha, alpha in {1,2};
        jacobFK(n_nod_elem, 2, (const float **) a_K, (const float **) Dw_hat_x_hat, JFk_x_hat);
        detJF = inv_2x2(JFk_x_hat,inv_JFk_x_hat);
        diff = fabsf(detJF)*weights[k]; // differential elem * weight for the quadrature formula
        for(j=0; j<2; j++){
            for(i=0; i<n_nod_elem; i++){
                Dw_x_quad[i][j] = Dw_hat_x_hat[i][0]*inv_JFk_x_hat[0][j] + Dw_hat_x_hat[i][1]*inv_JFk_x_hat[1][j];
            }
        }

        // pre-compute a_{alpha,beta}(x_quad)
        a_alpha_beta[0][0] = a11(x_quad); a_alpha_beta[0][1] = a12(x_quad);
        a_alpha_beta[1][0] = a12(x_quad); a_alpha_beta[1][1] = a22(x_quad);

        // Contribution of the quadrature point for the quadrature formula
        q_contrib_gW(n_nod_elem,w_hat_x_hat,diff,f_OMEGA(x_quad),l_K_elem);
        q_contrib_gWW(n_nod_elem,(const float*)w_hat_x_hat,diff,a00(x_quad),A_K_elem);
        q_contrib_gdWdW(n_nod_elem,(const float**)Dw_hat_x_hat,diff,(const float**)a_alpha_beta,A_K_elem);

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


// Evaluate the integrals on the edges of element K (K prime) with quadrature formula
void intSegment(
        int n_quad_pts,
        int n_nod_seg,
        const float **nodes_coords,
        const float **x_hat_quad, // (2) -> dim = 1 ?
        const float *weights,
        float (*bN)(const float*),
        float (*fN)(const float*),
        float** A_K_edg,
        float* l_K_edg
        )
{
}

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
        )
{

    int i,j,n_quad_pts_elem,n_quad_pts_edg;




    // ----------- Integral on the element K -----------

    n_quad_pts_elem = quad_order(t);
    // Memory
    float **x_quad_hat_elem = matF_alloc(n_quad_pts_elem,2);
    float *weights_elem = malloc(n_quad_pts_elem*sizeof(float));
    float **A_K_elem = matF_alloc0(n_nod_elem,n_nod_elem);  // Every value set to 0
    float *l_K_elem = calloc(n_nod_elem,sizeof(float));     // Every value set to 0

    wp_quad(t,x_quad_hat_elem,weights_elem); // evaluate weights and quadrature pts for the quadrature on the element
    intElem(t,n_quad_pts_elem, n_nod_elem,(const float **)a_K,(const float**)x_quad_hat_elem, weights_elem, a00, a11, a12, a22, f_OMEGA, A_K_elem, l_K_elem);

    // Update results :
    for(i=0; i<n_nod_elem; i++) {
        l_K[i] += l_K_elem[i];
        for(j=0; j<n_nod_elem; j++) A_K[i][j] += A_K_elem[i][j];
    }

    // free Memory (weights_elem will be reallocated)
    free_mat(x_quad_hat_elem);
    free_mat(A_K_elem);
    free(l_K_elem);


    // ----------- Integral on a potential edge element K

    // Memory update (clean and new) // (5) Should I define n_nod_edg = 2 and replace everywhere ?

    float **x_quad_hat_edg = matF_alloc(n_quad_pts_edg,2);
    float *weights_edg = realloc(weights_edg,n_quad_pts_edg*sizeof(float));
    int *edg_nodes = malloc(2*sizeof(int));     // 2 because we take order 1 for approximation of edges
    float **edg_nodes_coords = matF_alloc(2,2); // same (and we're in dim 2)
    float **A_K_edg = matF_alloc0(2,2);         // Every value set to 0 // (3) --> what is the correct size : n_nod_elem or n_nod_edg ?? Seems to me to be n_nod_edg because in intSegment we evaluate base fonction with type 3 (segment) ==> nb base fonction = 2 = n_nod_edg
    float *l_K_edg = calloc(2,sizeof(float));   // Every value set to 0

    for(i=0; i<n_edg_elem; i++){
        // Check if the edge is on the edge of the domain (i.e. is particular)
        for(j=0; j<n_Dh; j++) if(ref_edg_K[i]==ref_Dh[j]) {nodes_D[i]=0;}                          // Dirichlet homogeneous
        for(j=0; j<n_Dh; j++) if(ref_edg_K[i]==ref_Dnh[j]) {nodes_D[i]=-1; uD_aK[i] = uD(a_K[i]);} // Dirichlet non-homogeneous
        for(j=0; j<n_NF; j++) if(ref_edg_K[i]==ref_NF[j]) {                                        // Neumann or Fourier
            vertices_Edge(i+1,t,edg_nodes);
            selectPts(2,edg_nodes,(float**)a_K,edg_nodes_coords);       // Careful we lost the const on a_K here (2 because we take order 1 for approximation of edges)
            wp_quad(3,x_quad_hat_edg,weights_edg);
            intSegment(n_quad_pts_edg,2,(const float**)edg_nodes_coords,(const float**)x_quad_hat_edg, weights_edg, bN, fN, A_K_edg, l_K_edg);

            // Update results (we had the contribution for each node at the end of the edge (2 nodes)
            for(i=0; i<2; i++) {
                l_K[edg_nodes[i]] += l_K_edg[i];
                for(j=0; j<2; j++) A_K[i][j] += A_K_edg[edg_nodes[i]][edg_nodes[j]];
            }
        }
    }

    // free memory
    free_mat(x_quad_hat_edg);
    free(l_K_edg);
    free(edg_nodes);
    free_mat(edg_nodes_coords);
    free_mat(A_K_edg);
    free(weights_edg);
}



// (4) : Hypothese : We suppose that, dynamical tab of sum of contribution are initialized to 0
