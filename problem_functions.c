#include "problem_functions.h"
#include <math.h>
#include <stdio.h>

float a00(const float* x) { return 0;}

float a11(const float* x) { return 1;}

float a12(const float* x) { return 0;}

float a22(const float* x) { return 1;}

float bN(const float* x) { return 0.;}

float fN(const float* x) { return 0.;} // fN = 0 for pb 3 domain 1

// for DOMAIN 2
/*
float fN(const float* x) {
    if (x[0]==0.5){
        return -M_PI*cosf(M_PI*x[1]);
    }else if (x[1] == 0.5){
        return -M_PI*cosf(M_PI*x[0]);
    }
    else{
        printf("WARNING, wrong case of use for this fN, x = %f, y = %f\n",x[0],x[1]);
        return 0.;
    }

}
*/
//float uD(const float* x) { return 100*x[0] + x[1];}


// for domain 2

float uD(const float* x) {
    float val = 0.;
    switch(nucas) {
        case 1 :
            val = 16. * x[0] * x[1] * (1 - x[0]) * (1 - x[1]);
            break;
        case 2 :
            val = sinf(M_PI*x[0])*sinf(M_PI*x[1]);
            break;
        case 3 :
            val=cosf(M_PI*x[0])*cosf(M_PI*x[1]);
            break;
        default:
            printf("*** uD : unknown example .\n");
            break;
    }
    return val;
}


// old for first test (before results)
//float f_OMEGA(const float* x) { return 0;}

float f_OMEGA(const float* x)
{
    const float PI = M_PI;
    float val = 0.;

    switch (nucas) {
        case 1 :
            val=-32.*(x[0]*(x[0]-1)+x[1]*(x[1]-1));;
            break;
        case 2 :
            val= 2.*sinf(PI*x[0])*sinf(PI*x[1])*PI*PI;
            break;
        case 3 :
            val=(1+2*PI*PI)*cosf(PI*x[0])*cosf(PI*x[1]);
            break;
        default :
            printf("*** SOLEX : exemple non prevu.\n");
            break;
    }
    return(val);
}

/// @brief exact solution function
float exact_sol(const float *x)
{
    const float PI = M_PI;
    float val = 0.;

    switch (nucas) {
        case 1 :
            val=16.*x[0]*x[1]*(1-x[0])*(1-x[1]);
            break;
        case 2 :
            val=sinf(PI*x[0])*sinf(PI*x[1]);
            break;
        case 3 :
            val=cosf(PI*x[0])*cosf(PI*x[1]);
            break;
        default :
            printf("*** SOLEX : exemple non prevu.\n");
            break;
    }
    return(val);
}
