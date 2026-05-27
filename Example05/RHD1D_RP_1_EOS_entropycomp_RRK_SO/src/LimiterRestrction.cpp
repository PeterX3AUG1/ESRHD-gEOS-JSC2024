
#include <iostream>
#include <cmath>
#include "LimiterRestrction.h"

using namespace std;


/** next: WENO -- 4 point Gauss-Lobatto points
 ****************************************************
 *****************************************************
 */
/**
 *
 *
 * @param v5 cell average u_{i-2}
 * @param v4
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v2
 * @param v1
 *
 * @return  point value u_{i-1/2}
 */
REAL WENOGL41(REAL v5, REAL v4, REAL v3, REAL v2, REAL v1) {
    REAL s1, s2, s3;
    REAL sum;
    REAL tmp1, tmp2;
    REAL eps = 1e-6;
    REAL w1, w2, w3;

    tmp1 = v1 - 2 * v2 + v3;
    tmp2 = v1 - 4 * v2 + 3 * v3;
    s1 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s1 = 0.1 / (s1 * s1);

    tmp1 = v2 - 2 * v3 + v4;
    tmp2 = v2 - v4;
    s2 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s2 = 0.6 / (s2 * s2);

    tmp1 = v3 - 2 * v4 + v5;
    tmp2 = 3 * v3 - 4 * v4 + v5;
    s3 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s3 = 0.3 / (s3 * s3);

    //  s1=0.1;
    //s2=0.6;
    //s3=0.3;

    sum = s1 + s2 + s3;
    w1 = s1 / sum;
    w2 = s2 / sum;
    w3 = s3 / sum;

    return (w1 * (v1 / 3 - 7 * v2 / 6 + 11 * v3 / 6) + w2 * (-v2 / 6 + 5 * v3 / 6 + v4 / 3) +
            w3 * (v3 / 3 + 5 * v4 / 6 - v5 / 6));
}

/**
 *
 *
 * @param v5 cell average u_{i-2}
 * @param v4
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v2
 * @param v1
 *
 * @return  point value u_{i-sqrt(5)/10}
 */
REAL WENOGL42(REAL v5, REAL v4, REAL v3, REAL v2, REAL v1) {
    REAL s1, s2, s3;
    REAL sum;
    REAL tmp1, tmp2;
    REAL w1, w2, w3;
    REAL value;
    REAL sq5 = sqrt(5.0);
    REAL eps = 1e-14;

    tmp1 = v1 - 2 * v2 + v3;
    tmp2 = v1 - 4 * v2 + 3 * v3;
    s1 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;

    s1 = (91 + 9 * sq5) / (440 * s1 * s1);

    tmp1 = v2 - 2 * v3 + v4;
    tmp2 = v2 - v4;
    s2 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;

    s2 = 129. / (220 * s2 * s2);

    tmp1 = v3 - 2 * v4 + v5;
    tmp2 = 3 * v3 - 4 * v4 + v5;
    s3 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s3 = (91 - 9 * sq5) / (440 * s3 * s3);

    sum = s1 + s2 + s3;
    w1 = s1 / sum;
    w2 = s2 / sum;
    w3 = s3 / sum;

    value = w1 * ((-1. / 60 + 0.05 * sq5) * v1 + (1. / 30 - sq5 * 0.2) * v2 + (59. / 60 + 0.15 * sq5) * v3) +
            w2 * ((-1. / 60 - 0.05 * sq5) * v2 + 31. / 30 * v3 + (-1.0 / 60 + 0.05 * sq5) * v4) +
            w3 * ((59. / 60 - 0.15 * sq5) * v3 + (1. / 30 + 0.2 * sq5) * v4 + (-1. / 60 - 0.05 * sq5) * v5);
    return value;
}

/**
 *
 *
 * @param v1 cell average u_{i-2}
 * @param v2
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v4
 * @param v5
 *
 * @return  point value u_{i+sqrt(5)/10}
 */
REAL WENOGL43(REAL v1, REAL v2, REAL v3, REAL v4, REAL v5) {
    REAL s1, s2, s3;
    REAL sum;
    REAL tmp1, tmp2;
    REAL w1, w2, w3;
    REAL value;
    REAL sq5 = sqrt(5.0);
    REAL eps = 1e-14;

    tmp1 = v1 - 2 * v2 + v3;
    tmp2 = v1 - 4 * v2 + 3 * v3;
    s1 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;

    s1 = (91 + 9 * sq5) / (440 * s1 * s1);

    tmp1 = v2 - 2 * v3 + v4;
    tmp2 = v2 - v4;
    s2 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;

    s2 = 129. / (220 * s2 * s2);

    tmp1 = v3 - 2 * v4 + v5;
    tmp2 = 3 * v3 - 4 * v4 + v5;
    s3 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s3 = (91 - 9 * sq5) / (440 * s3 * s3);

    sum = s1 + s2 + s3;
    w1 = s1 / sum;
    w2 = s2 / sum;
    w3 = s3 / sum;

    value = w1 * ((-1. / 60 + 0.05 * sq5) * v1 + (1. / 30 - sq5 * 0.2) * v2 + (59. / 60 + 0.15 * sq5) * v3) +
            w2 * ((-1. / 60 - 0.05 * sq5) * v2 + 31. / 30 * v3 + (-1.0 / 60 + 0.05 * sq5) * v4) +
            w3 * ((59. / 60 - 0.15 * sq5) * v3 + (1. / 30 + 0.2 * sq5) * v4 + (-1. / 60 - 0.05 * sq5) * v5);
    return value;
}


/**
 *
 *
 * @param v1 cell average u_{i-2}
 * @param v2
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v4
 * @param v5
 *
 * @return  point value u_{i+1/2}
 */
REAL WENOGL44(REAL v1, REAL v2, REAL v3, REAL v4, REAL v5) {
    REAL s1, s2, s3;
    REAL sum;
    REAL tmp1, tmp2;
    REAL eps = 1e-6;
    REAL w1, w2, w3;

    tmp1 = v1 - 2 * v2 + v3;
    tmp2 = v1 - 4 * v2 + 3 * v3;
    s1 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s1 = 0.1 / (s1 * s1);

    tmp1 = v2 - 2 * v3 + v4;
    tmp2 = v2 - v4;
    s2 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s2 = 0.6 / (s2 * s2);

    tmp1 = v3 - 2 * v4 + v5;
    tmp2 = 3 * v3 - 4 * v4 + v5;
    s3 = 13.0 / 12 * tmp1 * tmp1 + 0.25 * tmp2 * tmp2 + eps;
    s3 = 0.3 / (s3 * s3);

    // s1=0.1;
    //s2=0.6;
    //s3=0.3;

    sum = s1 + s2 + s3;
    w1 = s1 / sum;
    w2 = s2 / sum;
    w3 = s3 / sum;

    return (w1 * (v1 / 3 - 7 * v2 / 6 + 11 * v3 / 6) + w2 * (-v2 / 6 + 5 * v3 / 6 + v4 / 3) +
            w3 * (v3 / 3 + 5 * v4 / 6 - v5 / 6));
}


/** next: WENO PW
 ****************************************************
 *****************************************************
 */
/**
 *
 *
 * @param v1 mid point values average u_{i-2}
 * @param v2
 * @param v3  mid point values u_i , cell i: [i-1/2,i+1/2]
 * @param v4
 * @param v5
 *
 * @return  point value u_{i+1/2}
 */
REAL WENOPW44(REAL v1, REAL v2, REAL v3, REAL v4, REAL v5) {
    REAL s1, s2, s3;
    REAL sum;
    REAL tmp1, tmp2;
    REAL eps = 1e-12;
    REAL w1, w2, w3;

    tmp1 = (4.0 * v1 * v1 - 19.0 * v1 * v2 + 25.0 * v2 * v2 + 11.0 * v1 * v3 - 31.0 * v2 * v3 + 10.0 * v3 * v3) / 3.0;
    s1 = fabs(tmp1) + eps;
    s1 = (1.0 / 16.0) / (s1 * s1);

    tmp1 = (4.0 * v2 * v2 - 13.0 * v2 * v3 + 13.0 * v3 * v3 + 5.0 * v2 * v4 - 13.0 * v3 * v4 + 4.0 * v4 * v4) / 3.0;
    s2 = fabs(tmp1) + eps;
    s2 = (5.0 / 8.0) / (s2 * s2);

    tmp1 = (10.0 * v3 * v3 - 31.0 * v3 * v4 + 25.0 * v4 * v4 + 11.0 * v3 * v5 - 19.0 * v4 * v5 + 4.0 * v5 * v5) / 3.0;
    s3 = fabs(tmp1) + eps;
    s3 = (5.0 / 16.0) / (s3 * s3);

    // s1=1.0/16.0;
    //s2=5.0/8.0;
    //s3=5.0/16.0;

    sum = s1 + s2 + s3;
    w1 = s1 / sum;
    w2 = s2 / sum;
    w3 = s3 / sum;


    return (w1 * (0.375 * v1 - 1.25 * v2 + 1.875 * v3) + w2 * (-0.125 * v2 + 0.75 * v3 + 0.375 * v4) +
            w3 * (0.375 * v3 + 0.75 * v4 - 0.125 * v5));
}

/**
 *
 *
 * @param v5 mid point values average u_{i-2}
 * @param v4
 * @param v3  mid point values u_i , cell i: [i-1/2,i+1/2]
 * @param v2
 * @param v1
 *
 * @return  point value u_{i-1/2}
 */
REAL WENOPW41(REAL v5, REAL v4, REAL v3, REAL v2, REAL v1) {
    REAL s1, s2, s3;
    REAL sum;
    REAL tmp1, tmp2;
    REAL eps = 1.0E-12;
    REAL w1, w2, w3;

    tmp1 = (4.0 * v1 * v1 - 19.0 * v1 * v2 + 25.0 * v2 * v2 + 11.0 * v1 * v3 - 31.0 * v2 * v3 + 10.0 * v3 * v3) / 3.0;
    s1 = fabs(tmp1) + eps;
    s1 = (1.0 / 16.0) / (s1 * s1);

    tmp1 = (4.0 * v2 * v2 - 13.0 * v2 * v3 + 13.0 * v3 * v3 + 5.0 * v2 * v4 - 13.0 * v3 * v4 + 4.0 * v4 * v4) / 3.0;
    s2 = fabs(tmp1) + eps;
    s2 = (5.0 / 8.0) / (s2 * s2);

    tmp1 = (10.0 * v3 * v3 - 31.0 * v3 * v4 + 25.0 * v4 * v4 + 11.0 * v3 * v5 - 19.0 * v4 * v5 + 4.0 * v5 * v5) / 3.0;
    s3 = fabs(tmp1) + eps;
    s3 = (5.0 / 16.0) / (s3 * s3);

    //s1=1.0/16.0;
    //s2=5.0/8.0;
    //s3=5.0/16.0;

    sum = s1 + s2 + s3;
    w1 = s1 / sum;
    w2 = s2 / sum;
    w3 = s3 / sum;


    return (w1 * (0.375 * v1 - 1.25 * v2 + 1.875 * v3) + w2 * (-0.125 * v2 + 0.75 * v3 + 0.375 * v4) +
            w3 * (0.375 * v3 + 0.75 * v4 - 0.125 * v5));
}


/** next: ENO4
 ****************************************************
 *****************************************************
 */
REAL ENO4UR(REAL u1, REAL u2, REAL u3, REAL u4,
            REAL u5, REAL u6, REAL u7) {
    //ENO coefficients of u_j for a uniform mesh
    // coefficients for point-wise ENO
    //const REAL Crj[4][5]={
    //{ 2.1875,  0.3125, -0.0625,  0.0625, -0.3125},
    //{ -2.1875,  0.9375,  0.5625, -0.3125,  1.3125},
    //{  1.3125, -0.3125,  0.5625,  0.9375, -2.1875},
//	{ -0.3125,  0.0625, -0.0625,  0.3125,  2.1875} };

    const REAL Crj[4][5] = {
            {2.0833333333333333333333333333333,  0.25,                                -0.083333333333333333333333333333333, 0.083333333333333333333333333333333, -0.25},
            {-1.9166666666666666666666666666667, 1.0833333333333333333333333333333,   0.58333333333333333333333333333333,   -0.41666666666666666666666666666667, 1.0833333333333333333333333333333},
            {1.0833333333333333333333333333333,  -0.41666666666666666666666666666667, 0.58333333333333333333333333333333,   1.0833333333333333333333333333333,   -1.9166666666666666666666666666667},
            {-0.25,                              0.083333333333333333333333333333333, -0.083333333333333333333333333333333, 0.25,                                2.0833333333333333333333333333333}};


    const int k = 4;
    const int n = 7;

    REAL u[7];
    REAL dd[4][7];

    u[0] = u1;
    u[1] = u2;
    u[2] = u3;
    u[3] = u4;
    u[4] = u5;
    u[5] = u6;
    u[6] = u7;


    int i, j;
    for (i = 0; i < k; i++)
        for (j = 0; j < n; j++)
            dd[i][j] = 0.0;

    // Compute undivided difference of the primitive of u
    for (j = 0; j < n; j++)
        dd[0][j] = u[j];

    for (j = 2; j <= k; j++)
        for (i = 1; i <= (n + 1 - j); i++)
            dd[j - 1][i - 1] = dd[j - 2][i] - dd[j - 2][i - 1];


    // Determine r, the left-displacement of index

    // Indices to each stencil
    int curR = 0;
    int l;

    for (l = 2; l <= k; l++) {
        if (fabs(dd[l - 1][k - curR - 2]) <= fabs(dd[l - 1][k - curR - 1]))
            curR++;
    }


    REAL ur = 0.0;

    for (i = 0; i < k; i++)
        ur += u[k - curR + i - 1] * Crj[i][curR + 1];

    return ur;

}


REAL ENO4UL(REAL u1, REAL u2, REAL u3, REAL u4,
            REAL u5, REAL u6, REAL u7) {
    REAL ul;

    ul = ENO4UR(u7, u6, u5, u4, u3, u2, u1);

    return ul;

#if 0
    //ENO coefficients of u_j for a uniform mesh
    const REAL Crj[4][5]={
    { 2.1875,  0.3125, -0.0625,  0.0625, -0.3125},
    { -2.1875,  0.9375,  0.5625, -0.3125,  1.3125},
    {  1.3125, -0.3125,  0.5625,  0.9375, -2.1875},
    { -0.3125,  0.0625, -0.0625,  0.3125,  2.1875} };

    const int k=4;
    const int n=7;

    REAL u[7];
    REAL dd[4][7];

    u[0]=u1;
    u[1]=u2;
    u[2]=u3;
    u[3]=u4;
    u[4]=u5;
    u[5]=u6;
    u[6]=u7;


    int i,j;
    for ( i=0; i<k; i++ )
        for ( j=0; j<n; j++ )
            dd[i][j] = 0.0;

    // Compute undivided difference of the primitive of u
    for ( j=0; j<n; j++ )
        dd[0][j] = u[j];

    for ( j=2; j<=k; j++ )
        for ( i=1; i<=(n+1-j); i++ )
            dd[j-1][i-1] = dd[j-2][i] - dd[j-2][i-1];


    // Determine r, the left-displacement of index

    // Indices to each stencil
    int curR = 0;
    int l;

    for ( l=2; l<=k; l++ )
        if (  fabs( dd[l-1][k-curR-2] ) <= fabs( dd[l-1][k-curR-1] )  )
            curR++;

    REAL ul=0.0;

    for ( i=0; i<k; i++ )
        ul += u[ k-curR + i - 1 ] *  Crj[i][curR];


    return ul;
#endif

}

REAL Minmod(REAL a, REAL b) {

    if (a > 0) {
        if (b > 0)
            return a > b ? b : a;
        else
            return 0.0;
    } else {
        if (b > 0)
            return 0.0;
        else
            return a > b ? a : b;
    }
}
