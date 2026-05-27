#ifndef DEFINEGLOBALVARIABLE_H
#define DEFINEGLOBALVARIABLE_H

#include <iostream>

using namespace std;

typedef double REAL;

#define VEC_DIM 4

#define DISSIPATION 0

#define RK_ORDER 3

#define XX_FLAG 1
#define YY_FLAG 0

const int RHO_IND = 0;

const int MM_XXIND = 1;
const int MM_YYIND = 2;

const int VV_XXIND = 1;
const int VV_YYIND = 2;

const int ENE_IND = 3;
const int PRE_IND = 3;


const int FF_RHO_IND = 0;
const int FF_MM_XXIND = 1;
const int FF_MM_YYIND = 2;
const int FF_ENE_IND = 3;
const int FF_PRE_IND = 3;

const REAL PI = 4.0 * atan(1.0);

#if RK_ORDER == 4
const REAL four_coe[4] = {0.5, 0.5, 1, 0.5};
#endif

#if RK_ORDER == 3
const int rk_ind[3] = {1, 2, 0};
const REAL rk_coe[3][2] = {{0,         1},
                           {0.75,      0.25},
                           {1.0 / 3.0, 2.0 / 3.0}};
#endif

#if RK_ORDER == 2
const int rk_ind[3]={1,0,0};
const REAL rk_coe[2][2]={{0,1},{0.5,0.5}};
#endif

#if RK_ORDER == 1
const int rk_ind[3] = {0, 2, 0};
const REAL rk_coe[1][2] = {{0, 1}};
#endif


#endif
