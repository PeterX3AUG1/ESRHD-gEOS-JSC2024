#ifndef EULEREIGVALVEC_H
#define EULEREIGVALVEC_H

#include "defineGlobalVariable.h"


REAL SetEigScalingAVEXX(const REAL *priL, const REAL *priR,
                      REAL &BigG, REAL (*MatR)[VEC_DIM], REAL *EigVec);

REAL SetEigScalingAVEYY(const REAL *priL, const REAL *priR,
                        REAL &BigG, REAL (*MatR)[VEC_DIM], REAL *EigVec);

REAL MaxEigValue(REAL &BigG, int flag, const REAL *priU);

void AllEigValue(REAL &BigG, int flag, const REAL *priU, REAL *eig_value);

/*******************************************/
REAL SetEigVecScalingXX(const REAL *priVL, const REAL *priVR, REAL &BigG,
                        REAL (*MatR)[VEC_DIM]);

REAL SetEigVecScalingYY(const REAL *priVL, const REAL *priVR, REAL &BigG,
                        REAL (*MatR)[VEC_DIM]);

/*******************************************/

REAL ave_log(const REAL aL, const REAL aR);

#endif