#ifndef EULEREIGVALVEC_H
#define EULEREIGVALVEC_H

#include "defineGlobalVariable.h"
#include <omp.h>

REAL SetEigScalingAVE(const REAL *priL, const REAL *priR,
                      REAL &BigG, REAL (*MatR)[VEC_DIM], REAL *EigVec);

REAL MaxEigValue(REAL &BigG, const REAL *priU);

void AllEigValue(REAL &BigG, const REAL *priU, REAL *eig_value);

/*******************************************/
REAL SetEigVecScalingXX(const REAL *priVL, const REAL *priVR, REAL &BigG,
                        REAL (*MatR)[VEC_DIM]);

REAL SetEigVecScalingXX_(const REAL *priU, REAL &BigG, REAL (*MatR)[VEC_DIM]);

/*******************************************/

REAL ave_log(const REAL aL, const REAL aR);

#endif