#ifndef CONVPRIVENTRV_H
#define CONVPRIVENTRV_H

#include <iostream>
#include <math.h>
#include "defineGlobalVariable.h"

using namespace std;

REAL Hybrid_Root_p(REAL D, REAL mx, REAL my, REAL E, REAL BigG);

int computePriFromCon(REAL BigG, REAL *conU, REAL *priU);

int computeConFromPri(REAL BigG, REAL *priU, REAL *conU);

int computeEntrFromPri(REAL BigG, REAL *priU, REAL *Entr);


#endif