#ifndef ALLOCATEFREEMEMORY_H
#define ALLOCATEFREEMEMORY_H

#include <cmath>
#include <string>

#include "defineGlobalVariable.h"

typedef REAL DP;


DP *Arr(int);

DP **Arr(int, int);

int **ArrInt(int, int);

DP ***Arr(int, int, int);

DP ****Arr(int, int, int, int);

void FreeArr(DP *);

void FreeArr(DP **);

void FreeArr(DP ***, int);

void FreeArr(DP ****, int, int);

int *IntArr(int);

void FreeIntArr(int *);


#endif