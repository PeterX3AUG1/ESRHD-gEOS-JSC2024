#ifndef SETFLUX_H
#define SETFLUX_H

#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <memory>
#include "LimiterRestrction.h"
#include "EulerEigValVec.h"
#include "ConVPriVEntrV.h"
#include "defineGlobalVariable.h"

#define USE_IDEN 0

using namespace std;


void setFluxXX(REAL &BigG, REAL *conU, REAL *flux);

void setFluxYY(REAL &BigG, REAL *conU, REAL *flux);


void setEntropyConserFluxXX(REAL &BigG, REAL *priVL, REAL *priVR, REAL *flux);

void setEntropyConserFluxYY(REAL &BigG, REAL *priVL, REAL *priVR, REAL *flux);

void dissipation(int flag, int &ES_order, REAL &BigG, const REAL *priL, const REAL *priR,
                 const REAL *Entr1, const REAL *Entr2, const REAL *Entr3, const REAL *Entr4,
                 const REAL *Entr5, const REAL *Entr6, const REAL *Entr7, const REAL *Entr8,
                 REAL *Dissip);

void setLFFluxXX(REAL &BigG, REAL *conVL, REAL *conVR, REAL *flux);

void setLFFluxYY(REAL &BigG, REAL *conVL, REAL *conVR, REAL *flux);

void setEntropyFlux(REAL &BigG, REAL *priVL, REAL *priVR, REAL *EntrVL, REAL *EntrVR,
                    REAL *flux, REAL &fluxQ);

int computeFES(REAL *Entr,
               REAL *fluxHL, REAL *fluxHR,
               REAL *fluxLL, REAL *fluxLR,
               REAL &QHL, REAL &QHR,
               REAL &QLL, REAL &QLR,
               REAL *FES);


#endif