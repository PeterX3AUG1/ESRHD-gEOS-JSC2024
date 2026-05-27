#ifndef EULER1D_H
#define EULER1D_H

#include <fstream>
#include <ostream>
#include <iostream>
#include <assert.h>
#include <memory>

#include "ConVPriVEntrV.h"
#include "allocateFreeMemory.h"
#include "setFlux.h"
#include "EulerEigValVec.h"
#include "defineGlobalVariable.h"

class Euler1D {
public:
    Euler1D(int nx);

    ~Euler1D();

    void Init();

    void Evolve(int &EC_order, int &ES_order);

    void Boundary();

    void OutputOrder();

    void OutputFile();

    void ComAlert();

private:
    const int npxx, gpn;
    REAL **urk[4], **u;
    REAL **flux, **upri, **fluxES, **fluxLF;

    int t_ite, stop_ite;
    REAL t, stop_t[5];
    REAL dxx;
    REAL cfl, BigG;

    REAL xxbl, xxbr;
    int tcase;
    char cell_name_pri[50];

    REAL **uEntr;

    REAL Global_MaxEig();

    void ComputeFlux(int &EC_order, int &ES_order);

    void computeSourceFlux(int &EC_order, int &index, REAL delta_s);

    void ComputFluxESAndLF(int &EC_order, int &ES_order);


};


#endif