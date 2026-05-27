
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <stdarg.h>
#include <sys/stat.h>
#include <string>
#include "Euler1D.h"


using namespace std;

void Euler1D::Init() {

    REAL conU[VEC_DIM], priU[VEC_DIM];
    REAL xx;

    u = urk[0];


    BigG = 4 / 3.0;

    xxbl = 0.0, xxbr = 1.0;

    stop_t[0] = 0.4;

    dxx = (xxbr - xxbl) / npxx;

    cfl = 0.2;

    stop_ite = 1;

    for (int i = 0; i < npxx + gpn * 2; i++) {

        xx = xxbl + (i - gpn + 0.5) * dxx;

        if (xx <= 0.5) {
            priU[RHO_IND] = 1.0;
            priU[VV_XXIND] = 0.9;
            priU[PRE_IND] = 1.0;
        }
        else {
            priU[RHO_IND] = 1;
            priU[VV_XXIND] = 0;
            priU[PRE_IND] = 10;
        }

        computeConFromPri(BigG, priU, conU);


        for (int k = 0; k < VEC_DIM; k++) {
            u[i][k] = conU[k];
        }
    }


    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int k = 0; k < VEC_DIM; k++)
            urk[1][i][k] = u[i][k];

    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int k = 0; k < VEC_DIM; k++)
            urk[2][i][k] = u[i][k];


    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int k = 0; k < VEC_DIM; k++)
            urk[3][i][k] = u[i][k];
}

void Euler1D::Boundary() {

    // outflow bc
    for (int i = 0; i < gpn; i++) {
        for (int k = 0; k < VEC_DIM; k++) {
            u[i][k] = u[gpn][k];
        }
    }
    for (int i = npxx + gpn; i < npxx + 2 * gpn; i++) {
        for (int k = 0; k < VEC_DIM; k++) {
            u[i][k] = u[npxx + gpn - 1][k];
        }
    }
}

void Euler1D::OutputOrder() {

    ofstream os("RP_3_priU_RK3_EOS_1_ES5_400.dat", ios::ate);

    REAL conU[VEC_DIM], priU[VEC_DIM];
    REAL rho, xx_vel, pre;

    for (int i = gpn; i < npxx + gpn; i++) {

        for (int k = 0; k < VEC_DIM; k++) {
            conU[k] = u[i][k];
        }

        computePriFromCon(BigG, conU, priU);
        rho  = priU[RHO_IND];
        xx_vel = priU[VV_XXIND];
        pre = priU[PRE_IND];

        os << std::setprecision(16) << scientific;
        os << rho << " " << xx_vel << " " <<  pre << "\n";

    }
    os.close();
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx;
    int EC_order = 6, ES_order = 5;

    nx = 400;

    Euler1D compute(nx);
    compute.Init();
    compute.Evolve(EC_order, ES_order);

    return 0;
}



