
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


    BigG = 5 / 3.0;

    xxbl = -0.4, xxbr = 2.0;

    stop_t[0] = 0.8;

    dxx = (xxbr - xxbl) / npxx;

    cfl = 0.2;

    stop_ite = 1;


    REAL L = 0.3, K = 100;

#pragma omp parallel for private (xx,priU,conU)
    for (int i = 0; i < npxx + gpn * 2; i++) {

        xx = xxbl + (i - gpn + 0.5) * dxx;

        if ( xx >= L || xx <= -L) {
            priU[RHO_IND] = 1.0;
            priU[PRE_IND] = K;
            priU[VV_XXIND] = 0.0;
        }
        else {
            priU[RHO_IND] = 1.0 + exp( -1./ (1.0 - xx * xx / L / L) );
            priU[PRE_IND] = K * pow( priU[RHO_IND], BigG);
            REAL tpA = sqrt(K*BigG*(BigG-1)/(K*BigG+BigG-1)), tpKGAMMA = BigG*K*pow(priU[RHO_IND],BigG-1);
            REAL cs = sqrt(tpKGAMMA / (1+tpKGAMMA/(BigG-1)));
            REAL A = (sqrt(BigG-1)+cs)/(sqrt(BigG-1)-cs) * (sqrt(BigG-1)-tpA)/(sqrt(BigG-1)+tpA);
            priU[VV_XXIND] = 1.0 - 2./ (1.0 + pow(A,2./ sqrt(BigG - 1.0)));
        }

        computeConFromPri(BigG, priU, conU);


        for (int k = 0; k < VEC_DIM; k++) {
            u[i][k] = conU[k];
        }
    }

#pragma omp parallel for
    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int k = 0; k < VEC_DIM; k++)
            urk[1][i][k] = u[i][k];

#pragma omp parallel for
    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int k = 0; k < VEC_DIM; k++)
            urk[2][i][k] = u[i][k];

#pragma omp parallel for
    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int k = 0; k < VEC_DIM; k++)
            urk[3][i][k] = u[i][k];
}

void Euler1D::Boundary() {

    // pbc
    int n_real;
    for (int i = 0; i < gpn; i++) {
        n_real = i + npxx;
        for (int k = 0; k < VEC_DIM; k++) {
            u[i][k] = u[n_real][k];
        }
    }
    for (int i = npxx + gpn; i < npxx + 2 * gpn; i++) {
        n_real = i - npxx;
        for (int k = 0; k < VEC_DIM; k++) {
            u[i][k] = u[n_real][k];
        }
    }

}

void Euler1D::OutputOrder() {

    ofstream os("1DIsentropic_priU_RK3_IDEOS_nonES5_200_T08.dat", ios::ate);

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

//    int n_grid[] = {10, 20, 40, 80, 160, 320, 640, 1280};
////    int n_grid[] = {8, 16, 32, 64, 128, 256, 512, 1024};
//
//    for (int p = 0; p <= 5; p++) {
        nx = 200;

        Euler1D compute(nx);
        compute.Init();
        compute.Evolve(EC_order, ES_order);

//    }


    return 0;
}



