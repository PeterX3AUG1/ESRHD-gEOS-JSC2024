
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <stdarg.h>
#include <sys/stat.h>
#include <string>
#include "Euler1D.h"

REAL linf_errold = 0.0, l1_errold = 0.0, l2_errold = 0.0;

using namespace std;

void Euler1D::Init() {

    REAL conU[VEC_DIM], priU[VEC_DIM];
    REAL xx;

    u = urk[0];


    BigG = 5 / 3.0;

    xxbl = 0.0, xxbr = 1.0;

    stop_t[0] = 0.376;

    dxx = (xxbr - xxbl) / npxx;

    cfl = 0.4;

    stop_ite = 1;

    for (int i = 0; i < npxx + gpn * 2; i++) {

        xx = xxbl + (i - gpn + 0.5) * dxx;

        if (xx <= 0.5) {
            priU[RHO_IND] = 5.0;
            priU[VV_XXIND] = 0.0;
            priU[PRE_IND] = 50.0;
        }
        else {
            priU[RHO_IND] = 2.0 + 0.3 * sin(50*xx);
            priU[VV_XXIND] = 0.0;
            priU[PRE_IND] = 5.0;
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

    // 计算终止时刻的离散熵
    REAL FS = 0.0, S, pre, den, vx;
    REAL conU[VEC_DIM];
    REAL xx;

    for (int i = gpn; i < npxx + gpn; i++) {

        xx = xxbl + (i - gpn + 0.5) * dxx;

        for (int k = 0; k < VEC_DIM; k++) {
            conU[k] = u[i][k];
        }
        if (computePriFromCon(BigG, conU, upri[i])) {
            cout << "fail to get con from pri in WENOLimit:cell average" << endl;
            cout << "cell = " << i << endl;
            exit(1);
        }
        den = upri[i][RHO_IND];
        vx = upri[i][VV_XXIND];
        pre = upri[i][PRE_IND];

        S = 1.5*log(3.*pre+2.*den) - 3.*den/(3.*pre+2.*den) + 1.5*log(pre) - 4.*log(den);
        FS += -den / sqrt(1-vx*vx) * S;
    }
    FS = FS * dxx;
    std::cout << "t=" << t << " " << u[npxx / 2 + gpn][RHO_IND] << ", tolEntr = " << FS << std::endl;
    err << std::setprecision(16) << scientific;
    err <<  t << " " << FS << "\n";
    err.close();
}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx;
    int EC_order = 6, ES_order = 5;

    nx = 100;

    Euler1D compute(nx);
    compute.Init();
    compute.Evolve(EC_order, ES_order);

    return 0;
}



