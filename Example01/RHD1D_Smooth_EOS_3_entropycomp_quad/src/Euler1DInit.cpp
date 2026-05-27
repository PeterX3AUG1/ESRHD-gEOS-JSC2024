
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


    BigG = 5 / 3.0Q;

    xxbl = 0.0Q, xxbr = 2.0Q * PI;

    stop_t[0] = 1.5Q;

    dxx = (xxbr - xxbl) / npxx;

    cfl = 0.4Q;

    stop_ite = 1;

    priU[VV_XXIND] = 0.2Q;
    priU[PRE_IND] = 1.0Q;

    for (int i = 0; i < npxx + gpn * 2; i++) {

        xx = xxbl + (i - gpn + 0.5Q) * dxx;

        priU[RHO_IND] = 1.0Q + 0.2Q * sinq(xx);

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

        REAL theta = pre / den;
        REAL h = 2.5Q * theta + sqrtq(1.0Q + 9.0Q*theta * theta / 4.0Q);
        S = 1.5Q * logq(pre) - 2.5Q * logq(den) + 1.5Q * logq(h - theta);
        FS += -den / sqrtq(1-vx*vx) * S;
    }
    FS = FS * dxx;
    printf("t = %.16Qf, %.16Qf, tolEntr = %.32Qf\n", t, u[npxx / 2 + gpn][RHO_IND], FS);
//    std::cout << "t=" << t << " " << u[npxx / 2 + gpn][RHO_IND] << ", tolEntr = " << FS << std::endl;
//    err << std::setprecision(16) << scientific;
//    err <<  t << " " << FS << " " << 1 << "\n";
    err.close();
}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx;
    int EC_order = 6, ES_order = 5;

    int n_grid[] = {10, 20, 40, 80, 160, 320, 640, 1280};
//    int n_grid[] = {8, 16, 32, 64, 128, 256, 512, 1024};

//    for (int p = 0; p <= 5; p++) {
    nx = n_grid[3]; // Nx = 80

    Euler1D compute(nx);
    compute.Init();
    compute.Evolve(EC_order, ES_order);

//    }


    return 0;
}



