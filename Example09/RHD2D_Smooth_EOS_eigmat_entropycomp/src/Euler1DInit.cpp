
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
    REAL xx, yy;

    u = urk[0];


    BigG = 5 / 3.0;

    xxbl = 0.0, xxbr = 2.0 * PI;
    yybl = 0.0, yybr = 2.0 * PI;

    stop_t[0] = 0.1;

    dxx = (xxbr - xxbl) / npxx;
    dyy = (yybr - yybl) / npyy;

    cfl = 0.4;

    stop_ite = 1;

    priU[VV_XXIND] = 0.2;
    priU[VV_YYIND] = 0.2;
    priU[PRE_IND] = 1.0;

    for (int i = 0; i < npxx + gpn * 2; i++) {
        for (int j = 0; j < npyy + gpn * 2; j++) {
            xx = xxbl + (i - gpn + 0.5) * dxx;
            yy = yybl + (j - gpn + 0.5) * dyy;

            priU[RHO_IND] = 1.0 + 0.2 * sin(xx + yy);

            computeConFromPri(BigG, priU, conU);


            for (int k = 0; k < VEC_DIM; k++) {
                u[i][j][k] = conU[k];
            }
        }
    }

    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int j = 0; j < npyy + 2 * gpn; j++)
            for (int k = 0; k < VEC_DIM; k++)
                urk[1][i][j][k] = u[i][j][k];


    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int j = 0; j < npyy + 2 * gpn; j++)
            for (int k = 0; k < VEC_DIM; k++)
                urk[2][i][j][k] = u[i][j][k];


    for (int i = 0; i < npxx + 2 * gpn; i++)
        for (int j = 0; j < npyy + 2 * gpn; j++)
            for (int k = 0; k < VEC_DIM; k++)
                urk[3][i][j][k] = u[i][j][k];

}

void Euler1D::Boundary() {

    // pbc
    /** pbc **/
    for (int j = gpn; j <= npyy + gpn - 1; j++) {
        for (int i = 1; i <= gpn; i++) {
            for (int k = 0; k < VEC_DIM; k++) {
                u[gpn - i][j][k] = u[npxx + gpn - i][j][k];
                u[npxx + gpn - 1 + i][j][k] = u[gpn + i - 1][j][k];
            }
        }
    }

    for (int i = 0; i < npxx + 2 * gpn; i++) {
        for (int j = 1; j <= gpn; j++) {
            for (int k = 0; k < VEC_DIM; k++) {
                u[i][gpn - j][k] = u[i][npyy + gpn - j][k];
                u[i][npyy + gpn + j - 1][k] = u[i][gpn + j - 1][k];
            }
        }
    }

}

void Euler1D::OutputOrder() {

    // 计算终止时刻的离散熵
    REAL FS = 0.0, S, pre, den, vx, vy;
    REAL conU[VEC_DIM];
    REAL xx, yy;

    for (int i = gpn; i < npxx + gpn; i++) {
        for (int j = gpn; j < npyy + gpn; j++) {
            xx = xxbl + (i - gpn + 0.5) * dxx;
            yy = yybl + (j - gpn + 0.5) * dyy;
            for (int k = 0; k < VEC_DIM; k++) {
                conU[k] = u[i][j][k];
            }
            if (computePriFromCon(BigG, conU, upri[i][j])) {
                cout << "fail to get con from pri in WENOLimit:cell average" << endl;
                cout << "cell = " << i << endl;
                exit(1);
            }
            den = upri[i][j][RHO_IND];
            vx = upri[i][j][VV_XXIND];
            vy = upri[i][j][VV_YYIND];
            pre = upri[i][j][PRE_IND];

            S = 1.5 * log(3. * pre + 2. * den) - 3. * den / (3. * pre + 2. * den) + 1.5 * log(pre) - 4. * log(den);
            FS += -den / sqrt(1 - vx * vx - vy * vy) * S;
        }
    }
    FS = FS * dxx * dyy;
    std::cout << "t=" << t << " " << u[npxx / 2 + gpn][npyy / 2 + gpn][RHO_IND] << ", tolEntr = " << FS << std::endl;
    err << std::setprecision(16) << scientific;
    err <<  t << " " << FS << "\n";
    err.close();
}




/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx, ny;
    int EC_order = 6, ES_order = 5;

    int n_grid[] = {10, 20, 40, 80, 160, 320, 640, 1280};
//    int n_grid[] = {8, 16, 32, 64, 128, 256, 512, 1024};

//    for (int p = 0; p <= 5; p++) {
    nx = n_grid[3]; // Nx = 80
    ny = nx;

    Euler1D compute(nx, ny);
    compute.Init();
    compute.Evolve(EC_order, ES_order);

//    }


    return 0;
}



