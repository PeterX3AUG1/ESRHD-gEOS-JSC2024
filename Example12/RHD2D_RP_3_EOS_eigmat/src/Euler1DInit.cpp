
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

    xxbl = 0.0, xxbr = 1.0;
    yybl = 0.0, yybr = 1.0;

    stop_t[0] = 0.4;

    dxx = (xxbr - xxbl) / npxx;
    dyy = (yybr - yybl) / npyy;

    cfl = 0.4;

    stop_ite = 1;

    REAL uu, vv, den, pre;

    for (int i = 0; i < npxx + gpn * 2; i++) {
        for (int j = 0; j < npyy + gpn * 2; j++) {
            xx = xxbl + (i - gpn + 0.5) * dxx;
            yy = yybl + (j - gpn + 0.5) * dyy;

            if (xx > 0.5 && yy > 0.5)
            {
                uu = 0.0;
                vv = 0.0;
                den = 0.035145216124503;
                pre = 0.162931056509027;
            }
            else if (xx < 0.5 && yy > 0.5)
            {
                uu = 0.7;
                vv = 0.0;
                den = 0.1;
                pre = 1.0;
            }
            else if (xx < 0.5 && yy < 0.5)
            {
                uu = 0.0;
                vv = 0.0;
                den = 0.5;
                pre = 1.0;
            }
            else
            {
                uu = 0.0;
                vv = 0.7;
                den = 0.1;
                pre = 1.0;
            }

            priU[0] = den;
            priU[1] = uu;
            priU[2] = vv;
            priU[3] = pre;



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

    /** outflow B.C. **/
    for (int j = gpn; j <= npyy + gpn - 1; j++) {
        for (int i = 1; i <= gpn; i++) {
            for (int k = 0; k < VEC_DIM; k++) {
                u[gpn - i][j][k] = u[gpn][j][k];
                u[npxx + gpn - 1 + i][j][k] = u[npxx + gpn - 1][j][k];
            }
        }
    }

    for (int i = 0; i < npxx + 2 * gpn; i++) {
        for (int j = 1; j <= gpn; j++) {
            for (int k = 0; k < VEC_DIM; k++) {
                u[i][gpn - j][k] = u[i][gpn][k];
                u[i][npyy + gpn - 1 + j][k] = u[i][npyy + gpn - 1][k];
            }
        }
    }

}

void Euler1D::OutputOrder() {

    ofstream os("2DRHD_RP1_numsol.dat", ios::ate);

    REAL conU[VEC_DIM], priU[VEC_DIM];
    REAL xx, yy, rho_min = 1000.0, pre_min = 1000.0;
    REAL rho, uu_x, uv_y, p;


    for (int i = gpn; i < npxx + gpn; i++) {
        for (int j = gpn; j < npyy + gpn; j++) {
            xx = xxbl + (i - gpn + 0.5) * dxx;
            yy = yybl + (j - gpn + 0.5) * dyy;
            for (int k = 0; k < VEC_DIM; k++) {
                conU[k] = u[i][j][k];
            }

            computePriFromCon(BigG, conU, priU);
            rho_min = fmin(rho_min, priU[RHO_IND]);
            pre_min = fmin(pre_min, priU[PRE_IND]);

            rho = priU[0];
            uu_x = priU[1];
            uv_y = priU[2];
            p = priU[3];

            os << std::setprecision(16) << scientific;
            os << rho << " " << uu_x << " " << uv_y << " " << p << "\n";
        }
    }
    os.close();

}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx, ny;
    int EC_order = 6, ES_order = 5;

    nx = 400;
    ny = 400;

    Euler1D compute(nx, ny);
    compute.Init();
    compute.Evolve(EC_order, ES_order);

    return 0;
}



