
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

    ofstream err("err_2DEC6_EOS_eigmat.txt", ios::app);

    REAL conU[VEC_DIM], priU[VEC_DIM];
    REAL xx, yy, rho_min = 1000.0, pre_min = 1000.0;

    REAL linf_error = 0.0;
    REAL l1_error = 0.0;
    REAL l2_error = 0.0;

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

            REAL rho_exact = 1.0 + 0.2 * sin(xx + yy - 0.4 * t);
            REAL tmp_rho = fabs(priU[RHO_IND] - rho_exact);
            linf_error = fmax(linf_error, tmp_rho);
            l1_error += dxx * dyy * tmp_rho;
            l2_error += dxx * dyy * tmp_rho * tmp_rho;
        }
    }
    l2_error = sqrt(l2_error);

    REAL linf_order = 0.0, l1_order = 0.0, l2_order = 0.0;

    linf_order = log(linf_errold / linf_error) / log(2.0);
    l1_order = log(l1_errold / l1_error) / log(2.0);
    l2_order = log(l2_errold / l2_error) / log(2.0);

    cout << scientific << "Nx=" << npxx << ", Ny=" << npyy << "\t";
    cout << setprecision(16) << "linf_error=" << linf_error << "\t";
    cout << setprecision(8) << "linf_order=" << linf_order << "\n";
    cout << setprecision(16) << "l1_error=" << l1_error << "\t";
    cout << setprecision(8) << "l1_order=" << l1_order << "\n";
    cout << setprecision(16) << "l2_error=" << l2_error << "\t";
    cout << setprecision(8) << "l2_order=" << l2_order << "\n";

    err << scientific << "Nx=" << npxx << ", Ny=" << npyy << "\t";
    err << setprecision(16) << "linf_error=" << linf_error << "\t";
    err << setprecision(8) << "linf_order=" << linf_order << "\t";
    err << setprecision(16) << "l1_error=" << l1_error << "\t";
    err << setprecision(8) << "l1_order=" << l1_order << "\t";
    err << setprecision(16) << "l2_error=" << l2_error << "\t";
    err << setprecision(8) << "l2_order=" << l2_order << "\n";

    linf_errold = linf_error;
    l1_errold = l1_error;
    l2_errold = l2_error;


}


/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx, ny;
    int EC_order = 6, ES_order = 5;

    int n_grid[] = {10, 20, 40, 80, 160, 320, 640, 1280};
//    int n_grid[] = {8, 16, 32, 64, 128, 256, 512, 1024};

    for (int p = 0; p <= 5; p++) {
        nx = n_grid[p];
        ny = n_grid[p];

        Euler1D compute(nx, ny);
        compute.Init();
        compute.Evolve(EC_order, ES_order);

    }


    return 0;
}



