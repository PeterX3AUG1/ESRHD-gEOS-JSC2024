
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


    BigG = 5 / 3.0Q;

    xxbl = 0.0Q, xxbr = 2.0Q * PI;

    stop_t[0] = 1.5Q;

    dxx = (xxbr - xxbl) / npxx;

    cfl = 0.4Q;

    stop_ite = 1;

    priU[VV_XXIND] = 0.2Q;
    priU[PRE_IND] = 1.0Q;

#pragma omp parallel for private (xx,priU,conU)
    for (int i = 0; i < npxx + gpn * 2; i++) {

        xx = xxbl + (i - gpn + 0.5) * dxx;

        priU[RHO_IND] = 1.0Q + 0.2Q * sinq(xx);

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

    ofstream err("err_1DEC6_EOS_3_quad.txt", ios::app);

    REAL conU[VEC_DIM], priU[VEC_DIM];
    REAL xx, rho_min = 1000.0, pre_min = 1000.0;

    REAL linf_error = 0.0Q;
    REAL l1_error = 0.0Q;
    REAL l2_error = 0.0Q;

    for (int i = gpn; i < npxx + gpn; i++) {

        xx = xxbl + (i - gpn + 0.5) * dxx;

        for (int k = 0; k < VEC_DIM; k++) {
            conU[k] = u[i][k];
        }

        computePriFromCon(BigG, conU, priU);
        rho_min = fminq(rho_min, priU[RHO_IND]);
        pre_min = fminq(pre_min, priU[PRE_IND]);

        REAL rho_exact = 1.0Q + 0.2Q * sinq(xx - 0.2Q * t);
        REAL tmp_rho = fabsq(priU[RHO_IND] - rho_exact);
        linf_error = fmaxq(linf_error, tmp_rho);
        l1_error += dxx * tmp_rho;
        l2_error += dxx * tmp_rho * tmp_rho;

    }
    l2_error = sqrtq(l2_error);

    REAL linf_order = 0.0Q, l1_order = 0.0Q, l2_order = 0.0Q;

    linf_order = logq(linf_errold / linf_error) / logq(2.0Q);
    l1_order = logq(l1_errold / l1_error) / logq(2.0Q);
    l2_order = logq(l2_errold / l2_error) / logq(2.0Q);

    printf("Nx = %d\t", npxx);
    printf("linf_error = %.16Qf\t", linf_error);
    printf("linf_order = %.8Qf\n", linf_order);
    printf("l1_error = %.16Qf\t", l1_error);
    printf("l1_order = %.8Qf\n", l1_order);
    printf("l2_error = %.16Qf\t", l2_error);
    printf("l2_order = %.8Qf\n", l2_order);


    linf_errold = linf_error;
    l1_errold = l1_error;
    l2_errold = l2_error;


}



/******************************************************************************************************/
/******************************************************************************************************/
/******************************************************************************************************/

int main(int argc, char *argv[]) {
    int nx;
    int EC_order = 6, ES_order = 5;

    int n_grid[] = {10, 20, 40, 80, 160, 320, 640, 1280};
//    int n_grid[] = {8, 16, 32, 64, 128, 256, 512, 1024};

    for (int p = 0; p <= 5; p++) {
        nx = n_grid[p];

        Euler1D compute(nx);
        compute.Init();
        compute.Evolve(EC_order, ES_order);

    }


    return 0;
}



