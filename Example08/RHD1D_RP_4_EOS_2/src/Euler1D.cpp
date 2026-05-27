
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Euler1D.h"

using namespace std;

Euler1D::Euler1D(int nx) : npxx(nx), gpn(4) {

    for (int i = 0; i < 4; i++) {
        urk[i] = Arr(npxx + 2 * gpn, VEC_DIM);
    }

//    u = Arr(npxx + 2 * gpn, VEC_DIM);

    upri = Arr(npxx + 2 * gpn, VEC_DIM);
    flux = Arr(npxx + 2 * gpn + 1, VEC_DIM);
    uEntr = Arr(npxx + 2 * gpn, VEC_DIM);

    fluxES = Arr(npxx + 2 * gpn + 1, VEC_DIM);
    fluxLF = Arr(npxx + 2 * gpn + 1, VEC_DIM);
}

Euler1D::~Euler1D() {
}

void Euler1D::Evolve(int &EC_order, int &ES_order) {

    int rk_ite, rk_to_update;
    int small_ite = 0;
    int out_flag = 0;
    REAL dt, lambda, finte, conU[VEC_DIM];

    t = 0.0;

    t_ite = 0;

    REAL max_eig;
    while (t_ite < stop_ite && small_ite < 1e8) {

        small_ite++;

        out_flag = 0;

        u = urk[0];

        max_eig = Global_MaxEig();


        dt = cfl * dxx / max_eig;
//        dt = pow(dt, 5.0 / RK_ORDER);


        if (t + dt > stop_t[t_ite]) {
            dt = stop_t[t_ite++] - t;
            out_flag = 1;
        }

        lambda = dt / dxx;

        for (rk_ite = 0; rk_ite < RK_ORDER; rk_ite++) {

            rk_to_update = (rk_ite + 1) % RK_ORDER;

            for (int i = 0; i < npxx + 2 * gpn; i++) {

                for (int k = 0; k < VEC_DIM; k++)
                    conU[k] = u[i][k];

                if (computePriFromCon(BigG, conU, upri[i])) {
                    std::cout << " rk_ite=" << rk_ite << std::endl;
                    cout << "fail to get con from pri in WENOLimit:cell average" << endl;
                    cout << "cell = " << i << endl;
                    exit(1);
                }

                computeEntrFromPri(BigG, upri[i], uEntr[i]);

            }

            //flux
            ComputeFlux(EC_order, ES_order);

            for (int i = gpn; i < npxx + gpn; i++) {
                for (int k = 0; k < VEC_DIM; k++) {
                    finte = -flux[i][k] + flux[i - 1][k];

#if RK_ORDER == 4

                    if (rk_ite < 3) {
                        urk[rk_to_update][i][k] = urk[0][i][k] + four_coe[rk_ite] * lambda * finte;
                    } else {
                        urk[rk_to_update][i][k] = (-urk[0][i][k] + urk[1][i][k] + 2 * urk[2][i][k] + urk[3][i][k] +
                                                   0.5 * lambda * finte) * 1.0 / 3;
                    }

#else
                    urk[rk_to_update][i][k] =
                            urk[0][i][k] * rk_coe[rk_ite][0] + rk_coe[rk_ite][1] * (u[i][k] + lambda * finte);

#endif

                }

            }


            u = urk[rk_to_update];


            Boundary();

        }//RK-time-step

        t = t + dt;
        cout << "t = " << t << ", dt = " << dt << endl;
//        std::cout << "t=" << t << "," << "dt=" << dt << " " << u[npxx / 2 + gpn][RHO_IND] << std::endl;


        if (out_flag) {
            OutputOrder();
        }

    }//while

}

void Euler1D::ComputeFlux(int &EC_order, int &ES_order) {

    //Flux high order Entropy Conservation and LF
    ComputFluxESAndLF(EC_order, ES_order);


    //flux couple
    for (int i = gpn - 1; i < npxx + gpn; i++) {
        for (int k = 0; k < VEC_DIM; k++) {
            flux[i][k] = fluxES[i][k];
//            flux[i][k] = fluxLF[i][k];
        }
    }
}

void Euler1D::ComputFluxESAndLF(int &EC_order, int &ES_order) {

    REAL conU_left[VEC_DIM], conU_right[VEC_DIM];
    REAL priU_left[VEC_DIM], priU_right[VEC_DIM];
    REAL Entr_left[VEC_DIM], Entr_right[VEC_DIM];


    REAL ESfiip1[VEC_DIM], ESfim1ip1[VEC_DIM], ESfiip2[VEC_DIM];

    REAL ESfim2ip1[VEC_DIM], ESfim1ip2[VEC_DIM], ESfiip3[VEC_DIM];


    for (int i = gpn - 1; i < npxx + gpn; i++) {

        // compute fiip1
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i][k];
            priU_right[k] = upri[i + 1][k];
            Entr_left[k] = uEntr[i][k];
            Entr_right[k] = uEntr[i + 1][k];
        }
        setEntropyConserFlux(BigG, priU_left, priU_right, ESfiip1);


#if DISSIPATION
        // dissipation
        REAL Entr1[VEC_DIM], Entr2[VEC_DIM], Entr3[VEC_DIM], Entr4[VEC_DIM];
        REAL Entr5[VEC_DIM], Entr6[VEC_DIM], Entr7[VEC_DIM], Entr8[VEC_DIM];

        REAL Dissip[VEC_DIM] = {0.0};

        for (int k = 0; k < VEC_DIM; k++) {
            Entr1[k] = uEntr[i - 3][k];
            Entr2[k] = uEntr[i - 2][k];
            Entr3[k] = uEntr[i - 1][k];
            Entr4[k] = uEntr[i][k];
            Entr5[k] = uEntr[i + 1][k];
            Entr6[k] = uEntr[i + 2][k];
            Entr7[k] = uEntr[i + 3][k];
            Entr8[k] = uEntr[i + 4][k];
        }

        dissipation(ES_order, BigG, priU_left, priU_right,
                    Entr1, Entr2, Entr3, Entr4,
                    Entr5, Entr6, Entr7, Entr8,
                    Dissip);

#endif

        // compute fim1ip1
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i - 1][k];
            priU_right[k] = upri[i + 1][k];
            Entr_left[k] = uEntr[i - 1][k];
            Entr_right[k] = uEntr[i + 1][k];
        }

        setEntropyConserFlux(BigG, priU_left, priU_right, ESfim1ip1);

        // compute fiip2
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i][k];
            priU_right[k] = upri[i + 2][k];
            Entr_left[k] = uEntr[i][k];
            Entr_right[k] = uEntr[i + 2][k];
        }

        setEntropyConserFlux(BigG, priU_left, priU_right, ESfiip2);

        // compute fim2ip1
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i - 2][k];
            priU_right[k] = upri[i + 1][k];
            Entr_left[k] = uEntr[i - 2][k];
            Entr_right[k] = uEntr[i + 1][k];
        }
        setEntropyConserFlux(BigG, priU_left, priU_right, ESfim2ip1);

        // compute fim1ip2
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i - 1][k];
            priU_right[k] = upri[i + 2][k];
            Entr_left[k] = uEntr[i - 1][k];
            Entr_right[k] = uEntr[i + 2][k];
        }
        setEntropyConserFlux(BigG, priU_left, priU_right, ESfim1ip2);

        // compute fiip3
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i][k];
            priU_right[k] = upri[i + 3][k];
            Entr_left[k] = uEntr[i][k];
            Entr_right[k] = uEntr[i + 3][k];
        }
        setEntropyConserFlux(BigG, priU_left, priU_right, ESfiip3);

        //flux
        for (int k = 0; k < VEC_DIM; k++) {
            if (EC_order == 2) {
                fluxES[i][k] = ESfiip1[k];
            } else if (EC_order == 4) {
                fluxES[i][k] = (4.0 * ESfiip1[k] / 3.0) - ((ESfim1ip1[k] + ESfiip2[k]) / 6.0);
            } else if (EC_order == 6) {
                fluxES[i][k] = (1.5 * ESfiip1[k]) - (0.3 * (ESfim1ip1[k] + ESfiip2[k]))
                               + (1.0 / 30.0) * (ESfim2ip1[k] + ESfim1ip2[k] + ESfiip3[k]);
            }
#if DISSIPATION
            fluxES[i][k] -= Dissip[k];
#endif
        }//


        //LF
        REAL flux_LF[VEC_DIM];
        for (int k = 0; k < VEC_DIM; k++) {
            conU_left[k] = u[i][k];
            conU_right[k] = u[i + 1][k];
            priU_left[k] = upri[i][k];
            priU_right[k] = upri[i + 1][k];
        }
        setLFFlux(BigG, conU_left, conU_right, flux_LF);

        for (int k = 0; k < VEC_DIM; k++) {
            fluxLF[i][k] = flux_LF[k];
        }

    } // i


}

REAL Euler1D::Global_MaxEig() {
    int flag = 1;

    REAL conU[VEC_DIM], priU[VEC_DIM];

    REAL max_eig = 1.0E-15;

    for (int i = 0; i < npxx + 2 * gpn; i++) {
        for (int k = 0; k < VEC_DIM; k++) {
            conU[k] = u[i][k];
        }
        computePriFromCon(BigG, conU, priU);
        max_eig = max(max_eig, MaxEigValue(BigG, priU));
    }

    return max_eig;
}

void Euler1D::ComAlert() {

}
