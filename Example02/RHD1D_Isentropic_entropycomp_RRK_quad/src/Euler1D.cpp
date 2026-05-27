
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
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
    uEntr1 = Arr(npxx + 2 * gpn, VEC_DIM);
    uEntr2 = Arr(npxx + 2 * gpn, VEC_DIM);
    uEntr3 = Arr(npxx + 2 * gpn, VEC_DIM);

    finte_1 = Arr(npxx + 2 * gpn, VEC_DIM);
    finte_2 = Arr(npxx + 2 * gpn, VEC_DIM);
    finte_3 = Arr(npxx + 2 * gpn, VEC_DIM);

    fluxES = Arr(npxx + 2 * gpn + 1, VEC_DIM);
    fluxLF = Arr(npxx + 2 * gpn + 1, VEC_DIM);

    err.open("tolEntr_1DIsentropic_200_RRK_EC6_T08_SO_dt_bu.txt", ios::ate);
}

Euler1D::~Euler1D() {
}

void Euler1D::Evolve(int &EC_order, int &ES_order) {

    int rk_ite, rk_to_update;
    int small_ite = 0;
    int out_flag = 0;
    REAL dt, lambda, finte, conU[VEC_DIM];

    t = 0.0Q;

    t_ite = 0;

    REAL max_eig;
    REAL gam_n = 1.0Q;
    REAL S;


    while (t_ite < stop_ite && small_ite < 1e8) {

        S = 0.0Q;

        small_ite++;

        out_flag = 0;

        u = urk[0];

        max_eig = Global_MaxEig();

        dt = cfl * dxx / max_eig;
        dt = powq(dt, 6.0Q / RK_ORDER);


        if (t + gam_n * dt > stop_t[t_ite]) {
            dt = (stop_t[t_ite++] - t) / gam_n;
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

                if (rk_ite == 0 && i >= gpn && i < npxx + gpn) {
//                    S += -log(upri[i][PRE_IND] * pow(upri[i][RHO_IND], -BigG))*upri[i][RHO_IND]/ sqrt(1-pow(upri[i][VV_XXIND],2)); // ideal EOS
                    S += -upri[i][RHO_IND]/ sqrtq(1-powq(upri[i][VV_XXIND],2))*(logq(upri[i][PRE_IND])-BigG*logq(upri[i][RHO_IND]));

                }
                computeEntrFromPri(BigG, upri[i], uEntr[i]);

#if RELAX_FLAG == 1
                if (rk_ite == 0) {
                    computeEntrFromPri(BigG, upri[i], uEntr1[i]); // w(u^n)
                }
                if (rk_ite == 1) {
                    computeEntrFromPri(BigG, upri[i], uEntr2[i]); // w(u^(1))
                }
                if (rk_ite == 2) {
                    computeEntrFromPri(BigG, upri[i], uEntr3[i]); // w(u^(2))
                }
#endif
            }

            //flux
            ComputeFlux(EC_order, ES_order);

            for (int i = gpn; i < npxx + gpn; i++) {
                for (int k = 0; k < VEC_DIM; k++) {
                    finte = -flux[i][k] + flux[i - 1][k];
                    if (rk_ite == 0) {
                        finte_1[i][k] = finte;
                        urk[rk_to_update][i][k] = urk[0][i][k] + lambda * finte;
                        if (rk_to_update != rk_ite + 1) {
                            cout << "error: rk_to_update != rk_ite + 1" << std::endl;
                            exit(1);
                        }
                    }
                    if (rk_ite == 1) {
                        finte_2[i][k] = finte;
                        urk[rk_to_update][i][k] = 0.75Q * urk[0][i][k] + 0.25Q * (u[i][k] + lambda * finte);
                        if (rk_to_update != rk_ite + 1) {
                            cout << "error: rk_to_update != rk_ite + 1" << std::endl;
                            exit(1);
                        }
                    }
                    if (rk_ite == 2) {
                        finte_3[i][k] = finte;
//                        urk[rk_to_update][i][k] =
//                            urk[0][i][k] + lambda * (finte_1[i][k] + finte_2[i][k] + 4 * finte_3[i][k]) / 6.0;
                        if (rk_to_update != rk_ite - 2) {
                            cout << "error: rk_to_update != rk_ite - 2" << std::endl;
                            exit(1);
                        }
                    }



                }
            }

            if (rk_ite == 2) {
#if RELAX_FLAG == 1
                gam_n = Find_relaxation_gamma(dt, urk[0], finte_1, finte_2, finte_3, uEntr1, uEntr2, uEntr3);
#else
                gam_n = 1;
#endif
                for (int i = gpn; i < npxx + gpn; i++) {
                    for (int k = 0; k < VEC_DIM; k++) {
                        urk[rk_to_update][i][k] =
                                urk[0][i][k] + gam_n * lambda * (finte_1[i][k] + finte_2[i][k] + 4 * finte_3[i][k]) / 6.0Q;
                    }
                }
                if (rk_to_update != rk_ite - 2) {
                    cout << "error: rk_to_update != rk_ite - 2" << std::endl;
                    exit(1);
                }
            }


            u = urk[rk_to_update];


            Boundary();

        }//RK-time-step

        S = S * dxx;

        printf("t = %.16Qf, %.16Qf, tolEntr = %.32Qf, gam_n = %.32Qf\n", t, u[npxx / 2 + gpn][RHO_IND], S, gam_n);
//        std::cout << "t=" << t << " " << u[npxx / 2 + gpn][RHO_IND] << ", tolEntr = " << std::setprecision(16) << scientific << S << ", gam_n = " << gam_n << std::endl;
//        err << std::setprecision(16) << scientific;
//        err <<  t << " " << S << " " << gam_n << "\n";
        t = t + gam_n * dt;

//        if (t < stop_t[t_ite]) {
//            out_flag = 0;
//        }

        if (out_flag) {
            OutputOrder();
        }

    }//while

}

void Euler1D::ComputeFlux(int &EC_order, int &ES_order) {

    //Flux high order Entropy Conservation and LF
    ComputFluxESAndLF(EC_order, ES_order);


    //flux couple
#pragma omp parallel for
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

#pragma omp parallel for private (priU_left,priU_right,conU_left,conU_right,Entr_left,Entr_right,ESfiip1,ESfim1ip1,ESfiip2,ESfim2ip1,ESfim1ip2,ESfiip3)
    for (int i = gpn - 1; i < npxx + gpn; i++) {

        // compute fiip1
        for (int k = 0; k < VEC_DIM; k++) {
            priU_left[k] = upri[i][k];
            priU_right[k] = upri[i + 1][k];
            Entr_left[k] = uEntr[i][k];
            Entr_right[k] = uEntr[i + 1][k];
        }
        setEntropyConserFlux(BigG, priU_left, priU_right, ESfiip1);


        // dissipation
        REAL Entr1[VEC_DIM], Entr2[VEC_DIM], Entr3[VEC_DIM], Entr4[VEC_DIM];
        REAL Entr5[VEC_DIM], Entr6[VEC_DIM], Entr7[VEC_DIM], Entr8[VEC_DIM];

//        REAL Dissip[VEC_DIM] = {0.0};
        REAL Dissip[VEC_DIM];
        for (int k = 0; k < VEC_DIM; k++) {
            Dissip[k] = 0.0;
        }

#if DISSIPATION
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
                fluxES[i][k] = (4.0Q * ESfiip1[k] / 3.0Q) - ((ESfim1ip1[k] + ESfiip2[k]) / 6.0Q);
            } else if (EC_order == 6) {
                fluxES[i][k] = (1.5Q * ESfiip1[k]) - (0.3Q * (ESfim1ip1[k] + ESfiip2[k]))
                               + (1.0Q / 30.0Q) * (ESfim2ip1[k] + ESfim1ip2[k] + ESfiip3[k]);
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

    REAL max_eig = 1.0E-32Q;

#pragma omp parallel for private (conU,priU)
    for (int i = 0; i < npxx + 2 * gpn; i++) {
        for (int k = 0; k < VEC_DIM; k++) {
            conU[k] = u[i][k];
        }
        computePriFromCon(BigG, conU, priU);
        max_eig = max(max_eig, MaxEigValue(BigG, priU));
    }

    return max_eig;
}

REAL Euler1D::Find_relaxation_gamma(REAL ddt, REAL** conU_RRK, REAL** finte_1_RRK, REAL** finte_2_RRK, REAL** finte_3_RRK,
                                    REAL** Entr1_RRK, REAL** Entr2_RRK, REAL** Entr3_RRK) {
    REAL gam, e, eta, eta_gam;
    REAL dn[npxx + 2 * gpn][VEC_DIM], conU_plus_dnx[npxx + 2 * gpn][VEC_DIM];
    REAL priU_RRK[npxx + 2 * gpn][VEC_DIM], priU_plus_dnx[npxx + 2 * gpn][VEC_DIM];

    REAL TOL = 3e-32Q;
    REAL a = 1e-5Q, b = 2.0Q;
    REAL x = (a+b) / 2.0Q;

    int k_iter = 0; // 当前迭代步数
    REAL N = logq(b*1e32) / logq(2.0Q) + 1; // 二分法预估迭代步数
    REAL f0 = 1.0 + TOL;
    REAL f;
    while (f0 > TOL && k_iter < N) {
        e = 0.0Q;
        for (int i = gpn; i < npxx + gpn; i++) {
            computePriFromCon(BigG, conU_RRK[i], priU_RRK[i]);
            for (int k = 0; k < VEC_DIM; k++) {
                dn[i][k] = ddt * (finte_1_RRK[i][k] + finte_2_RRK[i][k] + 4 * finte_3_RRK[i][k]) / 6.0 / dxx;
                conU_plus_dnx[i][k] = conU_RRK[i][k] + x * dn[i][k];
                e += (Entr1_RRK[i][k] * finte_1_RRK[i][k] + Entr2_RRK[i][k] * finte_2_RRK[i][k] + 4 * Entr3_RRK[i][k] * finte_3_RRK[i][k]) / dxx;
            }
            computePriFromCon(BigG, conU_plus_dnx[i], priU_plus_dnx[i]);
        }
        e = e * ddt * dxx / 6.0;

        eta = 0.0Q, eta_gam = 0.0Q;
        for (int i = gpn; i < npxx + gpn; i++) {
            eta += -conU_RRK[i][RHO_IND] * (logq(priU_RRK[i][PRE_IND]) - BigG * logq(priU_RRK[i][RHO_IND]));
            eta_gam += -conU_plus_dnx[i][RHO_IND] * (logq(priU_plus_dnx[i][PRE_IND]) - BigG * logq(priU_plus_dnx[i][RHO_IND]));
        }
        eta = eta * dxx;
        eta_gam = eta_gam * dxx;

        f = eta_gam - eta - x*e;
        if (f < 0) {
            a = x;
        }
        else {
            b = x;
        }
        x = (a + b) / 2.0;
        f0 = fabsq(f);
        k_iter += 1;
    }
    if (x <= 0) {
        cout << "negative root for gamma!" << endl;
    }
    if (isnanq(x)) {
        cout << "invalid root for gamma!" << endl;
    }
    return x;

}

void Euler1D::ComAlert() {

}
