
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "Euler1D.h"

using namespace std;

Euler1D::Euler1D(int nx, int ny) : npxx(nx), npyy(ny), gpn(4) {

    urk = Arr(4, npxx + 2 * gpn, npyy + 2 * gpn, VEC_DIM); // 10/10, RK_ORDER --> 4

//    u = Arr(npxx + 2 * gpn, VEC_DIM);

    upri = Arr(npxx + 2 * gpn, npyy + 2 * gpn, VEC_DIM);
    flux_xx = Arr(npxx + 2 * gpn + 1, npyy + 2 * gpn + 1, VEC_DIM);
    flux_yy = Arr(npxx + 2 * gpn + 1, npyy + 2 * gpn + 1, VEC_DIM);
    uEntr = Arr(npxx + 2 * gpn, npyy + 2* gpn, VEC_DIM);

    fluxES_xx = Arr(npxx + 2 * gpn + 1, npyy + 2 * gpn + 1, VEC_DIM);
    fluxES_yy = Arr(npxx + 2 * gpn + 1, npyy + 2 * gpn + 1, VEC_DIM);
    fluxLF_xx = Arr(npxx + 2 * gpn + 1, npyy + 2 * gpn + 1, VEC_DIM);
    fluxLF_yy = Arr(npxx + 2 * gpn + 1, npyy + 2 * gpn + 1, VEC_DIM);

    err.open("tolEntr_2DSmooth80_EOS_3_EC6.txt", ios::ate);
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

    REAL max_eig_xx, max_eig_yy;

    REAL S;

    while (t_ite < stop_ite && small_ite < 1e8) {

        S = 0.0;

        small_ite++;

        out_flag = 0;

        u = urk[0];

        max_eig_xx = Global_MaxEig(XX_FLAG);
        max_eig_yy = Global_MaxEig(YY_FLAG);

        dt = cfl / (max_eig_xx / dxx + max_eig_yy / dyy);
        dt = pow(dt, 6.0 / RK_ORDER);


        if (t + dt > stop_t[t_ite]) {
            dt = stop_t[t_ite++] - t;
            out_flag = 1;
        }

        lambda = dt / (dxx * dyy);

        for (rk_ite = 0; rk_ite < RK_ORDER; rk_ite++) {

            rk_to_update = (rk_ite + 1) % RK_ORDER;

            for (int i = 0; i < npxx + 2 * gpn; i++) {
                for (int j = 0; j < npyy + 2 * gpn; j++) {
                    for (int k = 0; k < VEC_DIM; k++)
                        conU[k] = u[i][j][k];

                    if (computePriFromCon(BigG, conU, upri[i][j])) {
                        std::cout << " rk_ite=" << rk_ite << std::endl;
                        cout << "fail to get con from pri in WENOLimit:cell average" << endl;
                        cout << "cell = (" << i << ", " << j << ")" << endl;
                        exit(1);
                    }

                    if (rk_ite == 0 && i >= gpn && i < npxx + gpn && j >= gpn && j < npyy + gpn) {
                        REAL theta = upri[i][j][PRE_IND] / upri[i][j][RHO_IND];
                        REAL h = 2.5 * theta + sqrt(1.0 + 9.*theta * theta / 4.0);
                        REAL SEntr = 1.5 * log(upri[i][j][PRE_IND]) - 2.5 * log(upri[i][j][RHO_IND]) + 1.5 * log(h - theta);
                        S += -upri[i][j][RHO_IND]/ sqrt(1-pow(upri[i][j][VV_XXIND],2)- pow(upri[i][j][VV_YYIND],2))*SEntr;
                    }

                    computeEntrFromPri(BigG, upri[i][j], uEntr[i][j]);

                }
            }

            //flux
            ComputeFlux(EC_order, ES_order);

            for (int i = gpn; i < npxx + gpn; i++) {
                for (int j = gpn; j < npyy + gpn; j++) {
                    for (int k = 0; k < VEC_DIM; k++) {
                        finte = (-flux_xx[i][j][k] + flux_xx[i - 1][j][k]) * dyy + (-flux_yy[i][j][k] + flux_yy[i][j - 1][k]) * dxx;


#if RK_ORDER == 4

                    if (rk_ite < 3) {
                        urk[rk_to_update][i][j][k] = urk[0][i][j][k] + four_coe[rk_ite] * lambda * finte;
                    } else {
                        urk[rk_to_update][i][j][k] = (-urk[0][i][j][k] + urk[1][i][j][k] + 2 * urk[2][i][j][k] + urk[3][i][j][k] +
                                                   0.5 * lambda * finte) * 1.0 / 3;
                    }

#else
                    urk[rk_to_update][i][j][k] =
                            urk[0][i][j][k] * rk_coe[rk_ite][0] + rk_coe[rk_ite][1] * (u[i][j][k] + lambda * finte);

#endif
                    }
                }

            }


            u = urk[rk_to_update];


            Boundary();

        }//RK-time-step

        S = S * dxx * dyy;

//        t = t + dt;
        std::cout << "t=" << t << "," << " " << u[npxx / 2 + gpn][npyy / 2 + gpn][RHO_IND] << ", tolEntr = " << S << std::endl;
        err << std::setprecision(16) << scientific;
        err <<  t <<  " " << S << "\n";
        t = t + dt;

        if (out_flag) {
            OutputOrder();
        }

    }//while

}

void Euler1D::ComputeFlux(int &EC_order, int &ES_order) {

    //Flux high order Entropy Conservation and LF
    ComputFluxESAndLFXX(EC_order, ES_order);
    ComputFluxESAndLFYY(EC_order, ES_order);


    //flux couple
    for (int i = gpn - 1; i < npxx + gpn; i++) {
        for (int j = gpn; j < npyy + gpn; j++) {
            for (int k = 0; k < VEC_DIM; k++) {
                flux_xx[i][j][k] = fluxES_xx[i][j][k];
//                flux_xx[i][j][k] = fluxLF_xx[i][j][k];
            }
        }
    }

    for (int j = gpn - 1; j < npyy + gpn; j++) {
        for (int i = gpn; i < npxx + gpn; i++) {
            for (int k = 0; k < VEC_DIM; k++) {
                flux_yy[i][j][k] = fluxES_yy[i][j][k];
//                flux_yy[i][j][k] = fluxLF_yy[i][j][k];
            }
        }
    }
}

void Euler1D::ComputFluxESAndLFXX(int &EC_order, int &ES_order) {

    REAL conU_left[VEC_DIM], conU_right[VEC_DIM];
    REAL priU_left[VEC_DIM], priU_right[VEC_DIM];
    REAL Entr_left[VEC_DIM], Entr_right[VEC_DIM];


    REAL ESfiip1[VEC_DIM], ESfim1ip1[VEC_DIM], ESfiip2[VEC_DIM];

    REAL ESfim2ip1[VEC_DIM], ESfim1ip2[VEC_DIM], ESfiip3[VEC_DIM];



    for (int j = gpn; j < npyy + gpn; j++) {
        for (int i = gpn - 1; i < npxx + gpn; i++) {
            // compute fiip1
            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i + 1][j][k];
                Entr_left[k] = uEntr[i][j][k];
                Entr_right[k] = uEntr[i + 1][j][k];
            }
            setEntropyConserFluxXX(BigG, priU_left, priU_right, ESfiip1);



#if DISSIPATION
        // dissipation
        REAL Entr1[VEC_DIM], Entr2[VEC_DIM], Entr3[VEC_DIM], Entr4[VEC_DIM];
        REAL Entr5[VEC_DIM], Entr6[VEC_DIM], Entr7[VEC_DIM], Entr8[VEC_DIM];

        REAL Dissip[VEC_DIM] = {0.0};


            for (int k = 0; k < VEC_DIM; k++) {
                Entr1[k] = uEntr[i - 3][j][k];
                Entr2[k] = uEntr[i - 2][j][k];
                Entr3[k] = uEntr[i - 1][j][k];
                Entr4[k] = uEntr[i][j][k];
                Entr5[k] = uEntr[i + 1][j][k];
                Entr6[k] = uEntr[i + 2][j][k];
                Entr7[k] = uEntr[i + 3][j][k];
                Entr8[k] = uEntr[i + 4][j][k];
            }

            dissipation(XX_FLAG,ES_order, BigG, priU_left, priU_right,
                        Entr1, Entr2, Entr3, Entr4,
                        Entr5, Entr6, Entr7, Entr8,
                        Dissip);


#endif

        // compute fim1ip1

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i - 1][j][k];
                priU_right[k] = upri[i + 1][j][k];
                Entr_left[k] = uEntr[i - 1][j][k];
                Entr_right[k] = uEntr[i + 1][j][k];
            }
            setEntropyConserFluxXX(BigG, priU_left, priU_right, ESfim1ip1);



        // compute fiip2

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i + 2][j][k];
                Entr_left[k] = uEntr[i][j][k];
                Entr_right[k] = uEntr[i + 2][j][k];
            }
            setEntropyConserFluxXX(BigG, priU_left, priU_right, ESfiip2);



        // compute fim2ip1

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i - 2][j][k];
                priU_right[k] = upri[i + 1][j][k];
                Entr_left[k] = uEntr[i - 2][j][k];
                Entr_right[k] = uEntr[i + 1][j][k];
            }
            setEntropyConserFluxXX(BigG, priU_left, priU_right, ESfim2ip1);



        // compute fim1ip2

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i - 1][j][k];
                priU_right[k] = upri[i + 2][j][k];
                Entr_left[k] = uEntr[i - 1][j][k];
                Entr_right[k] = uEntr[i + 2][j][k];
            }
            setEntropyConserFluxXX(BigG, priU_left, priU_right, ESfim1ip2);



        // compute fiip3

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i + 3][j][k];
                Entr_left[k] = uEntr[i][j][k];
                Entr_right[k] = uEntr[i + 3][j][k];
            }
            setEntropyConserFluxXX(BigG, priU_left, priU_right, ESfiip3);



        //flux

            for (int k = 0; k < VEC_DIM; k++) {
                if (EC_order == 2) {
                    fluxES_xx[i][j][k] = ESfiip1[k];
                } else if (EC_order == 4) {
                    fluxES_xx[i][j][k] = (4.0 * ESfiip1[k] / 3.0) - ((ESfim1ip1[k] + ESfiip2[k]) / 6.0);
                } else if (EC_order == 6) {
                    fluxES_xx[i][j][k] = (1.5 * ESfiip1[k]) - (0.3 * (ESfim1ip1[k] + ESfiip2[k]))
                                   + (1.0 / 30.0) * (ESfim2ip1[k] + ESfim1ip2[k] + ESfiip3[k]);
                }
#if DISSIPATION
                fluxES_xx[i][j][k] -= Dissip[k];
#endif
            }//
        } // i


        //LF
        REAL flux_LF[VEC_DIM];
        for (int i = gpn - 1; i < npxx + gpn; i++){
            for (int k = 0; k < VEC_DIM; k++) {
                conU_left[k] = u[i][j][k];
                conU_right[k] = u[i + 1][j][k];
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i + 1][j][k];
            }
            setLFFluxXX(BigG, conU_left, conU_right, flux_LF);

            for (int k = 0; k < VEC_DIM; k++) {
                fluxLF_xx[i][j][k] = flux_LF[k];
            }
        }

    } // j


}


void Euler1D::ComputFluxESAndLFYY(int &EC_order, int &ES_order) {

    REAL conU_left[VEC_DIM], conU_right[VEC_DIM];
    REAL priU_left[VEC_DIM], priU_right[VEC_DIM];
    REAL Entr_left[VEC_DIM], Entr_right[VEC_DIM];


    REAL ESfiip1[VEC_DIM], ESfim1ip1[VEC_DIM], ESfiip2[VEC_DIM];

    REAL ESfim2ip1[VEC_DIM], ESfim1ip2[VEC_DIM], ESfiip3[VEC_DIM];


    for (int i = gpn; i < npxx + gpn; i++) {
        for (int j = gpn - 1; j < npyy + gpn; j++) {
            // compute fiip1
            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i][j + 1][k];
                Entr_left[k] = uEntr[i][j][k];
                Entr_right[k] = uEntr[i][j + 1][k];
            }
            setEntropyConserFluxYY(BigG, priU_left, priU_right, ESfiip1);



#if DISSIPATION
        // dissipation
        REAL Entr1[VEC_DIM], Entr2[VEC_DIM], Entr3[VEC_DIM], Entr4[VEC_DIM];
        REAL Entr5[VEC_DIM], Entr6[VEC_DIM], Entr7[VEC_DIM], Entr8[VEC_DIM];

        REAL Dissip[VEC_DIM] = {0.0};


            for (int k = 0; k < VEC_DIM; k++) {
                Entr1[k] = uEntr[i][j - 3][k];
                Entr2[k] = uEntr[i][j - 2][k];
                Entr3[k] = uEntr[i][j - 1][k];
                Entr4[k] = uEntr[i][j][k];
                Entr5[k] = uEntr[i][j + 1][k];
                Entr6[k] = uEntr[i][j + 2][k];
                Entr7[k] = uEntr[i][j + 3][k];
                Entr8[k] = uEntr[i][j + 4][k];
            }

            dissipation(0, ES_order, BigG, priU_left, priU_right,
                        Entr1, Entr2, Entr3, Entr4,
                        Entr5, Entr6, Entr7, Entr8,
                        Dissip);


#endif

        // compute fim1ip1

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j - 1][k];
                priU_right[k] = upri[i][j + 1][k];
                Entr_left[k] = uEntr[i][j - 1][k];
                Entr_right[k] = uEntr[i][j + 1][k];
            }

            setEntropyConserFluxYY(BigG, priU_left, priU_right, ESfim1ip1);


        // compute fiip2

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i][j + 2][k];
                Entr_left[k] = uEntr[i][j][k];
                Entr_right[k] = uEntr[i][j + 2][k];
            }

            setEntropyConserFluxYY(BigG, priU_left, priU_right, ESfiip2);


        // compute fim2ip1

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j - 2][k];
                priU_right[k] = upri[i][j + 1][k];
                Entr_left[k] = uEntr[i][j - 2][k];
                Entr_right[k] = uEntr[i][j + 1][k];
            }
            setEntropyConserFluxYY(BigG, priU_left, priU_right, ESfim2ip1);


        // compute fim1ip2

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j - 1][k];
                priU_right[k] = upri[i][j + 2][k];
                Entr_left[k] = uEntr[i][j - 1][k];
                Entr_right[k] = uEntr[i][j + 2][k];
            }
            setEntropyConserFluxYY(BigG, priU_left, priU_right, ESfim1ip2);


        // compute fiip3

            for (int k = 0; k < VEC_DIM; k++) {
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i][j + 3][k];
                Entr_left[k] = uEntr[i][j][k];
                Entr_right[k] = uEntr[i][j + 3][k];
            }
            setEntropyConserFluxYY(BigG, priU_left, priU_right, ESfiip3);


        //flux

            for (int k = 0; k < VEC_DIM; k++) {
                if (EC_order == 2) {
                    fluxES_yy[i][j][k] = ESfiip1[k];
                } else if (EC_order == 4) {
                    fluxES_yy[i][j][k] = (4.0 * ESfiip1[k] / 3.0) - ((ESfim1ip1[k] + ESfiip2[k]) / 6.0);
                } else if (EC_order == 6) {
                    fluxES_yy[i][j][k] = (1.5 * ESfiip1[k]) - (0.3 * (ESfim1ip1[k] + ESfiip2[k]))
                                      + (1.0 / 30.0) * (ESfim2ip1[k] + ESfim1ip2[k] + ESfiip3[k]);
                }
#if DISSIPATION
                fluxES_yy[i][j][k] -= Dissip[k];
#endif
            }//
        } // j


        //LF
        REAL flux_LF[VEC_DIM];
        for (int j = gpn - 1; j < npyy + gpn; j++){
            for (int k = 0; k < VEC_DIM; k++) {
                conU_left[k] = u[i][j][k];
                conU_right[k] = u[i][j + 1][k];
                priU_left[k] = upri[i][j][k];
                priU_right[k] = upri[i][j + 1][k];
            }
            setLFFluxYY(BigG, conU_left, conU_right, flux_LF);

            for (int k = 0; k < VEC_DIM; k++) {
                fluxLF_yy[i][j][k] = flux_LF[k];
            }
        }

    } // i


}




REAL Euler1D::Global_MaxEig(int flag) {
    REAL conU[VEC_DIM], priU[VEC_DIM];

    REAL max_eig = 1.0E-15;

    for (int i = 0; i < npxx + 2 * gpn; i++) {
        for (int j = 0; j < npyy + 2 * gpn; j++) {
            for (int k = 0; k < VEC_DIM; k++) {
                conU[k] = u[i][j][k];
            }
        }
        computePriFromCon(BigG, conU, priU);
        max_eig = max(max_eig, MaxEigValue(BigG, flag, priU));
    }

    return max_eig;
}

void Euler1D::ComAlert() {

}
