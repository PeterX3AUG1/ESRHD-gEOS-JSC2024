
#include <cmath>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "setFlux.h"

using namespace std;

//dsr: x-direction-flux-function from conservative variables
void setFluxXX(REAL &BigG, REAL *conU, REAL *flux) {
    REAL D, m ,E;
    D = conU[RHO_IND];
    m = conU[MM_XXIND];
    E = conU[ENE_IND];

    REAL pre = Hybrid_Root_p(D,m,E,BigG);
    REAL xx_vel = m / (E + pre);
    REAL rho = D * sqrt(1 - pow(xx_vel,2));

    flux[RHO_IND] = conU[RHO_IND] * xx_vel;        /**< rho v_x */
    flux[MM_XXIND] = conU[MM_XXIND] * xx_vel + pre;
    flux[ENE_IND] = conU[MM_XXIND];
}

void setEntropyConserFlux(REAL &BigG, REAL *priVL, REAL *priVR, REAL *flux) {
    const REAL dL = priVL[RHO_IND];
    const REAL v1L = priVL[VV_XXIND];
    const REAL pL = priVL[PRE_IND];
    const REAL w_Lor_L = 1 / sqrt(1 - pow(v1L,2));

    const REAL dR = priVR[RHO_IND];
    const REAL v1R = priVR[VV_XXIND];
    const REAL pR = priVR[PRE_IND];
    const REAL w_Lor_R = 1 / sqrt(1 - pow(v1R,2));

    const REAL betaL = dL / pL;
    const REAL betaR = dR / pR;

    const REAL v1 = (v1L + v1R) / 2.0;
    const REAL w_Lor = (w_Lor_L + w_Lor_R) / 2.0;

    const REAL d = (dL + dR) / 2.0;
    const REAL beta = (betaL + betaR) / 2.0;
    const REAL p = d / beta;

    const REAL d_log = ave_log(dL, dR);
    const REAL beta_log = ave_log(betaL, betaR);
    const REAL p_log = d_log / beta_log;

//    const REAL vsq = v1L * v1R;

    const REAL z_1 = d; // 算术平均
    const REAL z_2 = (w_Lor_L * v1L + w_Lor_R * v1R) / 2.0; // 算术平均
    const REAL z_3 = beta; // 算术平均

    const REAL betam2p3_log = ave_log(2.*betaL+3.0,2.*betaR+3.0);

    // 熵守恒通量：
    flux[FF_RHO_IND] = d_log * z_2;
    REAL e_hat = 1 + 3./beta_log - 3./betam2p3_log;
    REAL dh_hat = (p + d_log * e_hat) / (pow(w_Lor,2) - pow(z_2,2));
    flux[FF_MM_XXIND] = dh_hat * pow(z_2,2) + p;

    flux[FF_ENE_IND] = dh_hat * w_Lor * z_2;
}

void dissipation(int &ES_order, REAL &BigG, const REAL *priL, const REAL *priR,
                 const REAL *Entr1, const REAL *Entr2, const REAL *Entr3, const REAL *Entr4,
                 const REAL *Entr5, const REAL *Entr6, const REAL *Entr7, const REAL *Entr8,
                 REAL *Dissip) {

    REAL max_eig;
    REAL SMatR[VEC_DIM][VEC_DIM], EigVec[VEC_DIM];

    max_eig = SetEigScalingAVE(priL, priR, BigG, SMatR, EigVec);

    REAL Z1[VEC_DIM], Z2[VEC_DIM], Z3[VEC_DIM], Z4[VEC_DIM], Z5[VEC_DIM], Z6[VEC_DIM], Z7[VEC_DIM], Z8[VEC_DIM];

    for (int k = 0; k < VEC_DIM; k++) {
        Z1[k] = 0.0;
        Z2[k] = 0.0;
        Z3[k] = 0.0;
        Z4[k] = 0.0;
        Z5[k] = 0.0;
        Z6[k] = 0.0;
        Z7[k] = 0.0;
        Z8[k] = 0.0;

        // 对熵变量做伸缩处理：w = R'^{\top} W
        for (int j = 0; j < VEC_DIM; j++) {
            Z1[k] += SMatR[j][k] * Entr1[j]; //i-3
            Z2[k] += SMatR[j][k] * Entr2[j]; //i-2
            Z3[k] += SMatR[j][k] * Entr3[j]; //i-1
            Z4[k] += SMatR[j][k] * Entr4[j]; //i

            Z5[k] += SMatR[j][k] * Entr5[j]; //i+1
            Z6[k] += SMatR[j][k] * Entr6[j]; //i+2
            Z7[k] += SMatR[j][k] * Entr7[j]; //i+3
            Z8[k] += SMatR[j][k] * Entr8[j]; //i+4
        }
    }

    REAL ZLk, ZRk, DeltaZ[VEC_DIM];
    for (int k = 0; k < VEC_DIM; k++) {
        if (ES_order == 1) {
            ZLk = Z4[k];
            ZRk = Z5[k];
        } else if (ES_order == 2) {
            ZLk = Z4[k] + 0.5 * Minmod(Z4[k] - Z3[k], Z5[k] - Z4[k]);
            ZRk = Z5[k] - 0.5 * Minmod(Z6[k] - Z5[k], Z5[k] - Z4[k]);
        } else if (ES_order == 4) {
            ZLk = ENO4UR(Z1[k], Z2[k], Z3[k], Z4[k], Z5[k], Z6[k], Z7[k]);
            ZRk = ENO4UL(Z2[k], Z3[k], Z4[k], Z5[k], Z6[k], Z7[k], Z8[k]);
        } else if (ES_order == 5) {
            ZLk = WENOPW44(Z2[k], Z3[k], Z4[k], Z5[k], Z6[k]);
            ZRk = WENOPW41(Z3[k], Z4[k], Z5[k], Z6[k], Z7[k]);
        }

        DeltaZ[k] = ZRk - ZLk;

        if (ES_order == 5) {
            const REAL DZ1st = Z5[k] - Z4[k];
            if (DeltaZ[k] * DZ1st <= 0) {
                DeltaZ[k] = 0.0;
            }
        }//ES_order
    }


    for (int k = 0; k < VEC_DIM; k++) {
        Dissip[k] = 0.0;
        for (int j = 0; j < VEC_DIM; j++)
            Dissip[k] += SMatR[k][j] * DeltaZ[j];

        Dissip[k] *= 0.5 * max_eig;
    }
}

void setLFFlux(REAL &BigG, REAL *conVL, REAL *conVR, REAL *flux) {
    REAL fluxL[VEC_DIM], fluxR[VEC_DIM], priVL[VEC_DIM], priVR[VEC_DIM];
    setFluxXX(BigG, conVL, fluxL);
    setFluxXX(BigG, conVR, fluxR);

    computePriFromCon(BigG, conVL, priVL);
    computePriFromCon(BigG, conVR, priVR);

    REAL rho_L = priVL[RHO_IND], rho_R = priVR[RHO_IND];
    REAL pre_L = priVL[ENE_IND], pre_R = priVR[ENE_IND];
    REAL theta_L = pre_L / rho_L, theta_R = pre_R / rho_R;
    REAL e_L = 3.*theta_L * (3.*theta_L + 1.0) / (3.*theta_L + 2.0), e_R = 3.*theta_R * (3.*theta_R + 1.0) / (3.*theta_R + 2.0);
    REAL h_L = 1 + e_L + pre_L / rho_L, h_R = 1 + e_R + pre_R / rho_R;

    REAL gas_csL = sqrt(theta_L * (3.*theta_L + 2.0) * (18.*theta_L * theta_L + 24.*theta_L + 5.0)
                        / 3.0 / (6.*theta_L * theta_L + 4.*theta_L + 1.0) / (9.*theta_L *theta_L + 12.*theta_L + 2.0));
    REAL gas_csR = sqrt(theta_R * (3.*theta_R + 2.0) * (18.*theta_R * theta_R + 24.*theta_R + 5.0)
                        / 3.0 / (6.*theta_R * theta_R + 4.*theta_R + 1.0) / (9.*theta_R *theta_R + 12.*theta_R + 2.0));

    REAL eigL = (fabs(priVL[VV_XXIND]) + gas_csL) / (1 + fabs(priVL[VV_XXIND]) * gas_csL);
    REAL eigR = (fabs(priVR[VV_XXIND]) + gas_csR) / (1 + fabs(priVR[VV_XXIND]) * gas_csR);
    REAL max_eig = fmax(eigL, eigR);



    for (int m = 0; m < VEC_DIM; m++) {
        flux[m] = 0.5 * ((fluxR[m] + fluxL[m]) - max_eig * (conVR[m] - conVL[m]));
    }
}

