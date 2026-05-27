
#include <cmath>
#include "EulerEigValVec.h"

// 比较平均意义下左右两个状态的特征值，并返回较大值：
REAL SetEigScalingAVE(const REAL *priL, const REAL *priR, REAL &BigG,
                      REAL (*MatR)[VEC_DIM], REAL *EigVec) {

    REAL MatL[VEC_DIM][VEC_DIM], priU[VEC_DIM];
    REAL max_ave_eig = SetEigVecScalingXX(priL, priR, BigG, MatR);


    REAL max_eig1 = MaxEigValue(BigG, priL);
    REAL max_eig2 = MaxEigValue(BigG, priR);

    if (max_eig1 > max_ave_eig)
        max_ave_eig = max_eig1;

    if (max_eig2 > max_ave_eig)
        max_ave_eig = max_eig2;


    return max_ave_eig;
}

// 设置特征向量矩阵中物理量的平均值: 返回最大的特征值
REAL SetEigVecScalingXX(const REAL *priVL, const REAL *priVR,
                        REAL &BigG, REAL (*MatR)[VEC_DIM]) {

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

    // 速度取算术平均: {{v}}_{i+1/2}
    const REAL v1 = (v1L + v1R) / 2.0;
    // 洛伦兹因子取算术平均：{{w}}_{i+1/2}
    const REAL w_Lor = (w_Lor_L + w_Lor_R) / 2.0;
    // 密度取算术平均: {{\rho}}_{i+1/2}
    const REAL d = (dL + dR) / 2.0;
    // beta取算术平均: {{\beta}}_{i+1/2}
    const REAL beta = (betaL + betaR) / 2.0;
    const REAL betasqp4 = (sqrt(betaL * betaL + 4.0) + sqrt(betaR * betaR + 4.0)) / 2.0;
    const REAL p = d / beta;

    // 密度取对数平均: {{\rho}}_{i+1/2}^{\ln}
    const REAL d_log = ave_log(dL, dR);
    // beta取对数平均: {{\beta}}_{i+1/2}^{\ln}
    const REAL beta_log = ave_log(betaL, betaR);
    const REAL p_log = d_log / beta_log;
    const REAL betasqp4p2_log = ave_log(2.0 + sqrt(betaL * betaL + 4.0), 2.0 + sqrt(betaR * betaR + 4.0));
    const REAL h_log = beta / betasqp4 * (1.0 - 2./betasqp4p2_log) + 4./beta_log;


    const REAL theta_log = 1./ beta_log;

    const REAL cp = 2.0 + 4.*theta_log / sqrt(1.0 + 4.*theta_log * theta_log);
    const REAL cv = 1.0 + 4.*theta_log / sqrt(1.0 + 4.*theta_log * theta_log);
    const REAL gas_cs = sqrt(cp / cv * theta_log / h_log);

    // 平均意义下的特征向量:
    int i = 0;
    MatR[0][i] = 1.0;
    MatR[1][i] = (v1 - gas_cs) * w_Lor * h_log;
    MatR[2][i] = (1 - v1 * gas_cs) * w_Lor * h_log;

    i = 1;
    REAL coeff = 1./ sqrt(1.0 + 4.*theta_log * theta_log); // 1+e+\rho*e_r
    MatR[0][i] = 1.0;
    MatR[1][i] = v1 * w_Lor * coeff;
    MatR[2][i] = w_Lor * coeff;

    i = 2;
    MatR[0][i] = 1.0;
    MatR[1][i] = (v1 + gas_cs) * w_Lor * h_log;
    MatR[2][i] = (1 + v1 * gas_cs) * w_Lor * h_log;

    //scaling: 每个特征向量匹配的系数 d_1,d_2,d_3
    REAL tau[VEC_DIM];
    REAL d2e_r = -d_log * (theta_log + 4.*theta_log * theta_log / sqrt(1.0 + 4.*theta_log * theta_log)); // rho^2 * e_r
    tau[0] = (d2e_r / 2.0 / (coeff - h_log) - p_log * v1 / 2.0 / gas_cs / h_log) * w_Lor;
    tau[1] = -p_log * w_Lor / (coeff - h_log);
    tau[2] = (d2e_r / 2.0 / (coeff - h_log) + p_log * v1 / 2.0 / gas_cs / h_log) * w_Lor;

    for (int m = 0; m < VEC_DIM; m++)
        tau[m] = pow(tau[m], 0.5);

    // scaling 之后的特征向量矩阵 R‘: R’= R\sqrt{diag{d_1,d_2,d_3}}
    for (int j = 0; j < VEC_DIM; j++) {
        for (int i = 0; i < VEC_DIM; i++) {
            MatR[i][j] *= tau[j];
        }
    }

    REAL max_ave_eig = (fabs(v1) + gas_cs) / (1 + fabs(v1) * gas_cs);

    return max_ave_eig;
}


// Maxmum Eigenvalue
REAL MaxEigValue(REAL &BigG, const REAL *priU) {
    REAL max_eig = 0;

#if 1
    REAL v_dir;
    REAL rho = priU[RHO_IND];
    REAL vel_xx = priU[MM_XXIND];
    REAL pre = priU[ENE_IND];

    if (pre < 0 || rho < 0) {
        std::cout << "Eig negative pressure or density: " << rho << ", " << pre;
        exit(1);
    }


    REAL theta = pre / rho;
    REAL gas_cs = sqrt(2.*theta / (4.*theta + sqrt(1.0 + 4.*theta * theta)));

    max_eig = (fabs(vel_xx) + gas_cs) / (1 + fabs(vel_xx) * gas_cs);;

#endif

#if 0
    REAL eig_value[VEC_DIM];
    AllEigValue(BigG, priU, eig_value);

    for (int i = 0; i < VEC_DIM; i++) {
        if (fabs(eig_value[i]) > max_eig)
            max_eig = fabs(eig_value[i]);
    }
#endif

    return max_eig;

}

//  求出所有特征值：
void AllEigValue(REAL &BigG, const REAL *priU, REAL *eig_value) {
    REAL v_dir;

    REAL rho = priU[RHO_IND];
    REAL vel_xx = priU[MM_XXIND];
    REAL pre = priU[ENE_IND];

    v_dir = vel_xx;

    if (pre < 0 || rho < 0) {
        std::cout << "Eig negative pressure or density: " << rho << ", " << pre;
        exit(1);
    }


    REAL theta = pre / rho;
    REAL gas_cs = sqrt(2.*theta / (4.*theta + sqrt(1.0 + 4.*theta * theta)));

    eig_value[0] = (v_dir - gas_cs) / (1 - v_dir * gas_cs);
    eig_value[1] = v_dir;
    eig_value[2] = (v_dir + gas_cs) / (1 + v_dir * gas_cs);

}

/*****************************************************************************/
REAL ave_log(const REAL aL, const REAL aR) {
    const REAL xi = aL / aR;
    const REAL f = (xi - 1) / (xi + 1);
    const REAL u = f * f;

    REAL F;
    // Taylor's expansion of F:
    if (u < 0.001)
        F = 1.0 + u / 3.0 + u * u / 5.0 + pow(u, 3.0) / 7.0 + pow(u, 4.0) / 9.0
            + pow(u, 5.0) / 11.0 + pow(u, 6.0) / 13.0;
    else
        F = 0.5 * log(xi) / f;

//    return (aL+aR)/2.0;
    return (0.5 * (aL + aR) / F);

}