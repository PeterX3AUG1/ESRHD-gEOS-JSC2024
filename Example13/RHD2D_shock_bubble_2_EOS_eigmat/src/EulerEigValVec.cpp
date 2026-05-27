
#include <cmath>
#include "EulerEigValVec.h"

// 比较平均意义下左右两个状态的特征值，并返回较大值：
REAL SetEigScalingAVEXX(const REAL *priL, const REAL *priR, REAL &BigG,
                      REAL (*MatR)[VEC_DIM], REAL *EigVec) {

    REAL MatL[VEC_DIM][VEC_DIM], priU[VEC_DIM];
    REAL max_ave_eig = SetEigVecScalingXX(priL, priR, BigG, MatR);



    REAL max_eig1 = MaxEigValue(BigG, 1, priL);
    REAL max_eig2 = MaxEigValue(BigG, 1, priR);

    if (max_eig1 > max_ave_eig)
        max_ave_eig = max_eig1;

    if (max_eig2 > max_ave_eig)
        max_ave_eig = max_eig2;


    return max_ave_eig;
}


REAL SetEigScalingAVEYY(const REAL *priL, const REAL *priR, REAL &BigG,
                        REAL (*MatR)[VEC_DIM], REAL *EigVec) {

    REAL MatL[VEC_DIM][VEC_DIM], priU[VEC_DIM];
    REAL max_ave_eig = SetEigVecScalingYY(priL, priR, BigG, MatR);



    REAL max_eig1 = MaxEigValue(BigG, 0, priL);
    REAL max_eig2 = MaxEigValue(BigG, 0, priR);

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
    const REAL v2L = priVL[VV_YYIND];
    const REAL pL = priVL[PRE_IND];
    const REAL vsqL = v1L * v1L + v2L * v2L;
    const REAL w_Lor_L = 1 / sqrt(1 - vsqL);

    const REAL dR = priVR[RHO_IND];
    const REAL v1R = priVR[VV_XXIND];
    const REAL v2R = priVR[VV_YYIND];
    const REAL pR = priVR[PRE_IND];
    const REAL vsqR = v1R * v1R + v2R * v2R;
    const REAL w_Lor_R = 1 / sqrt(1 - vsqR);

    const REAL betaL = dL / pL;
    const REAL betaR = dR / pR;

    // 速度取算术平均: {{v}}_{i+1/2}
    const REAL v1 = (v1L + v1R) / 2.0;
    const REAL v2 = (v2L + v2R) / 2.0;
    // 洛伦兹因子取算术平均：{{w}}_{i+1/2}
    const REAL w_Lor = (w_Lor_L + w_Lor_R) / 2.0;
    // 密度取算术平均: {{\rho}}_{i+1/2}
    const REAL d = (dL + dR) / 2.0;
    // beta取算术平均: {{\beta}}_{i+1/2}
    const REAL beta = (betaL + betaR) / 2.0;
    const REAL p = d / beta;

    // 密度取对数平均: {{\rho}}_{i+1/2}^{\ln}
    const REAL d_log = ave_log(dL, dR);
    // beta取对数平均: {{\beta}}_{i+1/2}^{\ln}
    const REAL beta_log = ave_log(betaL, betaR);
    const REAL p_log = d_log / beta_log;
    const REAL betam2p3_log = ave_log(2.*betaL + 3.0, 2.*betaR + 3.0);
    const REAL h_log = 1 + 4./beta_log - 3./betam2p3_log;


    const REAL theta_log = 1./ beta_log;

    const REAL cp = (36.*theta_log * theta_log + 48.*theta_log + 10.0) / pow(3.*theta_log + 2.0,2);
    const REAL cv = (27.*theta_log * theta_log + 36.*theta_log + 6.0) / pow(3.*theta_log + 2.0,2);
    const REAL gas_cs = sqrt(cp / cv * theta_log / h_log);

    // 平均意义下的特征向量:
    REAL cssq = gas_cs * gas_cs;
    REAL vsq = v1 * v1 + v2 * v2;
    int i = 0;
    REAL lam1 = (v1*(1-cssq)-gas_cs*sqrt((1-vsq)*(1-v1*v1-v2*v2*cssq))) / (1-vsq*cssq);
    MatR[0][i] = 1.0;
    MatR[1][i] = (1-v1*v1) / (1-v1*lam1) * lam1 * w_Lor * h_log;
    MatR[2][i] = v2 * w_Lor * h_log;
    MatR[3][i] = (1-v1*v1) / (1-v1*lam1) * w_Lor * h_log;

    i = 1;
    REAL coeff = 1 - pow(3.*theta_log / (3.*theta_log + 2.0),2); // 1+e+\rho*e_r
    MatR[0][i] = 1.0 / w_Lor;
    MatR[1][i] = v1*coeff;
    MatR[2][i] = v2*coeff;
    MatR[3][i] = coeff;

    i = 2;
    MatR[0][i] = w_Lor*v2;
    MatR[1][i] = 2.*h_log * v1 * w_Lor * w_Lor * v2;
    MatR[2][i] = h_log * (1+2.*v2 * v2 * w_Lor * w_Lor);
    MatR[3][i] = 2.*h_log * w_Lor * w_Lor * v2;

    i = 3;
    REAL lam4 = (v1*(1-cssq)+gas_cs*sqrt((1-vsq)*(1-v1*v1-v2*v2*cssq))) / (1-vsq*cssq);
    MatR[0][i] = 1.0;
    MatR[1][i] = (1-v1*v1) / (1-v1*lam4) * lam4 * w_Lor * h_log;
    MatR[2][i] = v2 * w_Lor * h_log;
    MatR[3][i] = (1-v1*v1) / (1-v1*lam4) * w_Lor * h_log;

    //scaling: 每个特征向量匹配的系数 d_1,d_2,d_3
    REAL tau[VEC_DIM];
    tau[2] = p_log / h_log / (1-v1 * v1) / w_Lor;
    tau[1] = p_log * pow(w_Lor,3) / (h_log - 1.0);
    REAL de_r = -(27.*theta_log * theta_log + 36.*theta_log + 6.0) / pow(3.*theta_log + 2.0,2) / beta_log; // rho * e_r
//    REAL coef = -1./ beta_log + de_r; // 1+e+r*e_r-h
    REAL Bx = d_log*w_Lor+p_log*w_Lor/(coeff-h_log)-tau[2]*w_Lor*w_Lor*v2*v2;
    REAL Cx = p_log*v1*sqrt(1-v1*v1-v2*v2*cssq)/h_log/gas_cs/(1-v1*v1) * (w_Lor*w_Lor*(1-v1*v1)-cssq*(v1*v1+w_Lor*w_Lor*v2*v2))/(w_Lor*w_Lor*(1-v1*v1))/(1-(v1*v1+v2*v2)*cssq);
    tau[0] = (Bx-Cx) / 2.0;
    tau[3] = (Bx+Cx) / 2.0;

    for (int m = 0; m < VEC_DIM; m++)
        tau[m] = pow(tau[m], 0.5);

    // scaling 之后的特征向量矩阵 R‘: R’= R\sqrt{diag{d_1,d_2,d_3}}
    for (int j = 0; j < VEC_DIM; j++) {
        for (int i = 0; i < VEC_DIM; i++) {
            MatR[i][j] *= tau[j];
        }
    }

    REAL max_ave_eig = (fabs(v1)*(1-cssq) + gas_cs* sqrt((1-vsq)*(1-(v1*v1+v2*v2*cssq)))) / (1-vsq*cssq);

    return max_ave_eig;
}


REAL SetEigVecScalingYY(const REAL *priVL, const REAL *priVR,
                        REAL &BigG, REAL (*MatR)[VEC_DIM]) {

    const REAL dL = priVL[RHO_IND];
    const REAL v1L = priVL[VV_XXIND];
    const REAL v2L = priVL[VV_YYIND];
    const REAL pL = priVL[PRE_IND];
    const REAL vsqL = v1L * v1L + v2L * v2L;
    const REAL w_Lor_L = 1 / sqrt(1 - vsqL);

    const REAL dR = priVR[RHO_IND];
    const REAL v1R = priVR[VV_XXIND];
    const REAL v2R = priVR[VV_YYIND];
    const REAL pR = priVR[PRE_IND];
    const REAL vsqR = v1R * v1R + v2R * v2R;
    const REAL w_Lor_R = 1 / sqrt(1 - vsqR);

    const REAL betaL = dL / pL;
    const REAL betaR = dR / pR;

    // 速度取算术平均: {{v}}_{i+1/2}
    const REAL v1 = (v1L + v1R) / 2.0;
    const REAL v2 = (v2L + v2R) / 2.0;
    // 洛伦兹因子取算术平均：{{w}}_{i+1/2}
    const REAL w_Lor = (w_Lor_L + w_Lor_R) / 2.0;
    // 密度取算术平均: {{\rho}}_{i+1/2}
    const REAL d = (dL + dR) / 2.0;
    // beta取算术平均: {{\beta}}_{i+1/2}
    const REAL beta = (betaL + betaR) / 2.0;
    const REAL p = d / beta;

    // 密度取对数平均: {{\rho}}_{i+1/2}^{\ln}
    const REAL d_log = ave_log(dL, dR);
    // beta取对数平均: {{\beta}}_{i+1/2}^{\ln}
    const REAL beta_log = ave_log(betaL, betaR);
    const REAL p_log = d_log / beta_log;
    const REAL betam2p3_log = ave_log(2.*betaL + 3.0, 2.*betaR + 3.0);
    const REAL h_log = 1 + 4./beta_log - 3./betam2p3_log;


    const REAL theta_log = 1./ beta_log;

    const REAL cp = (36.*theta_log * theta_log + 48.*theta_log + 10.0) / pow(3.*theta_log + 2.0,2);
    const REAL cv = (27.*theta_log * theta_log + 36.*theta_log + 6.0) / pow(3.*theta_log + 2.0,2);
    const REAL gas_cs = sqrt(cp / cv * theta_log / h_log);

    // 平均意义下的特征向量:
    REAL cssq = gas_cs * gas_cs;
    REAL vsq = v1 * v1 + v2 * v2;
    int i = 0;
    REAL lam1 = (v2*(1-cssq)-gas_cs*sqrt((1-vsq)*(1-v2*v2-v1*v1*cssq))) / (1-vsq*cssq);
    MatR[0][i] = 1.0;
    MatR[1][i] = v1 * w_Lor * h_log;
    MatR[2][i] = (1-v2*v2) / (1-v2*lam1) * lam1 * w_Lor * h_log;
    MatR[3][i] = (1-v2*v2) / (1-v2*lam1) * w_Lor * h_log;

    i = 2;
    REAL coeff = 1 - pow(3.*theta_log / (3.*theta_log + 2.0),2); // 1+e+\rho*e_r
    MatR[0][i] = 1.0 / w_Lor;
    MatR[1][i] = v1*coeff;
    MatR[2][i] = v2*coeff;
    MatR[3][i] = coeff;

    i = 1;
    MatR[0][i] = w_Lor*v1;
    MatR[1][i] = h_log * (1+2.*v1 * v1 * w_Lor * w_Lor);
    MatR[2][i] = 2.*h_log * v2 * w_Lor * w_Lor * v1;
    MatR[3][i] = 2.*h_log * w_Lor * w_Lor * v1;

    i = 3;
    REAL lam4 = (v2*(1-cssq)+gas_cs*sqrt((1-vsq)*(1-v2*v2-v1*v1*cssq))) / (1-vsq*cssq);
    MatR[0][i] = 1.0;
    MatR[1][i] = v1 * w_Lor * h_log;
    MatR[2][i] = (1-v2*v2) / (1-v2*lam4) * lam4 * w_Lor * h_log;
    MatR[3][i] = (1-v2*v2) / (1-v2*lam4) * w_Lor * h_log;

    //scaling: 每个特征向量匹配的系数 d_1,d_2,d_3
    REAL tau[VEC_DIM];
    tau[1] = p_log / h_log / (1-v2 * v2) / w_Lor;
    tau[2] = p_log * pow(w_Lor,3) / (h_log - 1.0);
    REAL de_r = -(27.*theta_log * theta_log + 36.*theta_log + 6.0) / pow(3.*theta_log + 2.0,2) / beta_log; // rho * e_r
//    REAL coef = -1./ beta_log + de_r; // 1+e+r*e_r-h
    REAL By = d_log*w_Lor+p_log*w_Lor/(coeff-h_log)-tau[1]*w_Lor*w_Lor*v1*v1;
    REAL Cy = p_log*v2*sqrt(1-v2*v2-v1*v1*cssq)/h_log/gas_cs/(1-v2*v2) * (w_Lor*w_Lor*(1-v2*v2)-cssq*(v2*v2+w_Lor*w_Lor*v1*v1))/(w_Lor*w_Lor*(1-v2*v2))/(1-(v1*v1+v2*v2)*cssq);
    tau[0] = (By-Cy) / 2.0;
    tau[3] = (By+Cy) / 2.0;


    for (int m = 0; m < VEC_DIM; m++)
        tau[m] = pow(tau[m], 0.5);

    // scaling 之后的特征向量矩阵 R‘: R’= R\sqrt{diag{d_1,d_2,d_3}}
    for (int j = 0; j < VEC_DIM; j++) {
        for (int i = 0; i < VEC_DIM; i++) {
            MatR[i][j] *= tau[j];
        }
    }

    REAL max_ave_eig = (fabs(v2)*(1-cssq) + gas_cs* sqrt((1-vsq)*(1-(v2*v2+v1*v1*cssq)))) / (1-vsq*cssq);

    return max_ave_eig;
}



// Maxmum Eigenvalue
REAL MaxEigValue(REAL &BigG, int flag, const REAL *priU) {
    REAL max_eig = 0;

#if 1
    REAL v_dir;
    REAL rho = priU[RHO_IND];
    REAL vel_xx = priU[MM_XXIND];
    REAL vel_yy = priU[MM_YYIND];
    REAL pre = priU[ENE_IND];

    if (pre < 0 || rho < 0) {
        std::cout << "Eig negative pressure or density: " << rho << ", " << pre;
        exit(1);
    }


    REAL theta = pre / rho;
    REAL gas_cs = sqrt(theta * (3.*theta + 2.0) * (18.*theta * theta + 24.*theta + 5.0)
                       / 3.0 / (6.*theta * theta + 4.*theta + 1.0) / (9.*theta * theta + 12.*theta + 2.0));
    REAL cssq = gas_cs * gas_cs;
    REAL xvsq = vel_xx * vel_xx, yvsq = vel_yy * vel_yy;
    REAL vsq = xvsq + yvsq;

    REAL p_vel, n_vel;

    if (flag == 1) {
        p_vel = vel_xx;
        n_vel = vel_yy;
    }
    if (flag == 0) {
        p_vel = vel_yy;
        n_vel = vel_xx;
    }

    max_eig = (fabs(p_vel)*(1-cssq) + gas_cs* sqrt((1-vsq)*(1-(p_vel*p_vel+n_vel*n_vel*cssq)))) / (1-vsq*cssq);

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
void AllEigValue(REAL &BigG, int flag, const REAL *priU, REAL *eig_value) {

    REAL rho = priU[RHO_IND];
    REAL vel_xx = priU[MM_XXIND];
    REAL vel_yy = priU[MM_YYIND];
    REAL pre = priU[ENE_IND];

    if (pre < 0 || rho < 0) {
        std::cout << "Eig negative pressure or density: " << rho << ", " << pre;
        exit(1);
    }


    REAL theta = pre / rho;
    REAL gas_cs = sqrt(theta * (3.*theta + 2.0) * (18.*theta * theta + 24.*theta + 5.0)
                       / 3.0 / (6.*theta * theta + 4.*theta + 1.0) / (9.*theta * theta + 12.*theta + 2.0));
    REAL cssq = gas_cs * gas_cs;
    REAL xvsq = vel_xx * vel_xx, yvsq = vel_yy * vel_yy;
    REAL vsq = xvsq + yvsq;

    REAL p_vel, n_vel;

    if (flag == 1) {
        p_vel = vel_xx;
        n_vel = vel_yy;
    }
    if (flag == 0) {
        p_vel = vel_yy;
        n_vel = vel_xx;
    }
    (fabs(p_vel)*(1-cssq) + gas_cs* sqrt((1-vsq)*(1-(p_vel*p_vel+n_vel*n_vel*cssq)))) / (1-vsq*cssq);
    eig_value[0] = (p_vel*(1-cssq) - gas_cs* sqrt((1-vsq)*(1-(p_vel*p_vel+n_vel*n_vel*cssq)))) / (1-vsq*cssq);
    eig_value[1] = p_vel;
    eig_value[2] = p_vel;
    eig_value[3] = (p_vel*(1-cssq) + gas_cs* sqrt((1-vsq)*(1-(p_vel*p_vel+n_vel*n_vel*cssq)))) / (1-vsq*cssq);

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