
#include "ConVPriVEntrV.h"

using namespace std;

// 二分法求解方程 f(x) = 0在区间（a,b）中的根：
REAL Hybrid_Root_p(REAL D, REAL mx, REAL my, REAL E, REAL BigG) {
    REAL msq = mx * mx + my * my;
    if (D <= 0 || E <= sqrt(pow(D,2) + msq)) {
        std::cout << "state not in G." << std::endl;
    }
    REAL TOL = 1e-15;
    REAL a = 0, b = 2.* E / 3.0; // EOS #2 二分法上界改变值
    REAL x = b / 2.0;
    REAL k = 0; // 当前迭代步数
    REAL N = log(b*1e15) / log(2) + 1; // 二分法预估迭代步数
    REAL f0 = 1.0 + TOL;
    REAL f, temp;
    while (f0 > TOL && k < N) {
        temp = sqrt(1 - msq/(E+x)/(E+x));
        f = x - E + msq/(E+x) + sqrt(4.* x * x + D * D * temp * temp);
        if (f < 0) {
            a = x;
        }
        else {
            b = x;
        }
        x = (a + b) / 2.0;
        f0 = fabs(f);
        k += 1;
    }
    if (x <= 0) {
        cout << "negative root for pressure!" << endl;
    }
    if (isnan(x)) {
        cout << "invalid root for pressure!" << endl;
    }
    return x;
}

// Con2Pri：由守恒变量计算原始变量
int computePriFromCon(REAL BigG, REAL *conU, REAL *priU) {
    REAL D, mx, my, E;
    D = conU[RHO_IND];
    mx = conU[MM_XXIND];
    my = conU[MM_YYIND];
    E = conU[ENE_IND];

    REAL msq = mx * mx + my * my;

    REAL pre = Hybrid_Root_p(D,mx,my,E,BigG);

    REAL xx_vel = mx / (E + pre);
    REAL yy_vel = my / (E + pre);
    REAL vel_sq = xx_vel * xx_vel + yy_vel * yy_vel;
    REAL rho = D * sqrt(1 - vel_sq);

    if (pre <= 0 || rho <= 0) {
        std::cout << "negative density or pressure : " << conU[RHO_IND] << ", " << pre;
        return 1;
    }

    priU[RHO_IND] = rho;

    priU[MM_XXIND] = xx_vel;

    priU[MM_YYIND] = yy_vel;

    priU[PRE_IND] = pre;

    return 0;
}

// Pri2Con：由原始变量计算守恒变量
int computeConFromPri(REAL BigG, REAL *priU, REAL *conU) {

    REAL rho = priU[RHO_IND];

    REAL xx_vel = priU[VV_XXIND];

    REAL yy_vel = priU[VV_YYIND];
    REAL pre = priU[PRE_IND];


    REAL vel_sq = xx_vel * xx_vel + yy_vel * yy_vel;

    REAL w_Lor = 1 / sqrt(1 - vel_sq);

    REAL theta = pre / rho; // 9/25: EOS #1
    REAL e = theta - 1.0 + sqrt(1.0 + 4.*theta * theta); // EOS #2

    REAL h = 1 + e + pre / rho;


    conU[RHO_IND] = rho * w_Lor;
    conU[MM_XXIND] = rho * h * pow(w_Lor,2) * xx_vel;
    conU[MM_YYIND] = rho * h * pow(w_Lor,2) * yy_vel;
    conU[ENE_IND] = rho * h * pow(w_Lor,2) - pre;

    return 0;
}

// 用原始变量表示熵变量
int computeEntrFromPri(REAL BigG, REAL *priU, REAL *Entr) {
    REAL rho = priU[RHO_IND];
    REAL xx_vel = priU[VV_XXIND];
    REAL yy_vel = priU[VV_YYIND];
    REAL pre = priU[PRE_IND];
    REAL vel_sq = xx_vel * xx_vel + yy_vel * yy_vel;
    REAL w_Lor = 1 / sqrt(1 - vel_sq);
    REAL beta = rho / pre; // z_3, z_1 = rho, z_2 = w_Lor * v, z_3 = rho / pre
    REAL z_1 = rho;
    REAL z_2 = w_Lor * xx_vel;
    REAL z_3 = w_Lor * yy_vel;

    REAL theta = pre / rho;
    REAL h = 2.*theta + sqrt(1.0 + 4.*theta * theta);
    REAL SEntr = log(pre) - 2.*log(rho) + 2.* log(h); // 9/27: EOS #2

    Entr[RHO_IND] = beta * h - SEntr;
    Entr[MM_XXIND] = z_2 * beta;
    Entr[MM_YYIND] = z_3 * beta;
    Entr[ENE_IND] = -w_Lor * beta;

    return 0;
}
