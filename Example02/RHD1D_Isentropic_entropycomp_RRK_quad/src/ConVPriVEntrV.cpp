
#include "ConVPriVEntrV.h"

using namespace std;

// 二分法求解方程 f(x) = 0在区间（a,b）中的根：
REAL Hybrid_Root_p(REAL D, REAL m, REAL E, REAL BigG) {
    REAL msq = m * m;
    if (D <= 0 || E <= sqrtq(powq(D,2) + msq)) {
        std::cout << "state not in G." << std::endl;
    }
    REAL TOL = 3e-32Q;
    REAL a = 0, b = (BigG - 1) * (E - D * sqrtq(1 - msq / powq(E,2)));
    REAL x = b / 2.0Q;
    REAL k = 0; // 当前迭代步数
    REAL N = logq(b*1e32) / logq(2.0Q) + 1; // 二分法预估迭代步数
    REAL f0 = 1.0 + TOL;
    REAL f;
    while (f0 > TOL && k < N) {
        f = x / (BigG - 1) - E + msq / (E + x) + D * sqrtq(1 - msq / powq(E+x,2));
        if (f < 0) {
            a = x;
        }
        else {
            b = x;
        }
        x = (a + b) / 2.0Q;
        f0 = fabsq(f);
        k += 1;
    }
    if (x <= 0) {
        cout << "negative root for pressure!" << endl;
    }
    if (isnanq(x)) {
        cout << "invalid root for pressure!" << endl;
    }
    return x;
}

// Con2Pri：由守恒变量计算原始变量
int computePriFromCon(REAL BigG, REAL *conU, REAL *priU) {
    REAL D, m ,E;
    D = conU[RHO_IND];
    m = conU[MM_XXIND];
    E = conU[ENE_IND];

    REAL pre = Hybrid_Root_p(D,m,E,BigG);

    REAL xx_vel = m / (E + pre);
    REAL rho = D * sqrtq(1 - powq(xx_vel,2));

    if (pre <= 0 || rho <= 0) {
        printf("negative density or pressure: %.16Qf, %.16Qf\n", rho, pre);
        return 1;
    }

    priU[RHO_IND] = rho;

    priU[MM_XXIND] = xx_vel;

    priU[PRE_IND] = pre;

    return 0;
}

// Pri2Con：由原始变量计算守恒变量
int computeConFromPri(REAL BigG, REAL *priU, REAL *conU) {

    REAL rho = priU[RHO_IND];

    REAL xx_vel = priU[VV_XXIND];

    REAL pre = priU[PRE_IND];


    REAL w_Lor = 1 / sqrtq(1 - powq(xx_vel,2));

    REAL e = pre / ((BigG - 1) * rho);

    REAL h = 1 + e + pre / rho;


    conU[RHO_IND] = rho * w_Lor;
    conU[MM_XXIND] = rho * h * powq(w_Lor,2) * xx_vel;
    conU[ENE_IND] = rho * h * powq(w_Lor,2) - pre;

    return 0;
}

// 用原始变量表示熵变量
int computeEntrFromPri(REAL BigG, REAL *priU, REAL *Entr) {
    REAL rho = priU[RHO_IND];
    REAL xx_vel = priU[VV_XXIND];
    REAL pre = priU[PRE_IND];
    REAL w_Lor = 1 / sqrtq(1 - powq(xx_vel,2));
    REAL beta = rho / pre; // z_3, z_1 = rho, z_2 = w_Lor * v, z_3 = rho / pre
    REAL z_1 = rho;
    REAL z_2 = w_Lor * xx_vel;

    REAL SEntr = logq(pre * powq(rho, -BigG)); // ideal EOS


    Entr[RHO_IND] = (BigG - SEntr) / (BigG - 1.0Q) + beta;
    Entr[MM_XXIND] = z_2 * beta;
    Entr[ENE_IND] = -w_Lor * beta;

    return 0;
}
