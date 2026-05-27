
#include "ConVPriVEntrV.h"

using namespace std;


// 二分法求解方程 f(x) = 0在区间（a,b）中的根：
REAL Hybrid_Root_p(REAL D, REAL m, REAL E, REAL BigG) {
    if (D <= 0 || E <= sqrt(pow(D,2) + pow(m,2))) {
        std::cout << "state not in G, D = " << D << ", m = " << m << ", E =  " << E << std::endl;
    }
    REAL TOL = 1e-15;
    REAL a = 0, b = E / 3.0; // EOS #3 二分法上界改变值
    REAL x = b / 2.0;
    REAL k = 0; // 当前迭代步数
    REAL N = log(b*1e15) / log(2) + 1;
    REAL f0 = 1.0 + TOL;
    REAL f, temp;
    while (f0 > TOL && k < N) {
        temp = sqrt(1 - m*m/(E+x)/(E+x));
        f = 1.5 * x - E + m*m/(E+x) + 0.5 * sqrt(9.* x * x + 4.* D * D * temp * temp);
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
        cout << "negative root for pressure, p = " << x << endl;
    }
    if (isnan(x)) {
        cout << "invalid root for pressure!" << endl;
        exit(1);
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
    REAL rho = D * sqrt(1 - pow(xx_vel,2));

    if (pre <= 0 || rho <= 0) {
        std::cout << "negative density or pressure : " << conU[RHO_IND] << ", " << pre;
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

    REAL w_Lor = 1 / sqrt(1 - pow(xx_vel,2));
    REAL theta = pre / rho; // 9/25: EOS #1
    REAL e = 1.5 * theta - 1.0 + sqrt(1.0 + 9.*theta * theta / 4.0); // 9/27: EOS #3
    REAL h = 1 + e + pre / rho;

    conU[RHO_IND] = rho * w_Lor;
    conU[MM_XXIND] = rho * h * pow(w_Lor,2) * xx_vel;
    conU[ENE_IND] = rho * h * pow(w_Lor,2) - pre;

    return 0;
}

// 用原始变量表示熵变量
int computeEntrFromPri(REAL BigG, REAL *priU, REAL *Entr) {
    REAL rho = priU[RHO_IND];
    REAL xx_vel = priU[VV_XXIND];
    REAL pre = priU[PRE_IND];
    REAL w_Lor = 1 / sqrt(1 - pow(xx_vel,2));
    REAL beta = rho / pre; // z_3, z_1 = rho, z_2 = w_Lor * v, z_3 = rho / pre
    REAL z_1 = rho;
    REAL z_2 = w_Lor * xx_vel;

    REAL theta = pre / rho;
    REAL h = 2.5 * theta + sqrt(1.0 + 9.*theta * theta / 4.0);
    REAL SEntr = 1.5 * log(pre) - 2.5 * log(rho) + 1.5 * log(h - theta); // EOS #3

    Entr[RHO_IND] = beta * h - SEntr;
    Entr[MM_XXIND] = z_2 * beta;
    Entr[ENE_IND] = -w_Lor * beta;

    return 0;
}
