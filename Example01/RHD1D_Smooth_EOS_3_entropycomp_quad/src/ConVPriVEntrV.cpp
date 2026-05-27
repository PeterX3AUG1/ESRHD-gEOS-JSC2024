
#include "ConVPriVEntrV.h"

using namespace std;


// 二分法求解方程 f(x) = 0在区间（a,b）中的根：
REAL Hybrid_Root_p(REAL D, REAL m, REAL E, REAL BigG) {
    if (D <= 0 || E <= sqrtq(D*D +m*m)) {
        printf("state not in G, D = %.32Qf, m = %.32Qf, E = %.32Qf\n", D, m, E);
    }
    REAL TOL = 3e-32Q;
    REAL a = 0, b = E / 3.0; // EOS #3 二分法上界改变值
    REAL x = b / 2.0;
    REAL k = 0; // 当前迭代步数
    REAL N = logq(b*1e32) / logq(2.0Q) + 1; // 二分法预估迭代步数
    REAL f0 = 1.0 + TOL;
    REAL f, temp;
    while (f0 > TOL && k < N) {
        temp = sqrtq(1 - m*m/(E+x)/(E+x));
        f = 1.5 * x - E + m*m/(E+x) + 0.5 * sqrtq(9.0Q * x * x + 4.0Q * D * D * temp * temp);
        if (f < 0) {
            a = x;
        }
        else {
            b = x;
        }
        x = (a + b) / 2.0;
        f0 = fabsq(f);
        k += 1;
    }
    if (x <= 0) {
        printf("negative root for pressure, p = %.32Qf\n", x);
    }
    if (isnanq(x)) {
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
    REAL rho = D * sqrtq(1 - xx_vel * xx_vel);

    if (pre <= 0 || rho <= 0) {
        printf("negative density or pressure: %.32Qf\n", rho, pre);
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


    REAL w_Lor = 1.0Q / sqrtq(1 - xx_vel * xx_vel);

    REAL theta = pre / rho; // EOS #1
    REAL e = 1.5Q * theta - 1.0Q + sqrtq(1.0Q + 9.0Q *theta * theta / 4.0Q); // EOS #3

    REAL h = 1 + e + pre / rho;


    conU[RHO_IND] = rho * w_Lor;
    conU[MM_XXIND] = rho * h * w_Lor * w_Lor * xx_vel;
    conU[ENE_IND] = rho * h * w_Lor * w_Lor - pre;

    return 0;
}

// 用原始变量表示熵变量
int computeEntrFromPri(REAL BigG, REAL *priU, REAL *Entr) {
    REAL rho = priU[RHO_IND];
    REAL xx_vel = priU[VV_XXIND];
    REAL pre = priU[PRE_IND];
    REAL w_Lor = 1.0Q / sqrtq(1 - xx_vel * xx_vel);
    REAL beta = rho / pre; // z_3, z_1 = rho, z_2 = w_Lor * v, z_3 = rho / pre
    REAL z_1 = rho;
    REAL z_2 = w_Lor * xx_vel;

    REAL theta = pre / rho;
    REAL h = 2.5Q * theta + sqrtq(1.0Q + 9.0Q *theta * theta / 4.0Q);
    REAL SEntr = 1.5Q * logq(pre) - 2.5Q * logq(rho) + 1.5Q * logq(h - theta); // 9/27: EOS #3

    Entr[RHO_IND] = beta * h - SEntr;
    Entr[MM_XXIND] = z_2 * beta;
    Entr[ENE_IND] = -w_Lor * beta;

    return 0;
}
