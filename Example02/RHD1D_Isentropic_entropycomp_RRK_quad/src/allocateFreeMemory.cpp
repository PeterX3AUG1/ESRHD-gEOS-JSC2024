#include "allocateFreeMemory.h"
#include <iostream>
#include <cstdlib>

DP *Arr(int m) {
    DP *a = new DP[m];
    return a;
}

DP **Arr(int m, int n) {
    DP **a = new DP *[m];
    a[0] = new DP[m * n];
    for (int i = 1; i < m; i++)
        a[i] = a[i - 1] + n;

    return a;
}

int **ArrInt(int m, int n) {
    int **a = new int *[m];
    a[0] = new int[m * n];
    for (int i = 1; i < m; i++)
        a[i] = a[i - 1] + n;

    return a;
}


DP ***Arr(int m, int n, int q) {
    int i, j;
    DP ***a = new DP **[m];
    for (i = 0; i < m; i++)
        a[i] = new DP *[n];
    a[0][0] = new DP[m * n * q];

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            a[i][j] = a[0][0] + (i * n + j) * q;;

    return a;


}

DP ****Arr(int m, int n, int p, int q) {
    DP ****a = new DP ***[m];
    int i, j, k;
    for (i = 0; i < m; i++)
        a[i] = new DP **[n];
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            a[i][j] = new DP *[p];

    a[0][0][0] = new DP[m * n * p * q];
    int inc3 = q, inc2 = p * q, inc1 = n * inc2;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            for (k = 0; k < p; k++)
                a[i][j][k] = a[0][0][0] + i * inc1 + j * inc2 + k * inc3;

    return a;
}


void FreeArr(DP *a) {
    delete[]a;
}

void FreeArr(DP **a) {
    delete[] (a[0]);
    delete[] a;
}

/**
 *
 *
 * @param a
 * @param m leading dimension a[0..m-1][][]
 */
void FreeArr(DP ***a, int m) {
    delete[](a[0][0]);
    for (int i = 0; i < m; i++)
        delete[](a[i]);

    delete[]a;
}

void FreeArr(DP ****a, int m, int n) {
    delete[](a[0][0][0]);
    int i, j;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            delete[](a[i][j]);

    for (i = 0; i < m; i++)
        delete[](a[i]);

    delete[]a;
}


int *IntArr(int m) {
    int *a = new int[m];
    return a;
}

void FreeIntArr(int *a) {
    delete[]a;
}