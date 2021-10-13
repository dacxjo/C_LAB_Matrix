#include "linear_func.h"
#include <math.h>
#include <stdlib.h>

int resoltrisup(int n, double **A, double *b, double *x, double tol) {
    for (int i = 0; i < n; ++i) {
        if (fabs(A[i][i]) < tol || A[i][i] == 0) {
            return 1;
        } else {
            continue;
        }
    }
    x[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j <= n - 1; j++) {
            sum += (A[i][j] * x[j]);
        }
        x[i] += (b[i] - sum) / A[i][i];
    }
    return 0;
}

int resoltriinf(int n, double **A, double *b, double *x, double tol) {
    x[0] = b[0] / A[0][0];
    for (int i = 1; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j <= (i-1); ++j) {
            sum += (A[i][j] * x[j]);
        }
        x[i] += (b[i] - sum) / A[i][i];
    }
    return 0;
}

_Bool checkupptri(double **A, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j < i && A[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

_Bool checklowertri(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (A[i][j] != 0) {
                return 0;
            }
        }
    }
    return 1;
}

double prod_esc(int n, double *x, double *y) {
    double result = 0.0;
    for (int i = 0; i < n; ++i) {
        result += x[i] * y[i];
    }
    return result;
}

double *prodMatVect(double **M, double *x, int n) {
    double *v = (double *) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        v[i] = 0;
        for (int j = 0; j < n; j++) {
            v[i] += M[i][j] * x[j];
        }
    }
    return v;
}