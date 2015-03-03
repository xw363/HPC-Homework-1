/*
 * Implementation of functions defined in laplace.h
 */
#include "laplace.h"

double* get_u_0(int N)
{
    double* u_0 = (double*) malloc(N * sizeof(double));

    int i;
    for (i = 0; i < N; ++i) {
        u_0[i] = 0.0;
    }

    return u_0;
}

double residual(double* u, int N, double h)
{
    double res, tmp;
    int i;

    res = 0.0;
    for (i = 0; i < N; ++i) {
        if (i == 0) {
            tmp = 1 - (2 * u[0] - u[1]) / (h * h);
            res += tmp * tmp;
        }
        else if (i == N - 1) {
            tmp = 1 - (2 * u[N - 1] - u[N - 2]) / (h * h);
            res += tmp * tmp;
        }
        else {
            tmp = 1 - (2 * u[i] - u[i - 1] - u[i + 1]) / (h * h);
            res += tmp * tmp;
        }
    }
    res = sqrt(res);

    return res;
}

