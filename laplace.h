/*
 * Functions and parameters for jacobi.c and gauss_seidel.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RESIDUAL_FACTOR 1e-6
#define MAX_ITERATION 1e5

double* get_u_0(int N);
double residual(double* u, int N, double h);
