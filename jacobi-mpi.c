/*
 * The Jacobi method using MPI
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define RESIDUAL_FACTOR 1e-6
#define MAX_ITERATION 10
#define ROOT 0

void jacobi(int N, int rank, int p, MPI_Status *status);
double residual(double* u, int N, double h);

int main(int argc, char** argv)
{
    int N;  /* N - size of u; */
    int p, rank;  /* p - number of processes */
    char* ptr;  /* Dummy pointer for strtol() */
    MPI_Status status;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Get command line argument */
    if (argc != 2) {
        printf("You need to enter exactly one positive integer as "
               "argument\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if ((N = strtol(argv[1], &ptr, 10)) <= 0) {
        printf("Your argument should be a positive integer.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (N % p) {
        printf("Your argument should be a multiple of number of"
               " available processor cores.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Run the Jacobi method */
    jacobi(N, rank, p, &status);

    MPI_Finalize();

    return 0;
}

/* The Jacobi method */
void jacobi(int N, int rank, int p, MPI_Status *status) {
    /* Variables for the Jacobian method */
    double h = 1.0 / (N + 1);
    int n = N / p;  /* Size of sub-arrays */
    double *u;  /* The complete array */
    double *u_i, *u_i_old;  /* The sub-array of current rank
                                     and its copy */
    double u_prev = 0.0;  /* The last element of the previous sub-array */
                          /* For the first sub-array, it is always 0.0 */
    double u_next = 0.0;  /* The first element of the next sub-array */
                          /* For the last sub-array, it is always 0.0 */
    double res, res_min;  /* The residual and the minimum allowed
                             residual */
    int i;  /* Dummy index */
    int k = 0;  /* Loop counter, also used as tag for MPI communcation */
    double aui = 0.0;  /* Placeholder for sum of a_ij*u_j */

    /* Initialize u, u_i, and u_i_old */
    u = (double*) malloc(N * sizeof(double));
    u_i = (double*) malloc(n * sizeof(double));
    u_i_old = (double*) malloc(n * sizeof(double));
    for (i = 0; i < N; ++i) u[i] = 0.0;
    for (i = 0; i < n; ++i) u_i[i] = 0.0;

    /* Compute the initial residual and the minimum allowed residual */
    res = residual(u, N, h);
    res_min = res * RESIDUAL_FACTOR;

    /* The Jacobi loop */
    while (res > res_min && k < MAX_ITERATION) {
        k++;

        memcpy(u_i_old, u_i, sizeof(double) * n);

        if (p == 1) {
            /* The unparalleled version, just in case */
            for (i = 0; i < n; ++i) {
                aui = 0.0;

                if (i == 0)
                    aui = -u_i_old[1] / (h * h);
                else if (i == N - 1)
                    aui = -u_i_old[N - 2] / (h * h);
                else
                    aui = -(u_i_old[i - 1] + u_i_old[i + 1]) / (h * h);

                u_i[i] = (1 - aui) / 2 * h * h;
            }

            res = residual(u_i, N, h);
        } else {
            /* Receive values of points for stencil computation */
            if (k > 1) {
                if (rank > 0)
                    MPI_Recv(&u_prev, 1, MPI_DOUBLE, rank - 1, k - 1,
                             MPI_COMM_WORLD, status);

                if (rank < p - 1)
                    MPI_Recv(&u_next, 1, MPI_DOUBLE, rank + 1, k - 1,
                             MPI_COMM_WORLD, status);
            }

            /* Update all elements */
            for (i = 0; i < n; ++i) {
                aui = 0.0;

                if (i == 0)
                    aui = -(u_prev + u_i_old[i + 1]) / (h * h);
                else if (i == n - 1)
                    aui = -(u_i_old[i - 1] + u_next) / (h * h);
                else
                    aui = -(u_i_old[i - 1] + u_i_old[i + 1]) / (h * h);

                u_i[i] = (1 - aui) / 2 * h * h;
            }

            /* At the root rank, gather the entire vector u,
             * and broadcast the residual*/
            MPI_Gather(u_i, n, MPI_DOUBLE, u, n, MPI_DOUBLE, ROOT,
                       MPI_COMM_WORLD);

            if (rank == ROOT)
                res = residual(u, N, h);

            MPI_Bcast(&res, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

            /* Send values of points for stencil computation */
            if (res > res_min && k < MAX_ITERATION) {
                if (rank > 0)
                    MPI_Send(&u_i[0], 1, MPI_DOUBLE, rank - 1, k,
                             MPI_COMM_WORLD);

                if (rank < p - 1)
                    MPI_Send(&u_i[n - 1], 1, MPI_DOUBLE, rank + 1, k,
                             MPI_COMM_WORLD);
            }
        }

        /* Print residual at the end of each loop */
        if (rank == ROOT) printf("%.8f\n", res);
    }

    /* Post-process in the root rank */
    if (rank == ROOT) {
        if (p == 1) {
            /* For the unparalleled case, copy u_i to u,
             * as all updates are done to u_i */
            memcpy(u, u_i, sizeof(double) * N);
        }

        if (res > res_min)
            printf("Initial residual decreased by: %.8f\n",
                    1.0 - res / (res_min / RESIDUAL_FACTOR));
        printf("Number of iterations: %d\n", k);
    }

    /* Clean up */
    free(u);
    free(u_i);
    free(u_i_old);
}

/* Compute residual ||A*u-f|| */
double residual(double* u, int N, double h)
{
    double res = 0.0, tmp;
    int i;

    for (i = 0; i < N; ++i) {
        if (i == 0) {
            tmp = 1 - (2 * u[0] - u[1]) / (h * h);
            res += tmp * tmp;
        } else if (i == N - 1) {
            tmp = 1 - (2 * u[N - 1] - u[N - 2]) / (h * h);
            res += tmp * tmp;
        } else {
            tmp = 1 - (2 * u[i] - u[i - 1] - u[i + 1]) / (h * h);
            res += tmp * tmp;
        }
    }

    res = sqrt(res);

    return res;
}

