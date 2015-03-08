/*
 * The integer version of MPI ring communication
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

int main(int argc, char **argv) {
    int rank, size, tag = 99;
    int origin, destination;  /* Origin and destination of message */
    long N;  /* Total number of times to pass message */
    long n;  /* Number of times to pass message in current rank */
    long message;  /* Message */
    struct timeval start, finish;  /* Times that the message communication
                                      starts and finishes */
    double total_time;  /* Total communication time */
    char* ptr;  /* Dummy pointer for strtol() */
    long i;  /* Dummy index */
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    /* Set number of times to pass message in current rank */
    n = N / size;
    if (rank < N % size)
        n++;

    /* Set origin and destination of message */
    origin = (rank + size - 1) % size;
    destination = (rank + 1) % size;

    /* Set the first message */
    if (rank == 0)
        message = 0;

    /* Mark start time */
    gettimeofday(&start, NULL);

    /* Receive and send message */
    for (i = 0; i < n; ++i) {
        if (rank == 0 && i == 0) {
            /* Send the first message */
            MPI_Send(&message, 1, MPI_LONG, destination, tag,
                     MPI_COMM_WORLD);
        } else if (i == n - 1 && (rank + 1) % size == (N % size)) {
            /* Receive the last message */
            MPI_Recv(&message, 1, MPI_LONG, origin, tag, MPI_COMM_WORLD,
                     &status);
            printf("Rank %d received %ld from rank %d\n", rank, message,
                   origin);
        } else {
            /* Regular cases */
            MPI_Recv(&message, 1, MPI_LONG, origin, tag, MPI_COMM_WORLD,
                     &status);
            printf("Rank %d received %ld from rank %d\n", rank, message,
                   origin);

            message += rank;

            MPI_Send(&message, 1, MPI_LONG, destination, tag,
                     MPI_COMM_WORLD);
        }
    }

    /* Mark finish time and get communication time */
    gettimeofday(&finish, NULL);
    total_time = finish.tv_sec - start.tv_sec
                 + (finish.tv_usec - start.tv_usec) / 1e6;
    printf("Rank %d communication time: %.8f seconds\n", rank,
           total_time);

    MPI_Finalize();

    return 0;
}
