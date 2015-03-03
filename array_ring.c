/*
 * The array version of MPI ring communication
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

int main(int argc, char **argv) {
    int rank, size, tag = 99;
    int origin, destination;  /* Origin and destination of message */
    long* message;  /* Message in the ring */
    int N;  /* Number of times to pass the message */
    struct timeval start, finish;  /* Times that the message communication
                                      starts and finishes */
    double total_time;  /* Total communication time */
    char* ptr;  /* Dummy pointer for strtol() */
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    /* Get command line argument */
    if (argc != 2) {
        printf("You need to enter exactly one positive integer as"
               " argument\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if ((N = strtol(argv[1], &ptr, 10)) <= 0) {
        printf("Your argument should be a positive integer.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    message = malloc(size * sizeof(long));

    /* Set origin and destination of message */
    origin = (rank + size - 1) % size;
    destination = (rank + 1) % size;

    /* Receive and send message */
    if (rank == 0) {
        /*
         * Rank 0 initiates and terminates the message communication.
         * It first sends 0 to the next rank, and then wait and receive the
         * message from the last rank.
         * Timing is also done here.
         */
        message[0] = 0;
        gettimeofday(&start, NULL);
        MPI_Send(message, size, MPI_LONG, destination, tag, MPI_COMM_WORLD);
        MPI_Recv(message, size, MPI_LONG, origin, tag, MPI_COMM_WORLD,
                 &status);
        gettimeofday(&finish, NULL);
    } else {
        MPI_Recv(message, size, MPI_LONG, origin, tag, MPI_COMM_WORLD,
                 &status);
        message[rank] = message[origin] + rank;
        MPI_Send(message, size, MPI_LONG, destination, tag, MPI_COMM_WORLD);
    }

    /* Print communication time after the last message is received */
    if (rank == 0) {
        printf("Final message array:\n");

        int i;
        for (i = 0; i < size; ++i)
            printf("%ld\n", message[i]);

        total_time = finish.tv_sec - start.tv_sec
                     + (finish.tv_usec - start.tv_usec) / 1e6;
        printf("Total communication time: %.8f seconds\n", total_time);
        printf("Average communication time: %.8f seconds\n", total_time / size);
    }

    free(message);

    MPI_Finalize();

    return 0;
}

