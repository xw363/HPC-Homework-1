# HPC Homework 1
Homework 1, High Performance Computing (NYU Spring 2015)

*(Still in progress)*

## Compile the programs
Simply enter `make` in a terminal.

## Run the programs
* For the integer version of ring communication, enter
```
mpirun -np N ./int_ring N
```
where N is a positive integer.
* For the array version of ring communication, enter
```
mpirun -np N ./array_ring N
```
where `N` is a positive integer.
* For the parallelized Jacobi solver, enter
```
mpirun -np p ./jacobi-mpi N
```
where `N` and `p` are positive intergers, and `p` divides `N`.

## Note
To change the maximum allowed number of iterations in the parallelized Jacobi
solver, modify the number in the `#define MAX_ITERATION` macro in the top
portion of `jacobi-mpi.c`.

