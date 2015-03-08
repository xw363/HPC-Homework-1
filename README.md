# HPC Homework 1
Homework 1, High Performance Computing (NYU Spring 2015)

## Compile the programs
Simply enter `make` in a terminal.

## Run the programs
* For the integer version of ring communication, enter
```
mpirun -np p ./int_ring N
```
where `p` and `N` are a positive integer.
* For the array version of ring communication, enter
```
mpirun -np p ./array_ring N
```
where `p` and `N` are a positive integer.
* For the parallelized Jacobi solver, enter
```
mpirun -np p ./jacobi-mpi N num_iteration
```
where `N`, `p`, `num_iteration` are positive intergers, and `p` divides `N`.

##Miscellaneous
The file `Answer.txt` contains answers to non-coding questions.
