all: int_ring array_ring
int_ring: int_ring.c
	mpicc int_ring.c -o int_ring
array_ring: array_ring.c
	mpicc array_ring.c -o array_ring
jacobi-mpi: jacobi-mpi.c
	mpicc jacobi-mpi.c -o jacobi-mpi -lm

