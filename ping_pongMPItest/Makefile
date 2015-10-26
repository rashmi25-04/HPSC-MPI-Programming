#Uncomment your choice of compiler
CC    =   mpicc
#CC    =   mpigcc
#CC    =   mpipgcc

FLAGS := 

all:	ping_pong	\
	transpose	\


ping_pong: ping_pong_mpi.c
	$(CC) $(FLAGS) ping_pong_mpi.c -o ping_pong

transpose: transpose.c
	$(CC) $(FLAGS) transpose.c -o transpose
