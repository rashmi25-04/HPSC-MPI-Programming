/* ping pong test mpi */
/* Rashmi Oak */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#define SIZE_COUNT 23
#define	ITERATIONS 20

int main(int argc,char* argv[]){
	int 		my_rank;
	int			p;
	int			tag = 0;
	double 		start,end;
	int 		count = 1;
	double 		elapsed_times[100];
	double 		Message_size[100];
	int 		k = 0;
	int 		size = 1;
	MPI_Status	status;
	int 		i,j;
	char   my_name[MPI_MAX_PROCESSOR_NAME];   /* name of machine               */
	int    my_name_len;                       /* length of returned my_name    */  

	/* Start up MPI */
MPI_Init(&argc, &argv);

/* Find out process rank */
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

/* Find out number of processes */
MPI_Comm_size(MPI_COMM_WORLD, &p);

MPI_Get_processor_name( my_name, &my_name_len );
printf( "Rank %i is running on %s\n", my_rank, my_name ); 


for(i=0; i<SIZE_COUNT; i++){

	double *msg = (double *)malloc(count*sizeof(char *));
	//printf("count = %d\n",count );
	//MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0){
	/* Take starttime */
		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();
		/* Timed iterations */
		for(j=0; j<ITERATIONS; j++){
		MPI_Send(msg, count, MPI_DOUBLE, 1, tag,MPI_COMM_WORLD);
		MPI_Recv(msg,count, MPI_DOUBLE, 1, tag,MPI_COMM_WORLD, &status);
		}
	
		/* Take endtime */
		end = MPI_Wtime();
		/* delta time = endtime - starttime */
		elapsed_times[k] = end - start;
		Message_size[k]	= count;

	}
	if(my_rank == 1){
			
			MPI_Barrier(MPI_COMM_WORLD);
			for(j=0; j<ITERATIONS; j++){
			MPI_Recv(msg, count, MPI_DOUBLE, 0, tag,MPI_COMM_WORLD, &status);
			MPI_Send(msg, count, MPI_DOUBLE, 0, tag,MPI_COMM_WORLD);
			}

	}
count = count * 2;
k = k + 1;
free(msg);

} /* end of SIZE_COUNT loop */

if(my_rank == 0){
	for(i=0; i< SIZE_COUNT; i++){
		printf("%d : Message Size is : %f and elapsed_time is : %f\n",i,Message_size[i],elapsed_times[i] );
	}
	
}

/* Shut down MPI */
MPI_Finalize();
return 0;

} /* End of Loop */