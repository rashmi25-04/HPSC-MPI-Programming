/* Matrix Transpose using MPI derived data types*/
/* Rashmi Oak */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

int main(int argc,char* argv[]){

int 		my_rank; 	/* rank of the processes */
int 		p; 			/* total number of processes */
int 		tag = 0;
MPI_Status	status;   	/* return status for MPI_recv */
int 		n; 			/* take input from user */
int 		i,j,k;		/* counters */
MPI_Datatype	mpi_Column;
int 		temp;
int** 		matrix;


temp = atoi(argv[1]);
if(temp < 1){
	printf("Not a valid input\n");
}
else{
	n = temp;
}
MPI_Init(&argc, &argv);
/* Find out process rank */
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
/* Find out number of processes */
MPI_Comm_size(MPI_COMM_WORLD, &p);

/* MPI derived data type */
MPI_Type_vector(n,1,n,MPI_INT,&mpi_Column);
MPI_Type_commit(&mpi_Column);

/* Create a Matrix */
matrix = (int **)malloc(sizeof(int*) * n);
matrix[0] = (int *)malloc(sizeof(int *) * n * n);

i = 0;

for(j=0; j<n; j++){
	matrix[j] = matrix[0] + (j * n);
	if(my_rank == 1){
		for(k=0; k<n; k++){
			i = i+1;
			matrix[j][k] = i;
		}
	}
}

if(my_rank == 0){
	for(i=0; i<n; i++){
		MPI_Recv(&matrix[i][0],n,MPI_INT,1,tag,MPI_COMM_WORLD,&status);
		printf("Received Matrix is : \n");
		for(j=0; j<n; j++){
			for(k=0; k<n; k++){
				printf("%d\t", matrix[j][k]);
				
			}
			printf("\n");
		}

	}
}
if(my_rank == 1){
	printf("Original matrix is : \n");
	for(j=0; j<n; j++){
			for(k=0; k<n; k++){
				printf("%d\t", matrix[j][k]);
				
			}
			printf("\n");
		}
	for(j=0; j<n; j++){
		MPI_Send(&(matrix[0][j]),1,mpi_Column,0,tag,MPI_COMM_WORLD);
	}
}



 /* Shut down MPI if p < 2 */
MPI_Finalize();
 /* and end program */
return 0;
}