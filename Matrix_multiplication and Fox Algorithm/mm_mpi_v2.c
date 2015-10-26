#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define SIZE 2
#define MAX_ORDER 10000

int main(int argc, char* argv[]) {
  int             my_rank;
  int             p;
  
  
  float           local_x[MAX_ORDER];
  float           local_y[MAX_ORDER];
  float           global_x[MAX_ORDER];
  int             m, n;
  int             i,j,k;
  MPI_Datatype  mpi_column;
  float A[SIZE][SIZE];
  float B[SIZE][SIZE];
  double start,end,elapsed_time,local_elapsed;

  MPI_Init(&argc, &argv);
  /* Find out process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  /* Find out number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  m = SIZE;
  n = SIZE;

  /* MPI derived data type */
  MPI_Type_vector(SIZE/p,1,SIZE,MPI_FLOAT,&mpi_column);
  MPI_Type_commit(&mpi_column);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
        static int c = 0;
        for (i = 0; i < SIZE; i++){
            for (j = 0; j < SIZE; j++){
                A[i][j] = c++;
                B[i][j] = c++;
            }
        }
/*********************** CODE for printing matrix *******************/
if (my_rank == 0) {
      for(i=0; i<SIZE; i++){
        for(j=0; j<SIZE; j++){
          printf("A[%d][%d] = %f\t",i,j,A[i][j] );
        }
        printf("\n");
      }
      printf("\n");

    for(i=0; i<SIZE; i++){
      for(j=0; j<SIZE; j++){
        printf("B[%d][%d] = %f\t",i,j,B[i][j] );
      }
      printf("\n");
    }
  }
/*********************** CODE for printing matrix *******************/
//printf("my_rank = %d\n",my_rank );

//for(i= my_rank * SIZE/p; i< (my_rank+1) * SIZE/p; i++){
  for(i=0; i<SIZE; i++){
  //printf("i = %d and my_rank = %d\n", i,my_rank);  
  
  MPI_Allgather(&B[my_rank * SIZE/p][i], 1 , mpi_column, global_x, SIZE/p, MPI_FLOAT, MPI_COMM_WORLD);
  
    
     for(j=my_rank * SIZE/p; j<(my_rank +1 )* SIZE/p; j++){
        local_y[j] = 0;
        for(k=0; k<SIZE; k++){
          local_y[j] = local_y[j] + A[j][k] * global_x[k];
          printf("iteration %d global_x[%d] = %f and my_rank %d\n",i,k,global_x[k],my_rank);
          printf("iteration %d local_y[%d] : %f and my_rank %d\n", i,j,local_y[j] ,my_rank);
          //printf("%f\t", local_y[j]);
        }
        printf("\n");
      }

}
MPI_Barrier(MPI_COMM_WORLD); 
    end = MPI_Wtime();

    local_elapsed = end - start;
    
    MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(my_rank == 0)
      printf("Time Elapsed = %f\n", elapsed_time);
 /* Shut down MPI if p < 2 */
    MPI_Finalize();
 /* and end program */
return 0;


 }