
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define REPS 1048576
#define TRIALS 10

int main(argc, argv)
int argc;
char **argv;
{

int     rank, size, i, j;
double  t1, t2, in=1.0, out, ttime=0.0;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

if (rank == 0) {
  

  for (j=0; j<TRIALS; j++) {
    MPI_Barrier(MPI_COMM_WORLD);   /* synchronize all tasks */

    t1 = MPI_Wtime();
    for (i=0; i<REPS; i++) 
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    t2 = MPI_Wtime();

    ttime = ttime + (t2-t1);
    printf("Trial %d: time= %f\n", j, (t2-t1)/REPS);
    }
  printf("\nAverage time= %f seconds\n", ttime/(REPS*TRIALS));
  }

else {
  for (j=0; j<TRIALS; j++) {
    MPI_Barrier(MPI_COMM_WORLD);   /* synchronize all tasks */
    for (i=0; i<REPS; i++)
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }

MPI_Finalize();
}
