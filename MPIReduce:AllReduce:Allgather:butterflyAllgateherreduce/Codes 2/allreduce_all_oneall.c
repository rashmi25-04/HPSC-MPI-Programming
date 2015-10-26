
/* allreduce one to all and all to one */
/* Rashmi Oak */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    double      a = 1.0;         /* Left endpoint             */
    float       b = 2.0;         /* Right endpoint            */
    int         c = 3.0;         /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;
    double       sum,recv_sum;
    
     double start,finish,elapsed_time,local_elapsed;
    
    float sumf (double a,float b);    /* Calculate local integral  */
    int Ceiling_log2(int  x);

    /* Let the system do what it needs to start up MPI */
    
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD); 
    start = MPI_Wtime();


    printf("total number of processes %d\n",p );

    for(source = 1;source < p; source++){
        MPI_Send(&a, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);

    }
    if(my_rank == 0){
        for(source = 1;source < p; source++){
            MPI_Recv(&sum, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
            printf("%d is receiving and received a  = %f\n",my_rank,sum );
            //total = total + sum;
            a = sum;   
        }
        for(source = 1;source < p; source++){
            MPI_Send(&a, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);

            printf("total %f and my rank %d\n", a , my_rank);
        }
      
    }
    else{
        MPI_Recv(&recv_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        a = recv_sum;
        printf("a received = %f and my_rank is %d\n",a,my_rank );
    }
    MPI_Barrier(MPI_COMM_WORLD); 
    finish = MPI_Wtime();

    local_elapsed = finish - start;
    
    MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
    printf("execution time = %f and my_rank %d\n",elapsed_time,my_rank );
    }
    MPI_Finalize();

    }


float sumf(double a,float b){
    return a+b ;
}






