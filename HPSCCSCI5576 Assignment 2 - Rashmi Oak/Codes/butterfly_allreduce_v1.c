
/* butterfly log(p) allreduce */
/* Rashmi Oak */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define SIZE 1048576
double *a;
int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    //double       a = 1000.0;         /* Left endpoint             */
    float       b = 2.0;         /* Right endpoint            */
    int         c = 3.0;         /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    
    int         source;    /* Process sending integral  */
    int         dest ; 
    int         flag; /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;
    float       sum,recv_sum;
    float       total = 0;
    int         count = 0;
    int         i;
    float       recv_a;
    int power_2_stage;
    int stage;
    double start,finish,elapsed_time,local_elapsed;

    
    int Ceiling_log2(int  x);

    /* Let the system do what it needs to start up MPI */
    
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD); 
    start = MPI_Wtime();
    a = (double *)malloc(SIZE) ; 
    if(my_rank ==0){
    //printf("my Rank %d and my value is %f\n",my_rank,a);
    
    }
    if(atoi(argv[1]) == 0){
        printf("******************** LOW to HIGH *****************************\n");
        for(stage=0; stage <Ceiling_log2(p); stage++){
            power_2_stage = 1 << stage;
            //for(i=0; i<p; i++){
                //printf("power_2_stage = %d\n",power_2_stage);
                flag = my_rank ^ power_2_stage;
                //printf("flag = %d\n",flag );
                
                if(flag > my_rank){
                    //printf("stage : %d : my_rank is %d and I am sending to %d\n",stage,my_rank,flag);
                    MPI_Send(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD);
                    MPI_Recv(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD, &status);
                    //printf("stage : %d my_rank is %d and I am receiving from %d and received value is %f\n",stage,my_rank,flag ,a);
                }
                else{
                    //printf("stage : %d : my_rank is %d and I am sending to %d\n",stage,my_rank,flag);
                    MPI_Recv(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD, &status);
                    MPI_Send(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD);
                    //printf("stage : %d : my_rank is %d and I am receiving from %d and received value is %f\n",stage,my_rank,flag,a);

                }

                
        }
    }
    if(atoi(argv[1]) == 1) {
        printf("******************** HIGH to LOW *****************************\n");
        for(stage=Ceiling_log2(p)-1; stage >=0; stage--){
            power_2_stage = 1 << stage;
            //for(i=0; i<p; i++){
                //printf("power_2_stage = %d\n",power_2_stage);
                flag = my_rank ^ power_2_stage;
                //printf("flag = %d\n",flag );
                
                if(flag > my_rank){
                    //printf("stage : %d : my_rank is %d and I am sending to %d\n",stage,my_rank,flag);
                    MPI_Send(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD);
                    MPI_Recv(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD, &status);
                    //printf("stage : %d my_rank is %d and I am receiving from %d and received value is %f\n",stage,my_rank,flag ,a);
                }
                else{
                    //printf("stage : %d : my_rank is %d and I am sending to %d\n",stage,my_rank,flag);
                    MPI_Recv(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD, &status);
                    MPI_Send(a, 1, MPI_DOUBLE, flag, 0, MPI_COMM_WORLD);
                    //printf("stage : %d : my_rank is %d and I am receiving from %d and received value is %f\n",stage,my_rank,flag,a);

                }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD); 
    finish = MPI_Wtime();

    local_elapsed = finish - start;
    
    MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
    printf("execution time = %f\n",elapsed_time );
    }

    MPI_Finalize();

    }

int Ceiling_log2(int  x  /* in */) {
    /* Use unsigned so that right shift will fill
     * leftmost bit with 0
     */
    unsigned temp = (unsigned) x - 1;
    int result = 0;

    while (temp != 0) {
         temp = temp >> 1;
         result = result + 1 ;
    }
    return result;
} /* Ceiling_log2 */




