
/* tree based reduce */
/* Rashmi Oak */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    float       a = 0.0;         /* Left endpoint             */
    float       b = 3.0;         /* Right endpoint            */
    int         n = 1024;         /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    float       local_a;   /* Left endpoint my process  */
    float       local_b;   /* Right endpoint my process */
    int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
    float       integral;
    float       recv_integral;  /* Integral over my interval */
    float       total = 0;
    float       total1 = 0;    /* Total integral            */
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;
    int power_2_stage;
    int         dest_ptr;
    int         source_ptr;
    int stage;
    double start,finish,elapsed_time,local_elapsed;

    
    float Trap(float local_a, float local_b, int local_n,
              float h);    /* Calculate local integral  */
    int Ceiling_log2(int  x);

    /* Let the system do what it needs to start up MPI */
    
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD); 
    start = MPI_Wtime();


    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */

    /* Length of each process' interval of 
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral = Trap(local_a, local_b, local_n, h);

    printf("my Rank %d and my integral value %f\n",my_rank,integral );
    printf("total number of processes %d\n",p );

    if(atoi(argv[1]) == 0){
        printf("******** LOW to HIGH *********\n");
        for(stage=Ceiling_log2(p)-1; stage >= 0; stage--){
            printf("stage = %d\n",stage );
            power_2_stage = 1 << stage;
            printf("power_2_stage %d\n",power_2_stage );
            if(my_rank < power_2_stage){
                source_ptr = my_rank + power_2_stage;
                printf("process %d is receiving from %d\n",my_rank,source_ptr );
                //total = total + integral;
                MPI_Recv(&recv_integral, 1, MPI_FLOAT, source_ptr, 0, MPI_COMM_WORLD, &status);
                integral = integral + recv_integral;
                //total = integral + recv_integral;

            }
            if ((power_2_stage <= my_rank) && 
                (my_rank < 2*power_2_stage)) {
                dest_ptr = my_rank - power_2_stage;
                printf("process %d is sending to %d\n",my_rank,dest_ptr );
                MPI_Send(&integral, 1, MPI_FLOAT, dest_ptr, 0, MPI_COMM_WORLD);
            }
            
            }
        
        if(my_rank == 0){
        printf("Total integral = %f\n",integral);}
        MPI_Barrier(MPI_COMM_WORLD); 
        finish = MPI_Wtime();

        local_elapsed = finish - start;
        
        MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(my_rank == 0){
        printf("execution time = %f\n",elapsed_time );
        }
        /* Shut down MPI */
        MPI_Finalize();
    }
    if(atoi(argv[1]) == 1){
        printf("******** HIGH to LOW *********\n");
        for(stage=Ceiling_log2(p)-1; stage >= 0; stage--){
            printf("stage = %d\n",stage );
            power_2_stage = 1 << stage;
            printf("power_2_stage %d\n",power_2_stage );
            if(my_rank > (p-1) - power_2_stage){
                source_ptr = my_rank - power_2_stage;
                printf("process %d is receiving from %d\n",my_rank,source_ptr );
                //total = total + integral;
                MPI_Recv(&recv_integral, 1, MPI_FLOAT, source_ptr, 0, MPI_COMM_WORLD, &status);
                integral = integral + recv_integral;
                //total = integral + recv_integral;

            }
            if (((p-1) - power_2_stage >= my_rank) && (my_rank > (p-1) - 2*power_2_stage)) {
                dest_ptr = my_rank + power_2_stage;
                printf("process %d is sending to %d\n",my_rank,dest_ptr );
                MPI_Send(&integral, 1, MPI_FLOAT, dest_ptr, 0, MPI_COMM_WORLD);
            }
            
            }
        
        if(my_rank == (p-1)){
        printf("Total integral = %f\n",integral);}
        /* Shut down MPI */
        MPI_Barrier(MPI_COMM_WORLD); 
        finish = MPI_Wtime();

        local_elapsed = finish - start;
        
        MPI_Reduce(&local_elapsed, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(my_rank == 0){
        printf("execution time = %f\n",elapsed_time );
        }
        MPI_Finalize();

    }
} /*  main  */



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



float Trap(
          float  local_a   /* in */, 
          float  local_b   /* in */, 
          int    local_n   /* in */, 
          float  h         /* in */) { 

    float integral;   /* Store result in integral  */ 
    float x; 
    int i; 

    float f(float x); /* function we're integrating */

    integral = (f(local_a) + f(local_b))/2.0; 
    x = local_a; 
    for (i = 1; i <= local_n-1; i++) { 
        x = x + h; 
        integral = integral + f(x); 
    } 
    integral = integral*h; 
    return integral;
} /*  Trap  */


/********************************************************************/
float f(float x) { 
    float return_val; 
    /* Calculate f(x). */
    /* Store calculation in return_val. */ 
    return_val = x;
    return return_val; 
} /* f */



