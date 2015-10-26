
/* custom all refduce calling custom broadcast and custom reduce*/
/* Rashmi Oak */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
double *a;
int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
            /* Left endpoint             */
    float       b = 2.0;         /* Right endpoint            */
    int         c = 3.0;         /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    
    double       sum,recv_sum;
    float total = 0;
    int count = 0;
    int arg;
    int i;
    //int i;
   //MPI_Status  status;
    double start,finish,elapsed_time,local_elapsed;
    
    //float sumf (float a,float b);    /* Calculate local integral  */
    int Ceiling_log2(int  x);
    void allreduce_broadcast(int p,double *a,int my_rank,int arg);
    void allreduce_reduce(int p,double *a,int my_rank,int arg);


    /* Let the system do what it needs to start up MPI */
    
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    MPI_Barrier(MPI_COMM_WORLD); 
    start = MPI_Wtime();
    a = (double *)malloc(512) ; 
    printf("sizeof a %lu\n", sizeof(*a)); 
    if(my_rank == 0){
        *a = 1.0;
    
    }
    if(my_rank == p-1)
    {
        *a = 3.0;
    }
    
    
    if(atoi(argv[1]) == 0){
    arg = 0;
    allreduce_reduce(p,a,my_rank,arg);
    //printf("result = %f\n",*result );
    //printf("value of a = %f and my_rank is %d\n",a,my_rank );
    //a = a + 5;
    //printf("new a %f for rank %d\n",a,my_rank );
    allreduce_broadcast(p,a,my_rank,arg);

    }
    if(atoi(argv[1]) == 1){
    arg = 1;
    allreduce_reduce(p,a,my_rank,arg);
    //printf("result = %f\n",*result );
    //printf("value of a = %f and my_rank is %d\n",a,my_rank );
    //a = a + 5;
    //printf("new a %f for rank %d\n",a,my_rank );
    allreduce_broadcast(p,a,my_rank,arg);
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


float sumf(float a,float b){
    return a+b ;
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
    //printf("result = %d\n",result);
    return result;
} /* Ceiling_log2 */

void allreduce_broadcast(int p,double *a,int my_rank,int arg){
    int Ceiling_log2(int  x);
    int stage;
    int power_2_stage;
    int source_ptr;
    int dest_ptr;
    float recv;
    MPI_Status status;
    //int i;
    double *recv_a = (double *)malloc(512);
    printf("inside allreduce broadcast  a = %f and arg = %d\n", *a,arg);
    
    
    if(arg == 0){
        printf("******** LOW to HIGH broadcast *********\n");
        for (stage = 0; stage < Ceiling_log2(p); stage++){

        power_2_stage = 1 << stage;
        //printf("*******************inside allreduce_broadcast ***************** rank is %d 2^stage: %d\n",my_rank,power_2_stage);


        if ((power_2_stage <= my_rank) && (my_rank < 2*power_2_stage)) {
            //printf("inside if condition of I_receive %d  stage = %d\n", my_rank,stage );
            source_ptr = my_rank - power_2_stage;
            //printf("receive : : : : : source pointer is %d and my_rank is %d\n",source_ptr,my_rank );
            MPI_Recv(recv_a, 1, MPI_DOUBLE, source_ptr, 0, MPI_COMM_WORLD, &status);
            a = recv_a;
            printf("stage %d received a value is %f and my_rank is %d\n",stage,*a,my_rank );
        }
        if (my_rank < power_2_stage ){
            dest_ptr = my_rank + power_2_stage;
            printf("stage %d and my_rank %d and dest_ptr %d\n",stage,my_rank,dest_ptr);
            MPI_Send(a, 1, MPI_DOUBLE, dest_ptr, 0, MPI_COMM_WORLD);

    }
   
}
//printf("a = %f value for my_rank = %d\n",a,my_rank );
}

    if(arg == 1){
        printf("******** HIGH to LOW broadcast *********\n");
        for (stage = 0; stage < Ceiling_log2(p); stage++){


        power_2_stage = 1 << stage;
        printf("*******************inside allreduce_broadcast ***************** rank is %d 2^stage: %d\n",my_rank,power_2_stage);

        if (((p-1) - power_2_stage >= my_rank) && (my_rank > (p-1) - 2*power_2_stage)) {
            printf("inside if condition of I_receive1  %d  stage = %d\n", my_rank,stage );
            source_ptr = my_rank + power_2_stage;
            printf("source pointer is %d and my_rank is %d\n",source_ptr,my_rank );
            MPI_Recv(recv_a, 1, MPI_DOUBLE, source_ptr, 0, MPI_COMM_WORLD, &status);
            a = recv_a;
            printf("stage %d received a value is %f and my_rank is %d\n",stage,*a,my_rank );


            
        }
        if(my_rank > (p-1) - power_2_stage){
            dest_ptr = my_rank - power_2_stage;
            printf("stage %d and my_rank %d\n",stage,my_rank);
            printf("stage %d and my_rank %d and dest_ptr %d\n",stage,my_rank,dest_ptr);
            MPI_Send(a, 1, MPI_DOUBLE, dest_ptr, 0, MPI_COMM_WORLD);
        }
        }
    }
printf("a = %f value for my_rank = %d\n",*a,my_rank );
//return a;
    
}

void allreduce_reduce(int p,double *a,int my_rank,int arg){
    int Ceiling_log2(int  x);
    int stage;
    int power_2_stage;
    int source_ptr;
    int dest_ptr;
    double *recv_a = (double *)malloc(512);
    MPI_Status status;
    //printf("sizeof recv_a %lu\n", );




    if(arg == 0){
        printf("******** LOW to HIGH reduce *********\n");
        for(stage=Ceiling_log2(p)-1; stage >= 0; stage--){
            printf("stage = %d\n",stage );
            power_2_stage = 1 << stage;
            printf("power_2_stage %d\n",power_2_stage );
            if(my_rank < power_2_stage){
                source_ptr = my_rank + power_2_stage;
                printf("process %d is receiving from %d\n",my_rank,source_ptr );
                //total = total + integral;
                MPI_Recv(recv_a, 1, MPI_DOUBLE, source_ptr, 0, MPI_COMM_WORLD, &status);
                a = recv_a;
                //a = 5 + recv_a;
                printf("received value is %f\n",*recv_a );;
                //total = integral + recv_integral;

            }
            if ((power_2_stage <= my_rank) && 
                (my_rank < 2*power_2_stage)) {
                dest_ptr = my_rank - power_2_stage;
                printf("process %d is sending to %d\n",my_rank,dest_ptr );
                //a = a + 5;
                MPI_Send(a, 1, MPI_DOUBLE, dest_ptr, 0, MPI_COMM_WORLD);
            }
            
            }
        
        if(my_rank == 0){
        printf("all reduce : a = %f my_rank = %d\n",*a,my_rank );
        //return a;
    }
        /* Shut down MPI */
        
    }
    if(arg == 1){
        printf("******** HIGH to LOW reduce *********\n");
        for(stage=Ceiling_log2(p)-1; stage >= 0; stage--){
            printf("stage = %d\n",stage );
            power_2_stage = 1 << stage;
            printf("power_2_stage %d\n",power_2_stage );
            if(my_rank > (p-1) - power_2_stage){
                source_ptr = my_rank - power_2_stage;
                printf("process %d is receiving from %d\n",my_rank,source_ptr );
                
                MPI_Recv(recv_a, 1, MPI_DOUBLE, source_ptr, 0, MPI_COMM_WORLD, &status);
                a = recv_a;
                printf("received value is %f\n",*recv_a);;
                //a = a + recv;

            }
            if (((p-1) - power_2_stage >= my_rank) && (my_rank > (p-1) - 2*power_2_stage)) {
                dest_ptr = my_rank + power_2_stage;
                printf("process %d is sending to %d\n",my_rank,dest_ptr );
                //a = a +5;
                MPI_Send(a, 1, MPI_DOUBLE, dest_ptr, 0, MPI_COMM_WORLD);
            }
            
            }
            if(my_rank == p-1){
         printf("all reduce : a = %f my_rank = %d\n",*a,my_rank );  
         
         //return a;
         }
        }
    
}