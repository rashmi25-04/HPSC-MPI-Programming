/* tree based broadcast */
/* Rashmi Oak */
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <math.h>

int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    float       a;         /* Left endpoint             */
    float       b;         /* Right endpoint            */
    int         n;         /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    float       local_a;   /* Left endpoint my process  */
    float       local_b;   /* Right endpoint my process */
    int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
    float       integral;  /* Integral over my interval */
    float       total;     /* Total integral            */
    int         source;    /* Process sending integral  */
    int         dest = 0;   /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;
    double start,finish,elapsed_time,local_elapsed;


    void Get_data1(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank, 
              int p);
    void Get_data2(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank, 
              int p);
    float Trap(float local_a, float local_b, int local_n,
              float h);    /* Calculate local integral  */

    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD); 
    start = MPI_Wtime();
    int dest1 =p-1;

    if(atoi(argv[1]) == 0){
        printf("low to high rank ordering is chosen\n");
        Get_data1(&a, &b, &n, my_rank, p);
        h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */

    /* Length of each process' interval of 
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral = Trap(local_a, local_b, local_n, h);

    /* Add up the integrals calculated by each process */
    if (my_rank == 0) {
        total = integral;
        for (source = 1; source < p; source++) {
            MPI_Recv(&integral, 1, MPI_FLOAT, source, tag, 
                MPI_COMM_WORLD, &status);
            total = total + integral;
            printf("process 0 is receiving from : %d\n",source);
        }
    } else {   
        MPI_Send(&integral, 1, MPI_FLOAT, dest, 
            tag, MPI_COMM_WORLD);
        printf("process %d is sending to : %d\n",my_rank,dest);

    }

    /* Print the result */
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n", n);
        printf("of the integral from %f to %f = %f\n", a, b, total); 
    }
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
        fflush(stdin);
        printf("high to low ordering is chosen\n");
        Get_data2(&a, &b, &n, my_rank, p);
        h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */
    /* Length of each process' interval of 
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral = Trap(local_a, local_b, local_n, h);
    /* Add up the integrals calculated by each process */
    if (my_rank == p-1) {
        total = integral;
        for (source = 0; source < p-1; source++) {
            MPI_Recv(&integral, 1, MPI_FLOAT, source, tag, 
                MPI_COMM_WORLD, &status);
            total = total + integral;
            printf("process %d is receiving from : %d\n",(p-1),source);
        }
    } else {   
        MPI_Send(&integral, 1, MPI_FLOAT, p-1, 
            tag, MPI_COMM_WORLD);
        printf("process %d is sending to : %d\n",my_rank,p-1);

    }

    /* Print the result */
    if (my_rank == p-1) {
        printf("With n = %d trapezoids, our estimate\n", 
            n);
        printf("of the integral from %f to %f = %f\n", 
            a, b, total); 
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


    
} /*  main  */



int Ceiling_log2(int  x  /* in */) {
    
    unsigned temp = (unsigned) x - 1;
    int result = 0;

    while (temp != 0) {
         temp = temp >> 1;
         result = result + 1 ;
    }
    //printf("result = %d\n",result);
    return result;
} /* Ceiling_log2 */


/********************************************************************/
int I_receive(
        int   stage       /* in  */,
        int   my_rank     /* in  */, 
        int   p,
        int*  source_ptr  /* out */) {
    int   power_2_stage;

    /* 2^stage = 1 << stage */
    power_2_stage = 1 << stage;
            printf("*******************inside I_receive ***************** rank is %d 2^stage: %d\n",my_rank,power_2_stage);

    if ((power_2_stage <= my_rank) && 
            (my_rank < 2*power_2_stage)) {
                printf("inside if condition of I_receive %d  stage = %d\n", my_rank,stage );

        *source_ptr = my_rank - power_2_stage;
                printf("receive : : : : : source pointer is %d and my_rank is %d\n",*source_ptr,my_rank );

        return 1;
    } else return 0;
} /* I_receive */
int I_receive1(
        int   stage       /* in  */,
        int   my_rank     /* in  */, 
        int   p,
        int*  source_ptr  /* out */) {
    int   power_2_stage;

    /* 2^stage = 1 << stage */
    power_2_stage = 1 << stage;
        printf("*******************inside I_receive1 ***************** rank is %d 2^stage: %d\n",my_rank,power_2_stage);
    if (((p-1) - power_2_stage >= my_rank) && (my_rank > (p-1) - 2*power_2_stage)) {
        printf("inside if condition of I_receive1  %d  stage = %d\n", my_rank,stage );
        *source_ptr = my_rank + power_2_stage;
        printf("source pointer is %d and my_rank is %d\n",*source_ptr,my_rank );

        return 1;
    }

    else return 0;
} /* I_receive */
/********************************************************************/
int I_send(
        int   stage     /* in  */,
        int   my_rank   /* in  */,
        int   p         /* in  */, 
        int*  dest_ptr  /* out */) {
    int power_2_stage;
    printf("*******************inside I_send ******************** %d\n", my_rank);

    /* 2^stage = 1 << stage */
    power_2_stage = 1 << stage;  
    if (my_rank < power_2_stage ){
        *dest_ptr = my_rank + power_2_stage;
                printf("stage %d and my_rank %d and dest_ptr %d\n",stage,my_rank,*dest_ptr);

        if (*dest_ptr >= p) return 0;

        else return 1;
    } else return 0;
} /* I_send */            
/********************************************************************/
int I_send1(
        int   stage     /* in  */,
        int   my_rank   /* in  */,
        int   p         /* in  */, 
        int*  dest_ptr  /* out */) {
    int power_2_stage;
    printf("*******************inside I_send1 ********************  %d\n",my_rank);
    /* 2^stage = 1 << stage */
    power_2_stage = 1 << stage;  
    //if (my_rank > power_2_stage && my_rank < (2*power_2_stage) + 2){
    if(my_rank > (p-1) - power_2_stage){
        *dest_ptr = my_rank - power_2_stage;
        
        printf("stage %d and my_rank %d\n",stage,my_rank);
        return 1;
    } else return 0;
} /* I_send */            
/********************************************************************/
void Send(
        float  a     /* in */,
        float  b     /* in */, 
        int    n     /* in */, 
        int    dest  /* in */) {

    MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
    MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
    printf("inside send function : current process is sending to: %d\n",dest);
} /* Send */
    

/********************************************************************/
void Receive(
        float*  a_ptr  /* out */, 
        float*  b_ptr  /* out */,
        int*    n_ptr  /* out */,
        int     source /* in  */) {

    MPI_Status status;
    printf("source == %d\n",source);
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0, 
        MPI_COMM_WORLD, &status);
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1, 
        MPI_COMM_WORLD, &status);
    MPI_Recv(n_ptr, 1, MPI_INT, source, 2, 
        MPI_COMM_WORLD, &status);
    printf("inside receive function : This process is receiving from : %d\n",source);


} /* Receive */




void Get_data1(
        float*  a_ptr    /* out */,
        float*  b_ptr    /* out */,
        int*    n_ptr    /* out */,
        int     my_rank  /* in  */, 
        int     p        /* in  */) {

    int source;
    int dest;
    int stage;

    int Ceiling_log2(int  x);
    int I_receive( int stage, int my_rank, int p,int*  source_ptr);
    int I_send(int stage, int my_rank, int p, int* dest_ptr);
    void Send(float a, float b, int n, int dest);
    void Receive(float* a_ptr, float* b_ptr, int* n_ptr, int source);

    if (my_rank == 0){
        //printf("Enter a, b, and n\n");
        *a_ptr = 0;
        *b_ptr = 1;
        *n_ptr = 1024;
    } 
for (stage = 0; stage < Ceiling_log2(p); stage++)
        if (I_receive(stage, my_rank, p, &source))
            Receive(a_ptr, b_ptr, n_ptr, source);
        else if (I_send(stage, my_rank, p, &dest))
            Send(*a_ptr, *b_ptr, *n_ptr, dest);

    
    
} /* Get_data1*/
void Get_data2(
        float*  a_ptr    /* out */,
        float*  b_ptr    /* out */,
        int*    n_ptr    /* out */,
        int     my_rank  /* in  */, 
        int     p        /* in  */) {

    int source;
    int dest;
    int stage;

    int Ceiling_log2(int  x);
    int I_receive1( int stage, int my_rank, int p,int*  source_ptr);
    int I_send1(int stage, int my_rank, int p, int* dest_ptr);
    void Send(float a, float b, int n, int dest);
    void Receive(float* a_ptr, float* b_ptr, int* n_ptr, int source);
    if (my_rank == p-1) {
        printf("Enter a, b, and n\n");
        *a_ptr = 0;
        *b_ptr = 1;
        *n_ptr = 1024;
    } 
    printf("Ceiling_log2(p) : %d\n", Ceiling_log2(p));
    for (stage = 0; stage < Ceiling_log2(p); stage++)
        if (I_receive1(stage, my_rank, p, &source))
            Receive(a_ptr, b_ptr, n_ptr, source);
        else if (I_send1(stage, my_rank, p, &dest))
            Send(*a_ptr, *b_ptr, *n_ptr, dest);

    

   
} /* Get_data2*/


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
    
    return_val = x*x;
    return return_val; 
} /* f */


