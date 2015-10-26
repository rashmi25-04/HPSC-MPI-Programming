#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
  
float simpson(float local_a, float local_b, float local_n, float h){
  int flag,k=1;                         
  double  sum;
  double f(double x)     ;  
  float local_h;  
  sum=f(local_a); 
  printf("f(local_a)%f\n",sum)   ;                  
  flag=4;
    while (k <= local_n-1)               
    {
                        
      sum = sum + flag*f(local_a+k*h); 
      flag=6-flag;    
      k++;                      
    }                   
    

  sum = ( sum + f(local_b) )*h/3 ;    
                                
  return sum;

 }

double f(double x)
    {
     double function;
     function = x*x;
     return(function);
    }

int main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    float       a = 0.0;   /* Left endpoint             */
    float       b = 2.0;   /* Right endpoint            */
    int         n ;  /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    float       local_a;   /* Left endpoint my process  */
    float       local_b;   /* Right endpoint my process */
    int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
    float       integral;  /* Integral over my interval */
    float       total;     /* Total integral            */
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;


    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    n = atoi(argv[1]);
    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */

    
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral = simpson(local_a, local_b, local_n, h);

    
    /* Add up the integrals calculated by each process */
    
    if (my_rank == 0) {
        total = integral;
        for (source = 1; source < p; source++) {
            MPI_Recv(&integral, 1, MPI_FLOAT, source, tag,
                MPI_COMM_WORLD, &status);
            total = total + integral;
        }
    } else {  

        MPI_Send(&integral, 1, MPI_FLOAT, dest,
            tag, MPI_COMM_WORLD);
    }
    

    /* Print the result */
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n",
            n);
        printf("of the integral from %f to %f = %f\n",
            a, b, total);
    }

    /* Shut down MPI */
    MPI_Finalize();
} /*  main  */


